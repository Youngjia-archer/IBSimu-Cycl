#!/usr/bin/env python3
"""独立验证 C++ VTK 导出器（``src/io/vtkwriter.cpp``）的输出。

分两部分：

1. **内建夹具自检**：手写的 .vti/.vtp 片段 → 校验读取器本身（无外部依赖，
   任何环境都能跑）。
2. **交叉验证**：若已构建 ``tests/cycl_vtk_export``，则运行它、再用本读取器
   读回产物并做**物理/内部一致性**检查。检查项刻意不复用 C++ 侧的常数，
   而是利用文件内部的自洽关系（例如 |B| 与 B 分量必须一致、z 与 t 成正比、
   轨道半径恒定），因此两侧任一出错都会被发现。

用法::

    python3 python/tests/test_vtk_io.py        # 仓库根目录或任意位置均可
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from ibsimu_cycl.vtk_io import ImageData, PolyData, read  # noqa: E402

FAILED: list[str] = []


def check(ok: bool, what: str) -> None:
    print(f"  [{'OK  ' if ok else 'FAIL'}] {what}")
    if not ok:
        FAILED.append(what)


def check_close(a, b, tol: float, what: str) -> None:
    ok = abs(float(a) - float(b)) <= tol * max(abs(float(a)), abs(float(b)), 1e-300)
    if not ok:
        print(f"       got {a!r}, expected {b!r} (tol {tol:g} rel)")
    check(ok, what)


# --------------------------------------------------------------------------
# 1. 内建夹具：校验读取器
# --------------------------------------------------------------------------

FIXTURE_VTI = """<?xml version="1.0"?>
<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">
  <ImageData WholeExtent="0 1 0 0 0 0" Origin="0 0 0" Spacing="0.5 0.5 0.5">
    <Piece Extent="0 1 0 0 0 0">
      <PointData Scalars="s">
        <DataArray type="Float64" Name="s" NumberOfComponents="1" format="ascii">
          1 2
        </DataArray>
      </PointData>
    </Piece>
  </ImageData>
</VTKFile>
"""

FIXTURE_VTP = """<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">
  <PolyData>
    <Piece NumberOfPoints="4" NumberOfVerts="0" NumberOfLines="2"
           NumberOfStrips="0" NumberOfPolys="0">
      <PointData Scalars="t">
        <DataArray type="Float64" Name="t" NumberOfComponents="1" format="ascii">
          0 1 2 3
        </DataArray>
      </PointData>
      <Points>
        <DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">
          0 0 0  1 0 0  5 5 5  6 5 5
        </DataArray>
      </Points>
      <Lines>
        <DataArray type="Int64" Name="connectivity" format="ascii">
          0 1 2 3
        </DataArray>
        <DataArray type="Int64" Name="offsets" format="ascii">
          2 4
        </DataArray>
      </Lines>
    </Piece>
  </PolyData>
</VTKFile>
"""


def test_fixtures(tmp: Path) -> None:
    print("\n[1] 读取器内建夹具自检")

    p = tmp / "fixture.vti"
    p.write_text(FIXTURE_VTI)
    fd = read(str(p))
    check(isinstance(fd, ImageData), "vti → ImageData")
    check(fd.dims == (2, 1, 1), "dims = (2,1,1)")
    check_close(fd.spacing[0], 0.5, 1e-12, "spacing = 0.5")
    check(np.allclose(fd.scalar("s"), [1.0, 2.0]), "标量值 = [1, 2]")
    check(fd.scalar_grid("s").shape == (1, 1, 2), "scalar_grid 形状 (nz,ny,nx)")

    p = tmp / "fixture.vtp"
    p.write_text(FIXTURE_VTP)
    pd = read(str(p))
    check(isinstance(pd, PolyData), "vtp → PolyData")
    check(pd.points.shape == (4, 3), "points 形状 (4,3)")
    check(len(pd.lines) == 2, "两条折线")
    check(list(pd.lines[0]) == [0, 1] and list(pd.lines[1]) == [2, 3],
          "connectivity/offsets 切分正确")
    check(np.allclose(pd.arrays["t"], [0, 1, 2, 3]), "点数据 t 正确")


# --------------------------------------------------------------------------
# 2. 交叉验证 C++ 导出产物
# --------------------------------------------------------------------------

def test_cpp_output(tmp: Path) -> None:
    exe = ROOT / "tests" / "cycl_vtk_export"
    if not exe.exists():
        print("\n[2] [SKIP] 未找到 tests/cycl_vtk_export（先 `make check` 或 "
              "`make -C tests cycl_vtk_export`）")
        return

    print("\n[2] 交叉验证 C++ 导出产物")
    r = subprocess.run([str(exe)], cwd=tmp, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-2000:])
    check(r.returncode == 0, f"tests/cycl_vtk_export 退出码 = 0（实为 {r.returncode}）")
    if r.returncode != 0:
        return

    # ---- 场图 .vti ----
    fd = read(str(tmp / "cycl_vtk_field.vti"))
    nx, ny, nz = fd.dims
    check((nx, ny, nz) == (9, 7, 5), "场图 dims = (9,7,5)")

    b = fd.vector("B")
    bmag = fd.scalar("Bmag")
    check(b.shape == (nx * ny * nz, 3), "矢量 B 形状 = (npts,3)")
    check(bmag.shape == (nx * ny * nz,), "标量 Bmag 长度 = npts")

    # 内部自洽：Bmag 必须等于 |B|（写出器把两者都写进同一文件）
    check(np.allclose(bmag, np.linalg.norm(b, axis=1), rtol=1e-12, atol=0),
          "Bmag == |B|（内部自洽）")
    # 该算例磁场只有 z 分量
    check(np.allclose(b[:, 0], 0.0) and np.allclose(b[:, 1], 0.0),
          "Bx == By == 0（仅轴向分量）")

    # 坐标轴与 origin/spacing 一致；并用场自身的对称性验证 i 最快排列：
    #   - 该算例场与 z 无关 → 各 k 层必须逐点相同；
    #   - cos 型分布在 x∈[-0.04,0.04] 上关于 x=0 对称 → 沿 i 轴必须镜像对称。
    # 若数组被转置（例如按 k 最快），这两条都会失败。
    x = fd.axis(0)
    bz_grid = b[:, 2].reshape(nz, ny, nx)
    check(np.all(np.diff(x) > 0), "x 轴单调递增")
    check(np.allclose(bz_grid, bz_grid[0][None, :, :], rtol=0, atol=0),
          "场与 z 无关 ⇒ 各 k 层逐点相同")
    check(np.allclose(bz_grid, bz_grid[:, :, ::-1], rtol=0, atol=0),
          "场关于 x=0 对称 ⇒ i 最快排列正确（无转置）")

    # ---- 轨道 .vtp ----
    pd = read(str(tmp / "cycl_vtk_orbit.vtp"))
    check(len(pd) == 801, f"轨道点数 = 801（实为 {len(pd)}）")
    check(len(pd.lines) == 1 and len(pd.lines[0]) == 801, "单条折线含全部 801 点")

    t = pd.arrays["t"]
    check(np.all(np.diff(t) > 0), "时刻严格递增")
    z = pd.points[:, 2]
    # 无 Bz 的 z 方向受力 → z 与 t 严格成正比（比值恒定）
    ratio = z[1:] / t[1:]
    check_close(ratio.min(), ratio.max(), 1e-9, "z/t 恒定 ⇒ z = v_parallel·t")

    # 已初始化的轨道半径应几乎恒定（相对变化 < 1e-3）
    r = np.hypot(pd.points[:, 0], pd.points[:, 1])
    span = (r.max() - r.min()) / r.mean()
    check(span < 1e-3, f"回旋半径恒定（相对振幅 {span:.2e} < 1e-3）")

    # 首末点：追了两整圈 ⇒ 横向位置应回到起点
    dx = pd.points[-1, 0] - pd.points[0, 0]
    dy = pd.points[-1, 1] - pd.points[0, 1]
    check(np.hypot(dx, dy) / r.mean() < 1e-3,
          f"两整圈后回到起点（{np.hypot(dx, dy) / r.mean():.2e} < 1e-3）")


def main() -> int:
    print("=== VTK 读取器 / 导出器交叉验证 ===")
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        test_fixtures(tmp)
        test_cpp_output(tmp)

    if FAILED:
        print(f"\n结果：*** {len(FAILED)} 项未通过 ***")
        for f in FAILED:
            print(f"  - {f}")
        return 1
    print("\n结果：全部通过")
    return 0


if __name__ == "__main__":
    sys.exit(main())
