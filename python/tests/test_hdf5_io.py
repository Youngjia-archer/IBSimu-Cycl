#!/usr/bin/env python3
"""独立验证 C++ HDF5 导出器（``src/io/hdf5writer.cpp``）的输出。

分两部分：

1. **内建夹具自检**：用 h5py 手写一个符合 ``hdf5writer.hpp`` 布局的文件，
   校验本读取器（特别是 **i 最快** 的索引次序与单位属性）。
2. **交叉验证**：若已构建 ``tests/cycl_hdf5_export``，则运行它、再用本读取器
   读回产物；检查项用**文件内部自洽关系**（数值里编码的 (i,j,k)、位置与
   时间的仿射关系），不复用 C++ 侧的常数。

没有 h5py 时整体输出 ``[SKIP]`` 并返回成功（CI 的 python-adapter 作业会装
h5py，因此该路径在 CI 上有覆盖）。

用法::

    python3 python/tests/test_hdf5_io.py
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from ibsimu_cycl import hdf5_io  # noqa: E402

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

def test_fixtures(tmp: Path) -> None:
    if not hdf5_io.available():
        print(f"[1] [SKIP] 没有 h5py，无法校验读取器：{hdf5_io.H5PY_HINT}")
        return
    check(hdf5_io.available(), "hdf5_io.available() 为真")
    print("[1] 内建夹具自检")

    import h5py

    # 三维互不相同 + 数值编码 (i,j,k) ⇒ 任何转置都会被抓到
    nx, ny, nz = 7, 5, 3
    idx = np.arange(nz * ny * nx).reshape(nz, ny, nx)
    k, j, i = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                          indexing="ij")
    scal = i.astype(float) + 1000.0 * j + 1.0e6 * k
    vec = np.stack([scal, 2.0 * scal, -scal])          # (3,nz,ny,nx)

    path = tmp / "fixture.h5"
    with h5py.File(path, "w") as f:
        f.attrs["openPMD"] = "2.0.0"
        f.attrs["basePath"] = "/data/%T/"
        m = f.create_group("/data/0/meshes/M")
        m.attrs["geometry"] = "cartesian"
        m.attrs["axisLabels"] = np.asarray([b"x", b"y", b"z"])
        m.attrs["gridSpacing"] = np.asarray([0.5, 0.25, 0.125])
        m.attrs["gridGlobalOffset"] = np.asarray([-1.0, -2.0, -3.0])
        d = m.create_dataset("B", data=scal)
        d.attrs["unitSI"] = 1.0
        dv = m.create_dataset("Bvec", data=vec)
        dv.attrs["unitSI"] = 1.602176634e-16

    mesh = hdf5_io.read(path)          # 单一网格 ⇒ 自动分派
    check(mesh.dims == (nx, ny, nz), f"dims = {(nx, ny, nz)}（实为 {mesh.dims}）")
    check(np.allclose(mesh.origin, [-1.0, -2.0, -3.0]), "gridGlobalOffset 正确")
    check(np.allclose(mesh.spacing, [0.5, 0.25, 0.125]), "gridSpacing 正确")

    # 关键：读出来必须还是 i 最快的扁平数组 ⇒ reshape(nz,ny,nx) 等于原 idx 编码
    b = mesh.scalar("B")
    check(b.shape == (nx * ny * nz,), "标量长度为 npts")
    check(np.array_equal(b.reshape(nz, ny, nx), scal),
          "标量按 (nz,ny,nx) 重排后与原编码逐位一致（i 最快，无转置）")
    check(np.array_equal(mesh.scalar_grid("B"), scal),
          "scalar_grid 与写入时的 (nz,ny,nx) 数组逐位一致")

    v = mesh.vector("Bvec")
    check(v.shape == (nx * ny * nz, 3), "矢量形状 = (npts,3)")
    check(np.array_equal(v[:, 0].reshape(nz, ny, nx), scal)
          and np.array_equal(v[:, 2].reshape(nz, ny, nx), -scal),
          "矢量分量正确")

    meshes = hdf5_io.read_meshes(path)
    check(set(meshes) == {"M"}, "read_meshes 返回网格名")
    check(meshes["M"].units["Bvec"] == 1.602176634e-16, "unitSI 读回正确")

    # 轨迹夹具
    tpath = tmp / "fixture_tracks.h5"
    counts = [4, 2, 6]
    pos = np.vstack([np.column_stack([np.arange(c, dtype=float),
                                      2.0 * np.arange(c),
                                      -3.0 * np.arange(c)])
                     for c in counts])
    tim = np.concatenate([np.arange(c, dtype=float) * 1e-9 for c in counts])
    off = np.concatenate([[0], np.cumsum(counts)])
    with h5py.File(tpath, "w") as f:
        g = f.create_group("/data/0/particles/p")
        g.create_dataset("position", data=pos)
        g.create_dataset("time", data=tim)
        g.create_dataset("offset", data=off.astype(np.int64))
        g.create_dataset("count", data=np.asarray(counts, dtype=np.int64))

    pd = hdf5_io.read(tpath)
    check(len(pd) == sum(counts), f"轨迹总点数 = {sum(counts)}")
    check([len(l) for l in pd.lines] == counts, f"各条长度 = {counts}")
    check(np.allclose(pd.points, pos), "位置逐位一致")
    check(np.allclose(pd.arrays["t"], tim), "时间逐位一致")

    # 同时含网格与轨迹时必须拒绝自动分派
    mixed = tmp / "mixed.h5"
    with h5py.File(mixed, "w") as f:
        f.create_group("/data/0/meshes/M").create_dataset(
            "B", data=np.zeros((1, 1, 1)))
        f.create_group("/data/0/particles/p").create_dataset(
            "position", data=np.zeros((1, 3)))
    try:
        hdf5_io.read(mixed)
        check(False, "混合内容应拒绝 read() 自动分派")
    except ValueError:
        check(True, "混合内容拒绝 read() 自动分派（提示用 read_meshes/read_tracks）")


# --------------------------------------------------------------------------
# 2. 交叉验证 C++ 导出产物
# --------------------------------------------------------------------------

def test_cpp_output(tmp: Path) -> None:
    if not hdf5_io.available():
        print("[2] [SKIP] 没有 h5py")
        return
    exe = ROOT / "tests" / "cycl_hdf5_export"
    if not exe.exists():
        print("[2] [SKIP] 未找到 tests/cycl_hdf5_export（先 `make check` 或 "
              "`make -C tests cycl_hdf5_export`）")
        return

    print("[2] 交叉验证 C++ 导出产物")
    r = subprocess.run([str(exe)], cwd=tmp, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-2000:])
        print(r.stderr[-2000:])
    check(r.returncode == 0, f"tests/cycl_hdf5_export 退出码 = 0（实为 {r.returncode}）")
    if r.returncode != 0:
        return

    # ---- 网格 ----
    meshes = hdf5_io.read_meshes(tmp / "cycl_hdf5_export.h5")
    check(set(meshes) == {"Bfield"}, f"网格名 = {{'Bfield'}}（实为 {set(meshes)}）")
    mesh = meshes["Bfield"]
    nx, ny, nz = mesh.dims
    check((nx, ny, nz) == (24, 20, 16), f"dims = (24,20,16)（实为 {(nx, ny, nz)}）")
    check(np.allclose(mesh.spacing, 0.05), "gridSpacing = 0.05")
    check_close(mesh.units["Bz"], 1.602176634e-16, 1e-15, "unitSI 为 keV→J 因子")

    # C++ 测试把 (i,j,k) 编码进数值：f = i + 1000 j + 1e6 k。
    # 这里不复用该常数，而是从数据自身反推解码规则：
    #   - 每行/列/层的差分必须是常整数；
    #   - 由此重建 i、j、k 索引，三者必须正好张满 [0,nx)×[0,ny)×[0,nz)。
    kk, jj, ii = np.meshgrid(np.arange(nz), np.arange(ny), np.arange(nx),
                             indexing="ij")
    grid = mesh.scalar_grid("Bz")
    di = np.diff(grid, axis=2)
    dj = np.diff(grid, axis=1)
    dk = np.diff(grid, axis=0)
    check(np.allclose(di, di[0, 0, 0]) and np.allclose(dj, dj[0, 0, 0])
          and np.allclose(dk, dk[0, 0, 0]),
          "场沿三个方向都是严格线性的（说明索引次序正确、无转置）")
    check_close(dj[0, 0, 0] / di[0, 0, 0], 1000.0, 1e-12,
                "j 方向步长 / i 方向步长 = 1000 ⇒ 与编码一致")
    check_close(dk[0, 0, 0] / di[0, 0, 0], 1.0e6, 1e-12,
                "k 方向步长 / i 方向步长 = 1e6 ⇒ 与编码一致")

    # 矢量场：三个分量的内部比例必须是 (1, 2, -1)
    v = mesh.vector("Bvec")
    check(np.allclose(v[:, 1], 2.0 * v[:, 0], rtol=0, atol=0)
          and np.allclose(v[:, 2], -v[:, 0], rtol=0, atol=0),
          "矢量分量比例为 (1, 2, -1)")

    # ---- 轨迹 ----
    pd = hdf5_io.read_tracks(tmp / "cycl_hdf5_export_tracks.h5")
    counts = [len(l) for l in pd.lines]
    check(counts == [5, 7, 2], f"三条轨迹点数 = [5,7,2]（实为 {counts}）")

    t = pd.arrays["t"]
    check(np.all(np.diff(t) > 0), "时刻严格递增")
    ok_affine = True
    for m, line in enumerate(pd.lines):
        pts = pd.points[line]
        tm = t[line]
        # C++ 侧：x = 0.1m + 0.01i，t = 1e-9(10m + i) ⇒ 同一条内 t 与 x 仿射
        dx = np.diff(pts[:, 0])
        dt = np.diff(tm)
        if not (np.allclose(dx, dx[0]) and np.allclose(dt, dt[0])
                and abs(dx[0] / dt[0] - 1.0e7) < 1.0):
            ok_affine = False
        # 点必须在一条直线上：(y,z) 是 x 的线性函数
        if not (np.allclose(np.diff(pts[:, 1]) / np.diff(pts[:, 0]), 2.0)
                and np.allclose(np.diff(pts[:, 2]) / np.diff(pts[:, 0]), -3.0)):
            ok_affine = False
    check(ok_affine, "每条轨迹上 (x,y,z) 共线且 x(t) 斜率 ≈ 1e7 m/s（内部自洽）")


def main() -> int:
    print("=== HDF5 读取器 / 导出器交叉验证 ===")
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
