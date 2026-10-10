#!/usr/bin/env python3
"""三维查看 IBSimu-Cycl 导出的 VTK 文件（.vti 场图 / .vtp 轨迹）。

这是 ParaView 之外的**脚本化**入口，便于批量出图与写进教程。

依赖 `pyvista`（``pip install pyvista``）。若未安装，会自动退化为
`ibsimu_cycl.vtk_io` 的文本摘要——仍然能检查文件内容与数值范围。

用法::

    # 交互查看（弹出窗口）
    python3 examples/cyclotron/view_3d.py cycl_vtk_field.vti cycl_vtk_orbit.vtp

    # 批量出图（无显示环境可用）
    python3 examples/cyclotron/view_3d.py --save orbit.png cycl_vtk_orbit.vtp

    # 只看某个标量、指定切片法向
    python3 examples/cyclotron/view_3d.py --field Bmag --normal x field.vti

用 ParaView 查看（功能最全）::

    paraview cycl_vtk_field.vti cycl_vtk_orbit.vtp
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from ibsimu_cycl import vtk_io  # noqa: E402


def summarize(paths) -> None:
    """不依赖 pyvista 的文本摘要。"""
    for p in paths:
        try:
            d = vtk_io.read(str(p))
        except Exception as e:  # noqa: BLE001
            print(f"{p}: 读取失败：{e}")
            continue

        if isinstance(d, vtk_io.ImageData):
            print(f"{p}\n  ImageData  dims={d.dims}  origin={d.origin}  "
                  f"spacing={d.spacing}")
            for name, a in d.arrays.items():
                rng = (a.min(), a.max()) if a.size else (float("nan"),) * 2
                print(f"    {name:12s} shape={str(a.shape):12s} "
                      f"range=[{rng[0]:.6g}, {rng[1]:.6g}]")
        else:
            print(f"{p}\n  PolyData  points={len(d.points)}  lines={len(d.lines)}")
            if d.points.size:
                lo, hi = d.points.min(axis=0), d.points.max(axis=0)
                print(f"    x∈[{lo[0]:.6g}, {hi[0]:.6g}]  "
                      f"y∈[{lo[1]:.6g}, {hi[1]:.6g}]  z∈[{lo[2]:.6g}, {hi[2]:.6g}]")
            for name, a in d.arrays.items():
                print(f"    {name:12s} shape={a.shape}  "
                      f"range=[{a.min():.6g}, {a.max():.6g}]")

    print("\n（安装 pyvista 可获得三维交互视图：pip install pyvista）")


def view(paths, args) -> int:
    try:
        import pyvista as pv
    except ImportError:
        print("未安装 pyvista，退化为文本摘要。\n")
        summarize(paths)
        return(0)

    if args.save:
        pv.OFF_SCREEN = True

    pl = pv.Plotter(window_size=list(args.size))

    for p in paths:
        mesh = pv.read(str(p))

        if isinstance(mesh, pv.ImageData):
            names = list(mesh.point_data.keys())
            if not names:
                continue
            name = args.field or names[0]
            if name not in names:
                print(f"警告：{p} 中没有数组 '{name}'，可选：{names}")
                continue
            dims = mesh.dimensions
            center = [mesh.origin[i] + (dims[i] - 1) * mesh.spacing[i] / 2.0
                      for i in range(3)]
            sl = mesh.slice(normal=args.normal, origin=center)
            if sl.n_points == 0:                    # 切片退化时直接铺点
                sl = mesh
            sl.set_active_scalars(name)
            pl.add_mesh(sl, cmap=args.cmap, show_scalar_bar=True,
                        scalar_bar_args={"title": name})
            pl.add_mesh(mesh.outline(), color="black", line_width=2)
        else:
            scalar = args.field if (args.field in mesh.point_data) else (
                "t" if "t" in mesh.point_data else None)
            pl.add_mesh(mesh, line_width=args.line_width,
                        render_lines_as_tubes=True,
                        scalars=scalar, cmap=args.cmap,
                        show_scalar_bar=scalar is not None,
                        scalar_bar_args={"title": scalar or ""})

    pl.add_axes()
    if args.save:
        pl.screenshot(args.save)
        print(f"已保存图像：{args.save}")
    else:
        pl.show()
    return(0)


def main() -> int:
    ap = argparse.ArgumentParser(description="查看 IBSimu-Cycl 导出的 VTK 文件")
    ap.add_argument("files", nargs="+", help=".vti / .vtp 文件")
    ap.add_argument("--save", metavar="PNG", help="保存图像而不是弹窗（可无显示环境）")
    ap.add_argument("--field", help="要显示的数组名（默认第一个）")
    ap.add_argument("--normal", default="z", choices=["x", "y", "z"],
                    help="场图切片法向（默认 z）")
    ap.add_argument("--cmap", default="viridis", help="颜色表")
    ap.add_argument("--line-width", type=float, default=2.0, help="轨迹线宽")
    ap.add_argument("--size", type=int, nargs=2, default=[1024, 768],
                    help="窗口尺寸")
    args = ap.parse_args()
    return(view(args.files, args))


if __name__ == "__main__":
    sys.exit(main())
