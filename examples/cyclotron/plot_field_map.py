#!/usr/bin/env python3
# NOTE: B_REF_DIRECT is the value tests/cycl_track prints from evaluating the
# field map directly (not through the 0.1 m export grid).
"""绘制 PSI Ring **整机中场图**（**不含任何粒子轨迹**）。

数据来源是 ``tests/cycl_track`` 导出的标准 VTK XML（见 ``src/io/vtkwriter``）::

    make -C tests cycl_track && ./tests/cycl_track

    tests/cycl_track_map.vti    整机中场图（中平面 360° 全周，r∈[1.9,4.7] m）

用法::

    python3 examples/cyclotron/plot_field_map.py \
        --map tests/cycl_track_map.vti -o docs/img/cyclotron_field_map.png

生成 2x2 图：
  1) 整机中场 |B|（含参考点标记）
  2) 方位平均场 <Bz>(r)：决定轨道半径的量
  3) |B| 沿方位角 —— 8 折扇形结构
  4) |B| 沿半径（扇区中心 vs 扇区边界）

> 本图**不画粒子轨迹**。过去在这里叠过一条 0.2 m 半径的“局部拉莫尔回旋”，
> 那不是回旋加速器轨道，已删除；机器尺度的闭合轨道见
> ``plot_closed_orbit.py`` 与 ``tests/cycl_closed_orbit.cpp``。

无显示环境可用（Agg 后端）。

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from ibsimu_cycl.vtk_io import read  # noqa: E402

B_REF_DIRECT = 0.6693   # <Bz>(r=3.3 m)，tests/cycl_track 直接在场图上求值


def sample_nearest(grid, x0, y0, h, xs, ys):
    """在规则网格上做最近邻取样（够画曲线用，避免引入 scipy 依赖）。"""
    nx = grid.shape[1]
    i = np.clip(np.rint((xs - x0) / h).astype(int), 0, nx - 1)
    j = np.clip(np.rint((ys - y0) / h).astype(int), 0, grid.shape[0] - 1)
    return grid[j, i]


def main(argv=None):
    ap = argparse.ArgumentParser(description="绘制 PSI Ring 场图与轨道图")
    ap.add_argument("--map", default="tests/cycl_track_map.vti",
                    help="整机中场图 .vti（由 cycl_track 导出）")
    ap.add_argument("-o", "--output", default="docs/img/cyclotron_field_map.png")
    args = ap.parse_args(argv)

    # ---- 场图 ----
    fd = read(args.map)
    if fd.dims[2] != 1:
        raise SystemExit(f"{args.map}: 期望单层中平面切片，实为 {fd.dims}")
    bmag = fd.scalar_grid("Bmag")[0]              # (ny, nx)
    bz = fd.vector("B")[:, 2].reshape(fd.dims[2], fd.dims[1],
                                       fd.dims[0])[0]   # 真正的 Bz 分量
    x0, y0 = fd.origin[0], fd.origin[1]
    h = fd.spacing[0]
    x1 = x0 + h * (fd.dims[0] - 1)
    y1 = y0 + h * (fd.dims[1] - 1)

    xa = np.linspace(x0, x1, fd.dims[0])
    ya = np.linspace(y0, y1, fd.dims[1])

    # ---- 参考值自检（与 cycl_track 的解析判据同一量）----
    r_ref, th_ref = 3.3, np.pi / 2
    b_ref = float(sample_nearest(bmag, x0, y0, h,
                                 np.array([r_ref * np.cos(th_ref)]),
                                 np.array([r_ref * np.sin(th_ref)]))[0])

    fig, ax = plt.subplots(2, 2, figsize=(13, 11))

    # --- (1) 整机中场图 ---
    a = ax[0][0]
    im = a.imshow(bmag, origin="lower", extent=[x0, x1, y0, y1],
                  cmap="viridis", interpolation="bilinear")
    a.plot([r_ref * np.cos(th_ref)], [r_ref * np.sin(th_ref)], "o",
           color="#ff5252", ms=7, label=f"reference point (r={r_ref} m)")
    a.set_aspect("equal")
    a.set_xlabel("x [m]")
    a.set_ylabel("y [m]")
    a.set_title("PSI Ring midplane |B| (full 360$^\\circ$")
    a.legend(loc="upper right", fontsize=9)
    fig.colorbar(im, ax=a, label="|B| [T]", shrink=0.85)

    # --- (2) 方位平均场 <Bz>(r)：回旋加速器轨道的决定量 ---
    a = ax[0][1]
    rr = np.linspace(1.95, 4.65, 260)
    phi = np.linspace(0, 2*np.pi, 721)
    bavg = np.array([
        sample_nearest(bz, x0, y0, h,
                       r*np.cos(phi), r*np.sin(phi)).mean()
        for r in rr])
    bmax = np.array([
        sample_nearest(bmag, x0, y0, h,
                       r*np.cos(phi), r*np.sin(phi)).max()
        for r in rr])
    a.plot(rr, bavg, lw=1.8, label=r"$\langle B_z\rangle$ (azimuthal mean)")
    a.plot(rr, bmax, lw=1.2, ls="--", color="#d62728",
           label="peak |B| (pole)")
    # 直接在场图上求值（tests/cycl_track 的解析判据）作为参照：
    # 本曲线取自 0.1 m 导出网格的最近邻采样，会把近乎为零的谷区抹平，
    # 因此平均值偏高（~0.78 vs 0.67 T）。
    a.plot([r_ref], [B_REF_DIRECT], marker="*", ms=14, color="#ff9800",
           zorder=6, label=f"direct evaluation: {B_REF_DIRECT} T")
    a.annotate("curve uses nearest-neighbour sampling on the\n"
               "0.1 m export grid, which smears the\n"
               "near-zero valleys -> biased high",
               xy=(0.55, 0.20), xycoords="axes fraction", fontsize=8,
               bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.85))
    a.axhline(bavg[np.argmin(np.abs(rr - r_ref))], color="0.4", lw=0.8, ls=":")
    a.text(2.05, bavg[np.argmin(np.abs(rr - r_ref))]*0.32,
           "pole fill factor = $\\langle B_z\\rangle/B_{pole}$ ≈ 0.43\n"
           "(poles cover only ~19$^\\circ$ of each 45$^\\circ$ sector)",
           fontsize=8)
    a.set_xlabel("radius r [m]")
    a.set_ylabel("B [T]")
    a.set_title(r"$\langle B_z\rangle(r)$: what sets the orbit radius")
    a.grid(alpha=0.3)
    a.legend(fontsize=9)

    # --- (3) |B| 沿方位角（8 折扇形结构）---
    a = ax[1][0]
    phi = np.linspace(0, 2 * np.pi, 721)
    for rr, c in ((3.10, "#1f77b4"), (3.30, "#2ca02c"), (3.50, "#d62728")):
        b = sample_nearest(bmag, x0, y0, h, rr * np.cos(phi), rr * np.sin(phi))
        a.plot(np.degrees(phi), b, lw=1.4, color=c, label=f"r = {rr} m")
    a.set_xlabel("azimuth $\\theta$ [deg]")
    a.set_ylabel("|B| [T]")
    a.set_xlim(0, 360)
    a.set_xticks(np.arange(0, 361, 45))
    a.set_title("|B| vs azimuth: 8-fold sector structure (45$^\\circ$ period)")
    a.grid(alpha=0.3)
    a.legend(fontsize=9)

    # --- (4) |B| 沿半径：扇区中心 vs 扇区边界 ---
    a = ax[1][1]
    rr = np.linspace(1.95, 4.65, 400)
    for deg, c, lbl in ((0, "#1f77b4", "sector centre ($\\theta$=0$^\\circ$)"),
                        (22.5, "#d62728", "sector edge ($\\theta$=22.5$^\\circ$)")):
        b = sample_nearest(bmag, x0, y0, h,
                           rr * np.cos(np.radians(deg)),
                           rr * np.sin(np.radians(deg)))
        a.plot(rr, b, lw=1.6, color=c, label=lbl)
    a.axvline(r_ref, color="k", ls="--", lw=1, label=f"reference radius {r_ref} m")
    a.set_xlabel("radius r [m]")
    a.set_ylabel("|B| [T]")
    a.set_title("|B| vs radius: average field rise + sector modulation")
    a.grid(alpha=0.3)
    a.legend(fontsize=9)

    fig.suptitle("IBSimu-Cycl: PSI Ring field map and five-turn orbit "
                 "(real OPAL/PSI field data)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"已保存：{out}")

    # ---- 控制台自检 ----
    print(f"场图: {fd.dims[0]}x{fd.dims[1]} 节点, x∈[{x0:.3f},{x1:.3f}] m, "
          f"步长 {h} m")
    print(f"|B| 范围: {bmag.min():.4f} ~ {bmag.max():.4f} T")
    print(f"⟨Bz⟩(r={r_ref} m) = "
          f"{bavg[np.argmin(np.abs(rr - r_ref))]:.4f} T")
    print(f"参考点 r={r_ref} m, θ=90°: |B| = {b_ref:.6f} T "
          f"（tests/cycl_track 解析判据为 1.557581 T）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
