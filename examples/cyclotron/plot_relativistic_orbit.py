#!/usr/bin/env python3
"""可视化真实 PSI Ring 场强下的相对论回旋轨道（含闭合轨道对比）。

数据来源::

    ./tests/cycl_relativistic     # -> cycl_relativistic.csv（朴素发射：不闭合）
    ./tests/cycl_closed_orbit     # -> cycl_closed_orbit.vtp（闭合轨道＝不动点）

用法::

    python3 plot_relativistic_orbit.py \\
        --naive cycl_relativistic.csv --closed cycl_closed_orbit.vtp \\
        -o docs/img/cyclotron_real_orbit.png

图的核心信息：**朴素发射（r=3.3 m 纯切向）不是闭合轨道**——径向摆动 0.62 m，
相邻两圈在同一方位角上的半径差达 0.24 m，每圈都走一条新曲线（看起来“轨道重叠”）；
而数值求出的闭合轨道（单圈映射的不动点）每圈重复，摆动只有真实的 8 折扇贝调制。

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python"))

from ibsimu_cycl.vtk_io import read  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--naive", default="cycl_relativistic.csv",
                    help="朴素发射轨迹 CSV（tests/cycl_relativistic 输出）")
    ap.add_argument("--closed", default="cycl_closed_orbit.vtp",
                    help="闭合轨道 VTP（tests/cycl_closed_orbit 输出）")
    ap.add_argument("-o", "--output", default="cyclotron_real_orbit.png")
    args = ap.parse_args(argv)

    if not os.path.exists(args.naive):
        print(f"{args.naive} not found; run ./tests/cycl_relativistic first")
        return 1

    d = np.loadtxt(args.naive, delimiter=",", comments="#", skiprows=1)
    c = {k: d[:, i] for i, k in
         enumerate(["t", "x", "y", "z", "gamma", "r", "Bmag"])}

    closed = None
    if os.path.exists(args.closed):
        pd = read(args.closed)
        closed = { "x": pd.points[:, 0], "y": pd.points[:, 1],
                   "t": pd.arrays["t"] }
    else:
        print(f"warning: {args.closed} not found; run ./tests/cycl_closed_orbit "
              "for the closed-orbit comparison")

    fig, ax = plt.subplots(2, 2, figsize=(12, 10))

    # --- (0,0) 两种轨道的直接对比 ---
    a = ax[0, 0]
    sc = a.scatter(c["x"], c["y"], c=c["Bmag"], s=4, cmap="viridis",
                   label="unmatched launch (turns do not repeat)")
    if closed is not None:
        a.plot(closed["x"], closed["y"], color="#ff3b30", lw=2.2,
               label="closed orbit = fixed point of the sector map")
        a.plot(closed["x"][0], closed["y"][0], "o", color="white",
               markeredgecolor="#ff3b30", ms=9, zorder=5)
    a.set_aspect("equal")
    a.set_xlabel("x [m]")
    a.set_ylabel("y [m]")
    a.set_title("Real PSI Ring field: unmatched launch vs closed orbit\n"
                "(marker = start point at r = 3.30 m)")
    a.legend(loc="upper right", fontsize=8)
    fig.colorbar(sc, ax=a, label="|B| [T]", shrink=0.85)

    # --- (0,1) 半径随时间：直观显示“逐圈不重合” ---
    a = ax[0, 1]
    a.plot(c["t"]*1e9, c["r"], lw=1.0, color="#1f77b4",
           label=f"unmatched launch (span {c['r'].max()-c['r'].min():.2f} m)")
    if closed is not None:
        rr = np.hypot(closed["x"], closed["y"])
        a.plot(closed["t"]*1e9, rr, lw=1.6, color="#ff3b30",
               label=f"closed orbit (span {rr.max()-rr.min():.2f} m)")
    a.set_xlabel("t [ns]")
    a.set_ylabel("r [m]")
    a.set_title("Orbit radius: the unmatched launch wanders over 0.62 m")
    a.grid(True, alpha=0.3)
    a.legend(fontsize=8)

    # --- (1,0) gamma 守恒 ---
    a = ax[1, 0]
    a.plot(c["t"]*1e9, c["gamma"], lw=1.2)
    a.set_xlabel("t [ns]")
    a.set_ylabel(r"$\gamma$")
    a.set_title(r"$\gamma$ conserved (relativistic Boris, drift {:.1e})".format(
        np.abs(c["gamma"]/c["gamma"][0] - 1).max()))
    a.grid(True, alpha=0.3)

    # --- (1,1) 沿轨道的场强 ---
    a = ax[1, 1]
    a.plot(c["t"]*1e9, c["Bmag"], lw=0.8)
    a.set_xlabel("t [ns]")
    a.set_ylabel("|B| [T]")
    a.set_title("Field along the orbit: 8 separated sectors\n"
                "(poles cover only ~19 deg of each 45 deg)")
    a.grid(True, alpha=0.3)

    fig.suptitle("IBSimu-Cycl: relativistic proton in the real PSI Ring field",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"saved: {out}")
    print(f"unmatched : r = {c['r'].min():.4f} .. {c['r'].max():.4f} m "
          f"(span {c['r'].max()-c['r'].min():.4f} m), "
          f"<|B|>={c['Bmag'].mean():.4f} T, gamma={c['gamma'][0]:.4f}")
    if closed is not None:
        rr = np.hypot(closed["x"], closed["y"])
        print(f"closed    : r = {rr.min():.4f} .. {rr.max():.4f} m "
              f"(span {rr.max()-rr.min():.4f} m)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
