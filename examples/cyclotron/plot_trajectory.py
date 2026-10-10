#!/usr/bin/env python3
"""可视化 IBSimu-Cycl 跟踪出的轨迹（cycl_track.csv）。

用法::

    python3 plot_trajectory.py [cycl_track.csv] [-o cyclotron_orbit.png]

生成 2x2 图：
  1) xy 平面轨迹（按 |B| 着色）
  2) |B| 随半径 r 的变化
  3) 速率 |v| 随时间（静磁场中应守恒）
  4) 速度方位角 phi_v 随时间（斜率即回旋角频率 omega_c）

无显示环境可用（使用 Agg 后端）。

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(path):
    # 文件含 1 行注释 + 1 行表头
    data = np.loadtxt(path, delimiter=",", comments="#", skiprows=2)
    cols = {name: data[:, i] for i, name in enumerate(
        ["t", "x", "y", "z", "vx", "vy", "vz", "Bmag"])}
    return cols


def unwrap_phase(vx, vy):
    phi = np.unwrap(np.arctan2(vy, vx))
    return phi


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default="cycl_track.csv")
    ap.add_argument("-o", "--output", default="cyclotron_orbit.png")
    args = ap.parse_args(argv)

    c = load(args.csv)
    t = c["t"]
    r = np.hypot(c["x"], c["y"])
    speed = np.sqrt(c["vx"]**2 + c["vy"]**2 + c["vz"]**2)
    phi = unwrap_phase(c["vx"], c["vy"])

    # 线性拟合回旋角频率 omega = dphi/dt
    A = np.vstack([t, np.ones_like(t)]).T
    slope, intercept = np.linalg.lstsq(A, phi, rcond=None)[0]
    omega = abs(slope)

    fig, ax = plt.subplots(2, 2, figsize=(11, 9))

    s0 = ax[0, 0].scatter(c["x"], c["y"], c=c["Bmag"], s=6, cmap="viridis")
    ax[0, 0].set_aspect("equal")
    ax[0, 0].set_xlabel("x [m]")
    ax[0, 0].set_ylabel("y [m]")
    ax[0, 0].set_title("Larmor gyration — NOT a cyclotron orbit (color = |B|)")
    # 说清楚这是“局部探针”，不是绕机器中心的轨道
    r_larmor = float(np.mean(speed))/omega
    ax[0, 0].annotate(
        "$r_L = v/\\omega_c = %.1f$ mm (diameter %.0f mm)\n"
        "0.1c proton in B = %.2f T.\n"
        "Deliberately local: the field is nearly\n"
        "uniform here, which is what lets the phase\n"
        "slope measure $\\omega_c=qB/m$ cleanly.\n"
        "(It still varies ~±25%% over this 0.4 m\n"
        "orbit, so the turns drift a little.)\n"
        "A real cyclotron orbit (r ≈ 3.3 m) needs\n"
        "a closed-orbit solve, not a tangential\n"
        "launch — see tests/cycl_closed_orbit."
        % (r_larmor*1e3, 2.0*r_larmor*1e3, float(np.mean(c["Bmag"]))),
        xy=(0.02, 0.02), xycoords="axes fraction", fontsize=7.5, va="bottom",
        bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.85))
    fig.colorbar(s0, ax=ax[0, 0], label="|B| [T]")

    ax[0, 1].plot(r, c["Bmag"], lw=1.5)
    ax[0, 1].set_xlabel("r [m]")
    ax[0, 1].set_ylabel("|B| [T]")
    ax[0, 1].set_title("Field magnitude vs radius")
    ax[0, 1].grid(True, alpha=0.3)

    ax[1, 0].plot(t*1e9, (speed/speed[0] - 1.0)*1e6, lw=1.0)
    ax[1, 0].set_xlabel("t [ns]")
    ax[1, 0].set_ylabel("relative speed drift [ppm]")
    ax[1, 0].set_title("Speed conservation (max |$\\Delta v$|/v = {:.1e})".format(
        np.abs(speed/speed[0] - 1.0).max()))
    ax[1, 0].grid(True, alpha=0.3)

    ax[1, 1].plot(t*1e9, phi, lw=1.5, label="unwrapped $\\varphi_v$")
    ax[1, 1].plot(t*1e9, slope*t + intercept, "--", lw=1.2, label="linear fit")
    ax[1, 1].set_xlabel("t [ns]")
    ax[1, 1].set_ylabel("$\\varphi_v$ [rad]")
    ax[1, 1].set_title(f"Velocity phase -> $\\omega_c$ = {omega:.4e} rad/s")
    ax[1, 1].legend(fontsize=9)
    ax[1, 1].grid(True, alpha=0.3)

    fig.suptitle("IBSimu-Cycl: Larmor gyration of a 0.1c proton — local probe "
                 "for $\\omega_c=qB/m$ (r = 3.3 m, real PSI Ring map)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.output, dpi=130)

    print(f"saved: {args.output}")
    print(f"points: {len(t)}  |B| range: {c['Bmag'].min():.4f}..{c['Bmag'].max():.4f} T")
    print(f"omega_c (from v-phase fit) = {omega:.6e} rad/s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
