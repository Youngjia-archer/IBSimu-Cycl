#!/usr/bin/env python3
"""可视化真实 PSI Ring 场强下的相对论回旋轨道（cycl_relativistic.csv）。

用法::

    ./tests/cycl_relativistic          # 写出 cycl_relativistic.csv
    python3 plot_relativistic_orbit.py [-o cyclotron_real_orbit.png]

列：t,x,y,z,gamma,r,Bmag

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="?", default="cycl_relativistic.csv")
    ap.add_argument("-o", "--output", default="cyclotron_real_orbit.png")
    args = ap.parse_args(argv)

    if not os.path.exists(args.csv):
        print(f"{args.csv} not found; run ./tests/cycl_relativistic first")
        return 1

    d = np.loadtxt(args.csv, delimiter=",", comments="#", skiprows=1)
    c = {k: d[:, i] for i, k in enumerate(["t", "x", "y", "z", "gamma", "r", "Bmag"])}

    fig, ax = plt.subplots(2, 2, figsize=(11, 9))

    s = ax[0, 0].scatter(c["x"], c["y"], c=c["Bmag"], s=4, cmap="viridis")
    ax[0, 0].set_aspect("equal")
    ax[0, 0].set_xlabel("x [m]")
    ax[0, 0].set_ylabel("y [m]")
    ax[0, 0].set_title("Real PSI Ring orbit (color = |B|)")
    fig.colorbar(s, ax=ax[0, 0], label="|B| [T]")

    ax[0, 1].plot(c["t"]*1e9, c["r"], lw=1.2)
    ax[0, 1].set_xlabel("t [ns]")
    ax[0, 1].set_ylabel("r [m]")
    ax[0, 1].set_title(f"Orbit radius (bounded, mean {c['r'].mean():.3f} m)")
    ax[0, 1].grid(True, alpha=0.3)

    ax[1, 0].plot(c["t"]*1e9, c["gamma"], lw=1.2)
    ax[1, 0].set_xlabel("t [ns]")
    ax[1, 0].set_ylabel(r"$\gamma$")
    ax[1, 0].set_title(r"$\gamma$ conserved (relativistic Boris, drift {:.1e})".format(
        np.abs(c["gamma"]/c["gamma"][0] - 1).max()))
    ax[1, 0].grid(True, alpha=0.3)

    ax[1, 1].plot(c["t"]*1e9, c["Bmag"], lw=0.8)
    ax[1, 1].set_xlabel("t [ns]")
    ax[1, 1].set_ylabel("|B| [T]")
    ax[1, 1].set_title("Field seen along the orbit (8-sector ripple)")
    ax[1, 1].grid(True, alpha=0.3)

    fig.suptitle("IBSimu-Cycl: relativistic proton in the real PSI Ring field", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.output, dpi=130)
    print(f"saved: {args.output}")
    print(f"points: {len(c['t'])}  <B>={c['Bmag'].mean():.4f} T  gamma={c['gamma'][0]:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
