#!/usr/bin/env python3
"""绘制真实 PSI Ring 场强下的**闭合轨道**（单圈映射的不动点）。

数据来源::

    ./tests/cycl_closed_orbit     # -> cycl_closed_orbit.vtp（3 圈，带时刻）
    ./tests/cycl_track            # -> cycl_track_map.vti（整机中场图，可选底图）

用法::

    python3 examples/cyclotron/plot_closed_orbit.py \
        --orbit tests/cycl_closed_orbit.vtp \
        --map   tests/cycl_track_map.vti \
        -o      docs/img/cyclotron_real_orbit.png

> **这是单能量稳态轨道，不含 RF 加速。** 真实回旋加速器的束流每圈获得能量、
> 轨道半径随之外扩（$r=\\sqrt{2KE/m}/\\omega_c$），是一条向外盘旋的螺旋；
> 本图画的是**固定能量下**的闭合轨道——加速暂未纳入（见 ROADMAP）。

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
    ap.add_argument("--orbit", default="cycl_closed_orbit.vtp")
    ap.add_argument("--map", default="cycl_track_map.vti",
                    help="整机中场图（可选底图）")
    ap.add_argument("-o", "--output", default="cyclotron_real_orbit.png")
    args = ap.parse_args(argv)

    if not os.path.exists(args.orbit):
        print(f"{args.orbit} not found; run ./tests/cycl_closed_orbit first")
        return 1

    pd = read(args.orbit)
    x, y = pd.points[:, 0], pd.points[:, 1]
    t = pd.arrays["t"]
    r = np.hypot(x, y)

    fd = read(args.map) if os.path.exists(args.map) else None

    fig, ax = plt.subplots(1, 2, figsize=(13, 6))

    # --- 左：闭合轨道（xy，叠加中场图作底图）---
    a = ax[0]
    if fd is not None:
        bmag = fd.scalar_grid("Bmag")[0]
        x0, y0 = fd.origin[0], fd.origin[1]
        h = fd.spacing[0]
        x1 = x0 + h*(fd.dims[0]-1)
        y1 = y0 + h*(fd.dims[1]-1)
        im = a.imshow(bmag, origin="lower", extent=[x0, x1, y0, y1],
                      cmap="Blues_r", interpolation="bilinear", alpha=0.55)
        fig.colorbar(im, ax=a, label="|B| [T]", shrink=0.85)
    a.plot(x, y, color="#d62728", lw=2.0,
           label="closed orbit (fixed point of the sector map)")
    a.plot(x[0], y[0], "o", color="white", markeredgecolor="#d62728", ms=9,
           zorder=5, label="start point")
    a.set_aspect("equal")
    a.set_xlabel("x [m]")
    a.set_ylabel("y [m]")
    a.set_title("Closed orbit at fixed energy ($\\gamma$=1.227)\n"
                "scalloped by the 8 separated sectors")
    a.legend(loc="upper right", fontsize=9)
    a.grid(alpha=0.20)

    # --- 右：半径随时间（8 折扇贝调制）---
    a = ax[1]
    a.plot(t*1e9, r*1e3, lw=1.4, color="#d62728")
    a.axhline(r.mean()*1e3, color="0.4", ls="--", lw=1,
              label=f"mean r = {r.mean():.4f} m")
    a.set_xlabel("t [ns]")
    a.set_ylabel("r [mm]")
    a.set_title("Orbit radius: only the physical scalloping\n"
                f"span = {1000*(r.max()-r.min()):.0f} mm "
                f"({100*(r.max()-r.min())/r.mean():.1f}%)")
    a.grid(alpha=0.3)
    a.legend(fontsize=9)

    fig.suptitle("IBSimu-Cycl: closed orbit of a 210 MeV proton in the real "
                 "PSI Ring field (single energy, no RF acceleration)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"saved: {out}")

    print(f"points: {len(t)}  duration: {t[-1]*1e9:.1f} ns")
    print(f"r = {r.min():.4f} .. {r.max():.4f} m  "
          f"(span {1e3*(r.max()-r.min()):.1f} mm, "
          f"{100*(r.max()-r.min())/r.mean():.2f}%)")
    print("note: single-energy closed orbit; RF acceleration (radius growth "
          "per turn) is NOT included")
    return 0


if __name__ == "__main__":
    sys.exit(main())
