#!/usr/bin/env python3
"""可视化回旋加速器多圈加速（cycl_accel_*.csv）。

用法::

    python3 plot_acceleration.py [-o cyclotron_acceleration.png]

读取（存在则用）：
  cycl_accel_locked.csv  —— 均匀场 + 双 RF 间隙（相位锁定，11 圈）
  cycl_accel_real.csv    —— 真实 PSI Ring 场图（缩放）+ 闭合轨道起步 + RF 加速
  cycl_accel_iso.csv     —— 同上，但场按等时性要求做了径向修正 c(r)=w0/w(r)
列：turn,t,KE_MeV,r,phi_rad

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


def load(path):
    if not os.path.exists(path):
        return None
    d = np.loadtxt(path, delimiter=",", comments="#", skiprows=1)
    return {k: d[:, i] for i, k in enumerate(["turn", "t", "KE", "r", "phi"])}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default="cyclotron_acceleration.png")
    args = ap.parse_args(argv)

    locked = load("cycl_accel_locked.csv")
    real = load("cycl_accel_real.csv")
    iso = load("cycl_accel_iso.csv")
    if locked is None:
        print("cycl_accel_locked.csv not found; run ./tests/cycl_accel first")
        return 1

    fig, ax = plt.subplots(2, 2, figsize=(11, 9))

    ax[0, 0].plot(locked["turn"], locked["KE"], "o-", label="uniform B (locked)")
    if real is not None:
        ax[0, 0].plot(real["turn"], real["KE"], "s--", label="real map (as-is)")
    if iso is not None:
        ax[0, 0].plot(iso["turn"], iso["KE"], "^-",
                      label="real map + isochronous correction")
    ax[0, 0].set_xlabel("turn")
    ax[0, 0].set_ylabel("kinetic energy [MeV]")
    ax[0, 0].set_title("Energy gain per turn  (design $2qV_0$ = 40 keV/turn)")
    ax[0, 0].legend(fontsize=9)
    ax[0, 0].grid(True, alpha=0.3)

    # 两条曲线的绝对半径差 17 倍（均匀场 r0≈0.2 m，真实场 r0≈3.3 m），
    # 用 r/r0 才能在同一坐标下看清「半径随能量外扩」这一共同规律。
    r0 = locked["r"][0]
    ax[0, 1].plot(locked["turn"], locked["r"]/r0, "o-",
                  label="uniform B ($r_0$ = {:.3f} m)".format(r0))
    if real is not None:
        ax[0, 1].plot(real["turn"], real["r"]/real["r"][0], "s--",
                      label="real map as-is ($r_0$ = {:.3f} m)".format(real["r"][0]))
    if iso is not None:
        ax[0, 1].plot(iso["turn"], iso["r"]/iso["r"][0], "^-",
                      label="isochronous ($r_0$ = {:.3f} m)".format(iso["r"][0]))
    ax[0, 1].set_xlabel("turn")
    ax[0, 1].set_ylabel(r"orbit radius at gap,  $r/r_0$")
    ax[0, 1].set_title("Orbit radius grows with energy")
    ax[0, 1].legend(fontsize=9)
    ax[0, 1].grid(True, alpha=0.3)

    ax[1, 0].plot(locked["turn"], locked["phi"], "o-")
    ax[1, 0].set_xlabel("turn")
    ax[1, 0].set_ylabel(r"RF phase at crossing [rad]")
    ax[1, 0].set_ylim(-0.05, 0.05)
    ax[1, 0].text(0.03, 0.90,
                  "max $|\\varphi|$ = {:.1e} rad".format(np.abs(locked["phi"]).max()),
                  transform=ax[1, 0].transAxes, fontsize=9)
    ax[1, 0].set_title(r"Phase locked when $\omega_{RF}=\omega_c$")
    ax[1, 0].grid(True, alpha=0.3)

    if real is None:
        ax[1, 1].axis("off")
    else:
        ax[1, 1].plot(real["turn"], real["phi"], "s--", label="real map (as-is)")
        nt = real["turn"][-1] - real["turn"][0]
        slip = (real["phi"][-1] - real["phi"][0])/nt if nt > 0 else 0.0
        txt = "real:  {:+.3f} rad/turn".format(slip)
        if iso is not None:
            ax[1, 1].plot(iso["turn"], iso["phi"], "^-",
                          label="isochronous field")
            nti = iso["turn"][-1] - iso["turn"][0]
            slipi = (iso["phi"][-1] - iso["phi"][0])/nti if nti > 0 else 0.0
            txt += "\nisochronous: {:+.4f} rad/turn".format(slipi)
        ax[1, 1].set_xlabel("turn")
        ax[1, 1].set_ylabel(r"RF phase [rad]")
        ax[1, 1].text(0.03, 0.88, txt, transform=ax[1, 1].transAxes, fontsize=9,
                      bbox=dict(facecolor="white", edgecolor="0.7", alpha=0.9))
        ax[1, 1].set_title("Phase slip: locking needs a field with "
                           "$\\langle B_z\\rangle(r)\\propto\\gamma(r)$")
        ax[1, 1].legend(fontsize=9, loc="lower left")
        ax[1, 1].grid(True, alpha=0.3)

    fig.suptitle("IBSimu-Cycl: multi-turn cyclotron acceleration "
                 "(real PSI Ring field, closed-orbit launch)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(args.output, dpi=130)
    print(f"saved: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
