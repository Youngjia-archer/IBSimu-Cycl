#!/usr/bin/env python3
"""解析圆柱腔（pillbox）TM010 模式生成器 —— RF 腔求解器的验证基准。

生成：
  - out/pillbox_E.txt : 电场图 (IBSIMU-CYCL-FIELDMAP)
  - out/pillbox_B.txt : 磁场图 (IBSIMU-CYCL-FIELDMAP)
  - out/pillbox_info.txt : 频率与几何信息

解析解（理想导体圆柱腔，半径 R）：
    k   = x01 / R              (x01 = J0 的第一个零点 = 2.4048255577)
    f   = c k / (2*pi)
    Ez(r)   = E0 * J0(k r)
    Bth(r)  = E0 / c * J1(k r)

用法::

    python3 pillbox_tm010.py --R 0.30 --E0 1.0e6 --nx 61 --ny 61 --nz 11 --outdir out

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import os
import sys
import numpy as np

C_LIGHT = 299792458.0
X01 = 2.4048255577          # J0 的第一个零点

MAGIC = "IBSIMU-CYCL-FIELDMAP"
VERSION = 1


# ---------------------------------------------------------------------------
# Bessel 函数（幂级数，机器精度）。
#   J0(x) = sum_m (-1)^m (x/2)^(2m) / (m!)^2
#   J1(x) = sum_m (-1)^m (x/2)^(2m+1) / (m! (m+1)!)
# 本问题自变量 k*r <= x01 ~ 2.40，级数收敛极快（数十项内到机器精度）。
# 使用精确级数而非 A&S 多项式近似，可保证 dJ0/dx = -J1 精确成立（Maxwell 校验需要）。
# ---------------------------------------------------------------------------
def j0(x):
    x = np.asarray(x, dtype=float)
    z = x*x/4.0
    term = np.ones_like(x)
    s = term.copy()
    for m in range(1, 120):
        term = term*(-z)/(m*m)
        s = s + term
        if np.all(np.abs(term) <= 1e-18*np.maximum(1.0, np.abs(s))):
            break
    return s


def j1(x):
    x = np.asarray(x, dtype=float)
    z = x*x/4.0
    term = x/2.0
    s = term.copy()
    for m in range(1, 120):
        term = term*(-z)/(m*(m+1))
        s = s + term
        if np.all(np.abs(term) <= 1e-18*np.maximum(1.0, np.abs(s))):
            break
    return s


def tm010_frequency(R):
    """TM010 谐振频率 [Hz]。"""
    return C_LIGHT*X01/(2.0*np.pi*R)


def build_modes(R, L, E0, nx, ny, nz):
    """返回 (x,y,z, E_grid, B_grid)，网格覆盖 [-R,R]x[-R,R]x[0,L]。"""
    x = np.linspace(-R, R, nx)
    y = np.linspace(-R, R, ny)
    z = np.linspace(0.0, L, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    r = np.hypot(X, Y)
    inside = r <= R
    k = X01/R

    Ez = np.where(inside, E0*j0(k*r), 0.0)
    Bt = np.where(inside, E0*j1(k*r)/C_LIGHT, 0.0)

    th = np.arctan2(Y, X)
    Bx = -np.sin(th)*Bt
    By = np.cos(th)*Bt

    E = np.stack([np.zeros_like(Ez), np.zeros_like(Ez), Ez], axis=-1)
    B = np.stack([Bx, By, np.zeros_like(Ez)], axis=-1)
    return x, y, z, E, B


def write_fieldmap(path, x, y, z, grid):
    nx, ny, nz = len(x), len(y), len(z)
    with open(path, "w") as f:
        f.write(f"{MAGIC} {VERSION}\n")
        f.write("# nx ny nz\n")
        f.write(f"{nx} {ny} {nz}\n")
        f.write("# x0 dx y0 dy z0 dz\n")
        f.write("{:.17g} {:.17g} {:.17g} {:.17g} {:.17g} {:.17g}\n".format(
            x[0], x[1]-x[0], y[0], y[1]-y[0], z[0], z[1]-z[0]))
        f.write("# data in order i(x), j(y), k(z): Fx Fy Fz\n")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    v = grid[i, j, k]
                    f.write("{:.9g} {:.9g} {:.9g}\n".format(v[0], v[1], v[2]))


def main(argv=None):
    ap = argparse.ArgumentParser(description="解析 pillbox TM010 模式生成器")
    ap.add_argument("--R", type=float, default=0.30, help="腔半径 [m]")
    ap.add_argument("--L", type=float, default=0.60, help="腔长度 [m]")
    ap.add_argument("--E0", type=float, default=1.0e6, help="轴上电场幅值 [V/m]")
    ap.add_argument("--nx", type=int, default=61)
    ap.add_argument("--ny", type=int, default=61)
    ap.add_argument("--nz", type=int, default=11)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)
    x, y, z, E, B = build_modes(args.R, args.L, args.E0, args.nx, args.ny, args.nz)

    write_fieldmap(os.path.join(args.outdir, "pillbox_E.txt"), x, y, z, E)
    write_fieldmap(os.path.join(args.outdir, "pillbox_B.txt"), x, y, z, B)

    f0 = tm010_frequency(args.R)
    with open(os.path.join(args.outdir, "pillbox_info.txt"), "w") as f:
        f.write(f"# analytic pillbox TM010\n")
        f.write(f"R = {args.R:.12g}\n")
        f.write(f"L = {args.L:.12g}\n")
        f.write(f"E0 = {args.E0:.12g}\n")
        f.write(f"x01 = {X01:.12g}\n")
        f.write(f"k = {X01/args.R:.12g}\n")
        f.write(f"f = {f0:.12g}\n")
        f.write(f"lambda = {C_LIGHT/f0:.12g}\n")

    print(f"R = {args.R} m  ->  f_TM010 = {f0/1e6:.6f} MHz  (lambda = {C_LIGHT/f0*100:.3f} cm)")
    print(f"grid: {args.nx}x{args.ny}x{args.nz}  ->  {args.outdir}/pillbox_{{E,B,info}}.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
