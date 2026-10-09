#!/usr/bin/env python3
"""解析 pillbox TM010 的自检（Maxwell 与频率）。

用法: python3 tests/test_pillbox_tm010.py

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # adapters/rfcavity3d

import pillbox_tm010 as pb                          # noqa: E402

FAILS = 0


def check(ok, msg):
    global FAILS
    print(("[ ok ]  " if ok else "[FAIL]  ") + msg)
    if not ok:
        FAILS += 1


def read_fmap(path):
    with open(path) as f:
        lines = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    assert lines[0].startswith("IBSIMU-CYCL-FIELDMAP")
    nx, ny, nz = (int(v) for v in lines[1].split())
    x0, dx, y0, dy, z0, dz = (float(v) for v in lines[2].split())
    data = np.array([[float(v) for v in ln.split()] for ln in lines[3:]])
    return (nx, ny, nz), (x0, dx, y0, dy, z0, dz), data.reshape(nx, ny, nz, 3)


def main():
    R, L, E0 = 0.30, 0.60, 1.0e6
    nx = ny = 121
    nz = 3

    with tempfile.TemporaryDirectory() as d:
        pb.main(["--R", str(R), "--L", str(L), "--E0", str(E0),
                 "--nx", str(nx), "--ny", str(ny), "--nz", str(nz), "--outdir", d])
        (nxr, nyr, nzr), (x0, dx, y0, dy, z0, dz), E = read_fmap(
            os.path.join(d, "pillbox_E.txt"))
        _, _, B = read_fmap(os.path.join(d, "pillbox_B.txt"))

    check((nxr, nyr, nzr) == (nx, ny, nz), "E map header dims")
    check(abs(dx - 2*R/(nx-1)) < 1e-15, "E map x step consistent with [-R,R]")

    # --- 频率 ---
    f0 = pb.tm010_frequency(R)
    f_exp = pb.C_LIGHT*pb.X01/(2*np.pi*R)
    check(abs(f0 - f_exp) < 1e-9, "frequency formula f = c*x01/(2*pi*R)")
    print(f"       f_TM010 = {f0/1e6:.6f} MHz for R = {R} m")
    check(abs(pb.j0(np.array([pb.X01]))[0]) < 2e-8, "J0(x01) == 0 (Bessel impl + first zero)")
    check(abs(pb.j1(np.array([0.0]))[0]) < 1e-12, "J1(0) == 0")

    # --- 场结构 ---
    # 中心节点（x=y=0 处，nx/ny 为奇数 -> 必有该节点）
    ic = nx//2
    Ez_center = E[ic, ic, 0, 2]
    check(abs(Ez_center - E0)/E0 < 1e-6, f"Ez(r=0) = E0 ({Ez_center:.3f} V/m)")

    # 半径方向: r > 0.99R 处 Ez 应远小于 E0
    x = x0 + dx*np.arange(nx)
    y = y0 + dy*np.arange(ny)
    X, Y = np.meshgrid(x, y, indexing="ij")
    r = np.hypot(X, Y)
    Ez = E[:, :, 0, 2]
    edge = r > 0.99*R
    check(np.abs(Ez[edge]).max()/E0 < 0.05,
          f"Ez near wall is small (max {np.abs(Ez[edge]).max()/E0:.3e} * E0)")

    # 沿 theta=0 (y=0, x>0) 校验 Maxwell: |dEz/dx| = omega*|Btheta|
    # 避开近轴(1/x^2 使差分误差放大)与近壁区域
    jy = ny//2                  # y = 0
    ipos = slice(nx//2 + 2, nx - 2)
    xs = x[ipos]
    ez = E[ipos, jy, 0, 2]
    by = B[ipos, jy, 0, 1]     # theta=0 时 phi_hat=+y，故 Btheta = By
    dEzdx = np.gradient(ez, xs)
    omega = 2*np.pi*f0
    m = (xs > 0.05*R) & (xs < 0.90*R)
    lhs = np.abs(dEzdx[m])
    rhs = omega*np.abs(by[m])
    rel = np.abs(lhs - rhs).max()/rhs.max()
    print(f"       Maxwell |dEz/dx| vs omega*|Bth|: rel err {rel:.2e}")
    check(rel < 0.02, "Maxwell: |curl E| == omega*|B| (amplitude, 2%)")

    # --- 电场图能被 C++ 端解析的格式检查 ---
    check(abs(B[ic, jy, 0, 0]) < 1e-12, "Bx = 0 on the theta=0 axis")

    print(f"\n{'FAILED' if FAILS else 'ALL TESTS PASSED'} ({FAILS} failure{'s' if FAILS != 1 else ''})")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
