#!/usr/bin/env python3
"""把散点/节点磁场数据重采样为规则笛卡尔网格，并写出 IBSimu-Cycl 场图格式。

输入：文本文件，每行 ``x y z Bx By Bz``（空白或逗号分隔，# 为注释）。
      —— 例如从 Elmer VTU 经 ParaView "Save Data (ASCII)" 导出的节点数据。
输出：IBSimu-Cycl 原生 ASCII 场图（``IBSIMU-CYCL-FIELDMAP 1``），
      可由 C++ 的 ``ibsimu_cycl::CFieldMap3D::load()`` 直接读取。

重采样使用反距离加权 (IDW)：对每个网格点取最近的 k 个源点，权重 1/d^power。

用法示例::

    python3 field_to_fieldmap3d.py nodes_B.txt out.txt \\
        --nx 141 --ny 141 --nz 61 --k 8 --power 2

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import sys
import numpy as np

MAGIC = "IBSIMU-CYCL-FIELDMAP"
VERSION = 1


def read_points(path):
    """读取 ``x y z Bx By Bz`` 文本，返回 (pts (N,3), vals (N,3))。"""
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            line = line.replace(",", " ")
            parts = line.split()
            if len(parts) < 6:
                continue
            rows.append([float(p) for p in parts[:6]])
    if not rows:
        raise ValueError(f"no data rows found in {path}")
    a = np.asarray(rows, dtype=float)
    return a[:, :3], a[:, 3:6]


def make_grid_axes(n, lo, hi, margin=0.0):
    """生成 n 个节点，范围 [lo, hi] 收缩/扩张 margin 比例。"""
    span = hi - lo
    lo2 = lo + margin * span
    hi2 = hi - margin * span
    return np.linspace(lo2, hi2, n)


def resample(pts, vals, x, y, z, k=8, power=2.0):
    """IDW 重采样到规则网格，返回 (nx, ny, nz, 3)。"""
    pts = np.asarray(pts, dtype=float)
    vals = np.asarray(vals, dtype=float)
    n = pts.shape[0]
    kk = int(min(k, n))

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    m = grid.shape[0]
    out = np.empty((m, 3), dtype=float)

    # 动态分块，限制 (chunk x n) 中间矩阵的内存
    chunk = max(1, int(4_000_000 // max(n, 1)))
    for s in range(0, m, chunk):
        g = grid[s:s + chunk]                                   # (c,3)
        d2 = ((g[:, None, :] - pts[None, :, :]) ** 2).sum(-1)    # (c,n)
        if kk < n:
            idx = np.argpartition(d2, kk - 1, axis=1)[:, :kk]    # (c,kk)
            dd2 = np.take_along_axis(d2, idx, axis=1)
        else:
            idx = np.tile(np.arange(n), (g.shape[0], 1))
            dd2 = d2
        exact = dd2 < 1.0e-18
        w = 1.0 / np.where(dd2 > 0.0, dd2, 1.0) ** (power / 2.0)
        w = np.where(exact, 1.0e300, w)                         # 命中源点 -> 取精确值
        wsum = w.sum(axis=1, keepdims=True)
        out[s:s + chunk] = np.einsum("cn,cni->ci", w, vals[idx]) / wsum

    return out.reshape(len(x), len(y), len(z), 3)


def write_fieldmap(path, x, y, z, fgrid):
    """写出 IBSimu-Cycl 原生场图（与 C++ CFieldMap3D::save_ascii 一致）。"""
    nx, ny, nz = len(x), len(y), len(z)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dz = z[1] - z[0]
    with open(path, "w") as f:
        f.write(f"{MAGIC} {VERSION}\n")
        f.write("# nx ny nz\n")
        f.write(f"{nx} {ny} {nz}\n")
        f.write("# x0 dx y0 dy z0 dz\n")
        f.write("{:.17g} {:.17g} {:.17g} {:.17g} {:.17g} {:.17g}\n".format(
            x[0], dx, y[0], dy, z[0], dz))
        f.write("# data in order i(x), j(y), k(z): Fx Fy Fz\n")
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    v = fgrid[i, j, k]
                    f.write("{:.9g} {:.9g} {:.9g}\n".format(v[0], v[1], v[2]))


def main(argv=None):
    ap = argparse.ArgumentParser(description="散点磁场 -> IBSimu-Cycl 规则网格场图")
    ap.add_argument("input", help="输入文本: x y z Bx By Bz")
    ap.add_argument("output", help="输出场图 (IBSIMU-CYCL-FIELDMAP)")
    ap.add_argument("--nx", type=int, required=True)
    ap.add_argument("--ny", type=int, required=True)
    ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--bounds", type=float, nargs=6, default=None,
                    metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
                    help="网格范围；缺省用数据包围盒")
    ap.add_argument("--margin", type=float, default=0.0,
                    help="范围收缩比例 (0..0.2)，避免外插到边缘")
    ap.add_argument("--k", type=int, default=8, help="IDW 邻居数")
    ap.add_argument("--power", type=float, default=2.0, help="IDW 指数")
    args = ap.parse_args(argv)

    pts, vals = read_points(args.input)

    if args.bounds is not None:
        xmin, xmax, ymin, ymax, zmin, zmax = args.bounds
    else:
        xmin, xmax = pts[:, 0].min(), pts[:, 0].max()
        ymin, ymax = pts[:, 1].min(), pts[:, 1].max()
        zmin, zmax = pts[:, 2].min(), pts[:, 2].max()

    x = make_grid_axes(args.nx, xmin, xmax, args.margin)
    y = make_grid_axes(args.ny, ymin, ymax, args.margin)
    z = make_grid_axes(args.nz, zmin, zmax, args.margin)

    fgrid = resample(pts, vals, x, y, z, k=args.k, power=args.power)
    write_fieldmap(args.output, x, y, z, fgrid)

    print(f"points: {pts.shape[0]}  grid: {args.nx}x{args.ny}x{args.nz}"
          f"  -> {args.output}")
    bmag = np.sqrt((fgrid ** 2).sum(-1))
    print(f"|B| range: {bmag.min():.6g} .. {bmag.max():.6g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
