#!/usr/bin/env python3
"""field_to_fieldmap3d.py 的自检（无需 Elmer）。

用法: python3 tests/test_field_to_fieldmap3d.py

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # adapters/magnet3d

import field_to_fieldmap3d as conv                 # noqa: E402

FAILS = 0


def check(ok, msg):
    global FAILS
    print(("[ ok ]  " if ok else "[FAIL]  ") + msg)
    if not ok:
        FAILS += 1


def read_fieldmap_header(path):
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f]
    assert lines[0].startswith("IBSIMU-CYCL-FIELDMAP")
    nx, ny, nz = (int(v) for v in lines[2].split())
    x0, dx, y0, dy, z0, dz = (float(v) for v in lines[4].split())
    return (nx, ny, nz), (x0, dx, y0, dy, z0, dz)


def hdr_ok(ok, msg):   # small helper to keep messages uniform
    check(ok, msg)


def test_uniform():
    rng = np.random.default_rng(0)
    pts = rng.uniform(-1, 1, size=(2000, 3))
    B = np.tile(np.array([0.1, -0.2, 1.5]), (pts.shape[0], 1))
    x = np.linspace(-1, 1, 9)
    y = np.linspace(-1, 1, 7)
    z = np.linspace(-1, 1, 5)
    g = conv.resample(pts, B, x, y, z, k=8, power=2)
    err = np.abs(g - np.array([0.1, -0.2, 1.5])).max()
    check(err < 1e-9, f"uniform field resampled exactly (max err {err:.2e})")

    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, "f.txt")
        conv.write_fieldmap(out, x, y, z, g)
        (nx, ny, nz), (x0, dx, y0, dy, z0, dz) = read_fieldmap_header(out)
        check((nx, ny, nz) == (9, 7, 5), "header dims match")
        check(abs(dx - (x[1] - x[0])) < 1e-15 and abs(dz - (z[1] - z[0])) < 1e-15,
              "header steps match")
        # 读取若干行做数值校验
        vals = np.loadtxt(out, comments="#", skiprows=5)
        check(vals.shape == (9 * 7 * 5, 3), "data row count == nx*ny*nz")
        check(np.abs(vals - np.array([0.1, -0.2, 1.5])).max() < 1e-8,
              "written values are correct")


def test_linear():
    rng = np.random.default_rng(1)
    pts = rng.uniform(-1, 1, size=(8000, 3))

    def field(p):
        return np.stack([0.1 * p[:, 0], 0.2 * p[:, 1], 1.0 - 0.05 * p[:, 2]], axis=1)

    B = field(pts)
    x = np.linspace(-0.8, 0.8, 11)
    y = np.linspace(-0.8, 0.8, 9)
    z = np.linspace(-0.8, 0.8, 7)
    g = conv.resample(pts, B, x, y, z, k=12, power=2)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    ref = np.stack([0.1 * X, 0.2 * Y, 1.0 - 0.05 * Z], axis=-1)
    rel = np.abs(g - ref).max() / max(1.0, np.abs(ref).max())
    check(rel < 1e-2, f"linear field IDW error small (rel {rel:.2e})")


def test_cli():
    with tempfile.TemporaryDirectory() as d:
        inp = os.path.join(d, "pts.txt")
        out = os.path.join(d, "fm.txt")
        rng = np.random.default_rng(2)
        pts = rng.uniform(-1, 1, size=(500, 3))
        with open(inp, "w") as f:
            f.write("# x y z Bx By Bz\n")
            for p in pts:
                f.write(f"{p[0]} {p[1]} {p[2]} 0 0 1.0\n")
        rc = conv.main([inp, out, "--nx", "5", "--ny", "5", "--nz", "5", "--k", "4"])
        check(rc == 0 and os.path.exists(out), "CLI runs and writes output")
        (nx, ny, nz), _ = read_fieldmap_header(out)
        check((nx, ny, nz) == (5, 5, 5), "CLI output header dims ok")


if __name__ == "__main__":
    test_uniform()
    test_linear()
    test_cli()
    print(f"\n{'FAILED' if FAILS else 'ALL TESTS PASSED'} ({FAILS} failure{'s' if FAILS != 1 else ''})")
    sys.exit(1 if FAILS else 0)
