#!/usr/bin/env python3
"""Bunched-beam space-charge study summary figure (PSI Ring injection, RF off).

Data sources (produced by ``./tests/cycl_sc_beam --full``)::

    cycl_sc_beam_sample_init_t0.csv     launch samples (xi,xi',eta,eta',zeta)
    cycl_sc_beam_sample_nosc_t<K>.csv   reference run, per-turn samples
    cycl_sc_beam_sample_sc_t<K>.csv     space-charge run, per-turn samples
    cycl_sc_beam_envelope_{nosc,sc}.csv intra-turn envelope (200 rows/turn)
    cycl_sc_beam_convergence.csv        |dE|peak/|E|peak per Vlasov iteration
    cycl_sc_beam_scprofile.csv          self-field profiles E_xi(xi), E_zeta(zeta)
    cycl_sc_beam_map_probe.csv          tunes (sector-section differential map)
    cycl_sc_beam_summary.csv            run parameters / headline numbers

Usage::

    python3 examples/cyclotron/plot_sc_beam.py -o docs/img/cyclotron_sc_beam.png

The figure covers a **single 1 mA bunch** (Qb = I/f_RF) in the PSI Ring field
map with the RF switched off, i.e. an injection-region space-charge study.
RF acceleration (190 turns) is out of scope here; see docs/WORK_LOG.md.

IBSimu-Cycl is a derivative work of IBSimu, licensed under GPL-3.0-or-later.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]


def find(name_candidates, dirs):
    """Locate the first existing file over candidate names / directories."""
    for d in dirs:
        for n in name_candidates:
            p = Path(d) / n
            if p.exists():
                return p
    return None


def load_csv(path):
    return np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))


def load_summary(path):
    kv = {}
    if path is None:
        return kv
    with open(path, newline="") as f:
        for row in csv.reader(f):
            if len(row) == 2 and row[0] != "key":
                try:
                    kv[row[0]] = float(row[1])
                except ValueError:
                    kv[row[0]] = row[1]
    return kv


def rollmean(x, w):
    if len(x) < w:
        return x
    k = np.ones(w) / w
    return np.convolve(x, k, mode="same")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default=None,
                    help="directory with cycl_sc_beam_* outputs "
                         "(default: CWD, then ./tests, then ../tests)")
    ap.add_argument("-o", "--out", default="docs/img/cyclotron_sc_beam.png")
    ap.add_argument("--dpi", type=int, default=150)
    args = ap.parse_args(argv)

    dirs = []
    if args.data_dir:
        dirs.append(args.data_dir)
    dirs += [os.getcwd(), os.path.join(os.getcwd(), "tests"),
             os.path.join(os.getcwd(), "..", "tests")]

    def f(names):
        return find(names, dirs)

    p_init = f(["cycl_sc_beam_sample_init_t0.csv"])
    p_sum = f(["cycl_sc_beam_summary.csv"])
    p_env_n = f(["cycl_sc_beam_envelope_nosc.csv"])
    p_env_s = f(["cycl_sc_beam_envelope_sc.csv"])
    p_conv = f(["cycl_sc_beam_convergence.csv"])
    p_prof = f(["cycl_sc_beam_scprofile.csv"])
    p_probe = f(["cycl_sc_beam_map_probe.csv"])

    def last_turn_csv(tag):
        best = None
        best_turn = -1
        for d in dirs:
            for cand in glob.glob(os.path.join(d, f"cycl_sc_beam_sample_{tag}_t*.csv")):
                m = re.search(r"_t(\d+)\.csv$", cand)
                if m and int(m.group(1)) > best_turn:
                    best_turn = int(m.group(1))
                    best = Path(cand)
        return best

    p_sc = last_turn_csv("sc")
    p_nosc = last_turn_csv("nosc")

    if p_init is None or p_sc is None:
        raise SystemExit("cycl_sc_beam sample CSVs not found — run "
                         "./tests/cycl_sc_beam --full first")

    kv = load_summary(p_sum)
    d0 = load_csv(p_init)
    dsc = load_csv(p_sc)
    dnosc = load_csv(p_nosc) if p_nosc else None

    fig, ax = plt.subplots(2, 4, figsize=(17.5, 8.2))

    # ---------------------------------------------------------------- (1)
    a = ax[0, 0]
    a.plot(d0[:, 0]*1e3, d0[:, 2]*1e3, ".", ms=1, color="0.75", label="launch")
    a.plot(dsc[:, 0]*1e3, dsc[:, 2]*1e3, ".", ms=1, color="#1f77b4",
           label="final (with SC)")
    th = np.linspace(0, 2*np.pi, 200)
    r0 = 1e3*kv.get("beam_r_m", 1e-3)
    a.plot(r0*np.cos(th), r0*np.sin(th), "--", color="k", lw=0.8,
           label=f"emission disk r={r0:.1f} mm")
    a.set_xlabel(r"$\xi$ (mm)"); a.set_ylabel(r"$\eta$ (mm)")
    a.set_title("(a) transverse distribution")
    a.set_aspect("equal"); a.legend(fontsize=7, loc="upper right")

    # ---------------------------------------------------------------- (2)
    a = ax[0, 1]
    a.plot(d0[:, 4]*1e3, d0[:, 0]*1e3, ".", ms=1, color="0.75")
    a.plot(dsc[:, 4]*1e3, dsc[:, 0]*1e3, ".", ms=1, color="#1f77b4")
    a.axvline(0, color="k", lw=0.5)
    a.set_xlabel(r"$\zeta$ (mm)"); a.set_ylabel(r"$\xi$ (mm)")
    a.set_title(r"(b) $\zeta$–$\xi$ (frame-0 projection)")

    # ---------------------------------------------------------------- (3)
    a = ax[0, 2]
    for p, tag, c in ((p_env_n, "no SC", "0.55"), (p_env_s, "with SC", "#1f77b4")):
        if p is None:
            continue
        e = load_csv(p)
        a.plot(e[:, 1], e[:, 2], color=c, label=tag + r" $\sigma_\xi$")
        a.plot(e[:, 1], e[:, 3], color=c, ls="--", label=tag + r" $\sigma_\eta$")
    a.set_xlabel("t (ns)  [turn 1]"); a.set_ylabel(r"$\sigma$ (mm)")
    a.set_title(r"(c) transverse envelope")
    a.legend(fontsize=7)

    # ---------------------------------------------------------------- (4)
    a = ax[0, 3]
    for p, tag, c in ((p_env_n, "no SC", "0.55"), (p_env_s, "with SC", "#1f77b4")):
        if p is None:
            continue
        e = load_csv(p)
        a.plot(e[:, 1], e[:, 4], color=c, label=tag)
    a.set_xlabel("t (ns)  [turn 1]"); a.set_ylabel(r"$\sigma_\zeta$ (mm)")
    a.set_title(r"(d) longitudinal envelope")
    a.legend(fontsize=7)

    # ---------------------------------------------------------------- (5)
    a = ax[1, 0]
    for p, tag, c in ((p_env_n, "no SC", "0.55"), (p_env_s, "with SC", "#1f77b4")):
        if p is None:
            continue
        e = load_csv(p)
        a.semilogy(e[:, 1], e[:, 5], color=c, lw=0.5, alpha=0.4)
        a.semilogy(e[:, 1], rollmean(e[:, 5], 21), color=c, lw=1.6,
                   label=tag + r" $\varepsilon_r$")
        a.semilogy(e[:, 1], np.maximum(e[:, 6], 1e-12), color=c, lw=0.5,
                   ls=":", alpha=0.5)
        a.semilogy(e[:, 1], np.maximum(rollmean(e[:, 6], 21), 1e-12), color=c,
                   lw=1.4, ls="--", label=tag + r" $\varepsilon_z$")
    a.set_xlabel("t (ns)  [turn 1]"); a.set_ylabel(r"$\varepsilon$ (m$\cdot$rad)")
    a.set_title("(e) RMS emittance (raw + 21-sample average)")
    a.legend(fontsize=7)

    # ---------------------------------------------------------------- (6)
    a = ax[1, 1]
    if p_conv is not None:
        c = load_csv(p_conv)
        a.semilogy(c[:, 0], c[:, 1], "o-", color="#d62728",
                   label=r"$|\Delta E|_{peak}/|E|_{peak}$")
        a.semilogy(c[:, 0], c[:, 2], "s--", color="#2ca02c", label="RMS")
    a.set_xlabel("Vlasov iteration"); a.set_ylabel("SC field change")
    a.set_title("(f) self-consistency iteration")
    a.legend(fontsize=7)

    # ---------------------------------------------------------------- (7)
    a = ax[1, 2]
    if p_prof is not None:
        with open(p_prof) as fh:
            rd = list(csv.reader(fh))
        hdr = rd[0]
        xi = [(float(r[1]), float(r[2])) for r in rd[1:] if r[0] == "xi"]
        ze = [(float(r[1]), float(r[2])) for r in rd[1:] if r[0] == "zeta"]
        if xi:
            x, y = np.array(xi).T
            a.plot(x, y, color="#1f77b4", label=r"$E_\xi(\xi)$")
        if ze:
            x, y = np.array(ze).T
            a.plot(x, y, color="#d62728", label=r"$E_\zeta(\zeta)$")
        a.legend(fontsize=7)
    a.set_xlabel("coordinate (mm)"); a.set_ylabel("E (V/m)")
    a.set_title("(g) self-field profiles (segment 0)")

    # ---------------------------------------------------------------- (8)
    a = ax[1, 3]
    labels, vals, colors = [], [], []
    if p_probe is not None:
        with open(p_probe) as fh:
            rd = list(csv.reader(fh))
        for r in rd[1:]:
            if len(r) < 3:
                continue
            run, nr, nz = r[0], float(r[1]), float(r[2])
            labels += [f"{run}\n" + r"$\nu_r$", f"{run}\n" + r"$\nu_z$"]
            vals += [nr, nz]
            colors += ["#1f77b4", "#2ca02c"]
    else:
        vals = [kv.get("nu_r_nosc", 0), kv.get("nu_z_nosc", 0),
                kv.get("nu_r_sc", 0), kv.get("nu_z_sc", 0)]
        labels = ["no SC\n" + r"$\nu_r$", "no SC\n" + r"$\nu_z$",
                  "with SC\n" + r"$\nu_r$", "with SC\n" + r"$\nu_z$"]
        colors = ["#1f77b4", "#2ca02c", "#1f77b4", "#2ca02c"]
    x = np.arange(len(vals))
    a.bar(x, vals, color=colors, alpha=0.8)
    a.set_xticks(x); a.set_xticklabels(labels, fontsize=7)
    a.set_ylim(0, max(1.25, 1.08*max(vals)))
    dn_r = kv.get("dnu_r", float("nan"))
    dn_z = kv.get("dnu_z", float("nan"))
    a.set_title(f"(h) tunes (mod 1);  "
                fr"$\Delta\nu_r$={dn_r:+.4f}, $\Delta\nu_z$={dn_z:+.4f}")

    # ------------------------------------------------------------- header
    n = kv.get("N", float("nan"))
    qb = kv.get("Qb_C", float("nan"))
    L = 1e3*kv.get("bunch_len_m", float("nan"))
    rbeam = 1e3*kv.get("beam_r_m", float("nan"))
    rco = kv.get("r_co_m", float("nan"))
    fig.suptitle(
        "PSI Ring 72 MeV bunched-beam space-charge study "
        f"(RF off) —  N={n:.3g} macro-particles, Qb={qb*1e12:.1f} pC "
        f"(1 mA / 50.65 MHz), beam r={rbeam:.1f} mm, L={L:.0f} mm, "
        f"closed orbit r={rco:.3f} m",
        fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
