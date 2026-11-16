#!/usr/bin/env python3
"""Turn results/*.json into the report figures in figures/*.png."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
R = lambda name: json.load(open(os.path.join(ROOT, "results", name + ".json")))
ENC = ["jordan_wigner", "parity", "bravyi_kitaev", "ternary_tree"]
NAME = {"jordan_wigner": "Jordan-Wigner", "parity": "Parity", "bravyi_kitaev": "Bravyi-Kitaev", "ternary_tree": "Ternary tree"}
COL = {"jordan_wigner": "#2a78d6", "parity": "#eb6834", "bravyi_kitaev": "#1baf7a", "ternary_tree": "#eda100"}
MARK = {"jordan_wigner": "o", "parity": "s", "bravyi_kitaev": "^", "ternary_tree": "D"}
PCOL = {"X": "#e34948", "Y": "#6250d6", "Z": "#008300"}
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
                     "axes.titlesize": 9.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
                     "figure.dpi": 160, "savefig.bbox": "tight", "lines.linewidth": 1.8, "lines.markersize": 4.5})


def out(fig, name):
    os.makedirs(os.path.join(ROOT, "figures"), exist_ok=True)
    fig.savefig(os.path.join(ROOT, "figures", name + ".png"))
    plt.close(fig)
    print("  -> figures/" + name + ".png")


def fig1():
    d = R("e1_validation")
    s = d["sweep"]
    mu = np.array(s["mu"])
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    ax[0].plot(mu, s["bdg"]["splitting"], color=INK, lw=1.2, label="exact free-fermion (BdG)")
    ax[0].plot(mu, s["bdg"]["gap"], color=INK, lw=1.2, ls="--")
    for k, e in enumerate(ENC):
        sl = slice(k, None, 4)
        ax[0].plot(mu[sl], np.array(s["ed"][e]["splitting"])[sl], MARK[e], color=COL[e], label=NAME[e] + " (ED)")
        ax[0].plot(mu[sl], np.array(s["ed"][e]["gap"])[sl], MARK[e], color=COL[e], mfc="white")
    ax[0].set_yscale("log"); ax[0].set_ylim(1e-13, 10)
    ax[0].axvline(2, color=MUTED, lw=0.8, ls=":")
    ax[0].text(2.05, 2e-5, "phase boundary\n$\\mu = 2t$", color=MUTED, fontsize=7.5)
    ax[0].text(0.3, 2.6, "bulk gap (open markers)", color=MUTED, fontsize=7.5)
    ax[0].set_xlabel("chemical potential $\\mu / t$"); ax[0].set_ylabel("energy above ground state")
    ax[0].set_title(f"(a) Ground-state splitting, N = {s['n']}")
    ax[0].legend(fontsize=7, loc="lower right")

    ax[1].plot(mu, s["bdg"]["edge"], color=INK, lw=1.2, label="BdG")
    for k, e in enumerate(ENC):
        sl = slice(k, None, 4)
        ax[1].plot(mu[sl], np.array(s["ed"][e]["edge"])[sl], MARK[e], color=COL[e], label=NAME[e])
    ax[1].plot(mu[mu < 2], 1 - (mu[mu < 2] / 2) ** 2, color=MUTED, lw=0.8, ls=":", label="$1-(\\mu/2t)^2$")
    ax[1].axvline(2, color=MUTED, lw=0.8, ls=":")
    ax[1].set_xlabel("chemical potential $\\mu / t$"); ax[1].set_ylabel("$|\\langle i\\gamma_0\\gamma_{2N-1}\\rangle|$")
    ax[1].set_title("(b) Edge-to-edge Majorana correlator"); ax[1].legend(fontsize=7)

    sc = d["scaling"]
    ax[2].plot(sc["n"], sc["bdg_splitting"], color=INK, lw=1.2, label="BdG")
    ax[2].plot(sc["n"], sc["ed_splitting"], "o", color=COL["jordan_wigner"], label="ED (Jordan-Wigner)")
    n = np.array(sc["n"], float)
    ax[2].plot(n, 1.5 * 0.5 ** n, color=MUTED, lw=0.8, ls=":", label="$\\propto (\\mu/2t)^N$")
    ax[2].set_yscale("log"); ax[2].set_xlabel("chain length N"); ax[2].set_ylabel("splitting $E_1 - E_0$")
    ax[2].set_title("(c) Exponential protection, $\\mu = t$"); ax[2].legend(fontsize=7)
    fig.tight_layout(); out(fig, "fig1_validation")


def fig2():
    d = R("e2_error_dictionary")
    n = d["n"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    run = d["runs"]["delta=1.0"]
    panels = [("Jordan-Wigner", [("JW %s (all qubits)", "-", "o", "")]),
              ("Parity", [("Parity %s (qubits 0..N-2)", "-", "o", " bulk"), ("Parity %s (qubit N-1)", "--", "s", " on qubit N-1")]),
              ("Bravyi-Kitaev", [("BK %s (all qubits)", "-", "o", "")])]
    for ax, (title, specs) in zip(axes, panels):
        ax.plot(n, run["unperturbed"], color=INK, lw=1.0, ls=":", label="no error")
        for fmt, ls, mk, suffix in specs:
            for p in "XYZ":
                ax.plot(n, run[fmt % p], ls=ls, marker=mk, color=PCOL[p], mfc=PCOL[p] if ls == "-" else "white",
                        label=f"{p} errors{suffix}")
        ax.set_yscale("log"); ax.set_ylim(3e-7, 1.5); ax.set_xlabel("chain length N"); ax.set_title(title)
        ax.legend(fontsize=6.8, loc="lower left", ncol=2)
    axes[0].set_ylabel("ground-state splitting $E_1-E_0$")
    fig.suptitle(f"Splitting under a static error field $\\varepsilon\\sum_j P_j$,  $\\varepsilon={d['eps']}$, $\\mu=t=\\Delta$",
                 x=0.01, ha="left", fontsize=10, fontweight="bold")
    fig.tight_layout(); out(fig, "fig2_splitting")


def fig2b():
    d = R("e2_error_dictionary")
    sc = d["eps_scan"]
    fig, ax = plt.subplots(1, 2, figsize=(8.4, 3.3))
    eps = np.array(sc["eps"])
    sel = [("JW X (all qubits)", PCOL["X"], "-", "JW X"), ("JW Y (all qubits)", PCOL["Y"], "-", "JW Y"),
           ("JW Z (all qubits)", PCOL["Z"], "-", "JW Z"), ("Parity X (qubit N-1)", PCOL["X"], "--", "Parity X, qubit N-1"),
           ("Parity Z (qubits 0..N-2)", PCOL["Z"], "--", "Parity Z, bulk")]
    for key, c, ls, lab in sel:
        ax[0].plot(eps, sc["curves"][key], ls=ls, marker="o", color=c, mfc=c if ls == "-" else "white", label=lab)
    ax[0].axhline(sc["unperturbed"], color=INK, lw=1, ls=":", label="no error")
    ax[0].plot(eps, 2 * eps * sc["n"], color=MUTED, lw=0.8, ls=":")
    ax[0].text(eps[2], 2 * eps[2] * sc["n"] * 1.5, "$2\\varepsilon N$", color=MUTED, fontsize=8)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("error strength $\\varepsilon / t$")
    ax[0].set_ylabel("splitting"); ax[0].set_title(f"(a) First order or not?  N = {sc['n']}"); ax[0].legend(fontsize=7)
    n = d["n"]
    run = d["runs"]["delta=1.0"]
    ax[1].plot(n, run["JW Y (all qubits)"], "o", color=PCOL["Y"], label="ED: JW Y errors")
    ax[1].plot(n, d["predictions"]["JW Y (all qubits)"], color=PCOL["Y"], lw=1, label="$\\mu_{\\rm eff}=2\\sqrt{(\\mu/2)^2+\\varepsilon^2}$")
    ax[1].plot(n, run["JW Z (all qubits)"], "s", color=PCOL["Z"], label="ED: JW Z errors")
    ax[1].plot(n, d["predictions"]["JW Z (all qubits)"], color=PCOL["Z"], lw=1, label="$\\mu_{\\rm eff}=\\mu+2\\varepsilon$")
    ax[1].set_yscale("log"); ax[1].set_xlabel("chain length N"); ax[1].set_ylabel("splitting")
    ax[1].set_title("(b) Protected errors = renormalised $\\mu$"); ax[1].legend(fontsize=7)
    fig.tight_layout(); out(fig, "fig2b_first_order")


def fig3():
    d = R("e2_error_dictionary")["dictionary"]
    n = d["n"]
    fig, axes = plt.subplots(1, 4, figsize=(11, 2.5))
    for ax, e in zip(axes, ENC):
        rows = d["encodings"][e]["rows"]
        m = np.array([[r["first_order"] for r in rows if r["pauli"] == p] for p in "XYZ"])
        im = ax.imshow(np.clip(m, 0, 2), cmap="Blues", vmin=0, vmax=2, aspect="auto")
        for i, p in enumerate("XYZ"):
            for q in range(n):
                r = next(x for x in rows if x["pauli"] == p and x["qubit"] == q)
                ax.text(q, i, f"{m[i, q]:.2f}\n{'odd' if r['parity'] == 'odd' else 'even'}", ha="center", va="center",
                        fontsize=5.6, color="white" if m[i, q] > 1.1 else INK)
        ax.set_xticks(range(n)); ax.set_yticks(range(3)); ax.set_yticklabels(list("XYZ"))
        ax.set_xlabel("qubit"); ax.set_title(NAME[e]); ax.grid(False)
        ax.set_xticks(np.arange(-.5, n), minor=True); ax.set_yticks(np.arange(-.5, 3), minor=True)
        ax.grid(which="minor", color="white", lw=2); ax.tick_params(which="minor", length=0)
    fig.colorbar(im, ax=axes, shrink=0.8, pad=0.01, label="first-order splitting / $\\varepsilon$")
    out(fig, "fig3_dictionary")


def fig4():
    d = R("e3_cost")
    rows = d["rows"]
    geoms = []
    for r in rows:
        if r["geometry"] not in geoms:
            geoms.append(r["geometry"])
    conns = ["line", "heavy_hex", "grid", "all_to_all"]
    fig, axes = plt.subplots(1, len(geoms) + 1, figsize=(12.5, 3.3))
    for ax, g in zip(axes, geoms):
        for k, e in enumerate(ENC):
            best = []
            for c in conns:
                cand = [r for r in rows if r["geometry"] == g and r["encoding"] == e and r["connectivity"] == c]
                b = min(cand, key=lambda r: r["cx_mean"])
                best.append(b)
            x = np.arange(len(conns)) + (k - 1.5) * 0.2
            y = [b["cx_mean"] for b in best]
            err = [[b["cx_mean"] - b["cx_min"] for b in best], [b["cx_max"] - b["cx_mean"] for b in best]]
            ax.bar(x, y, 0.18, color=COL[e], label=NAME[e], yerr=err, error_kw={"lw": 0.8, "ecolor": INK})
        ax.set_xticks(range(len(conns))); ax.set_xticklabels(["line", "heavy-hex", "grid", "all-to-all"])
        ax.set_title(g); ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("CNOTs per Trotter step (best ordering)"); axes[0].legend(fontsize=7)
    ws = d["weight_scaling"]
    ax = axes[-1]
    for e in ENC:
        ax.plot(ws["n"], ws["longest"][e], marker=MARK[e], color=COL[e], label=NAME[e])
    ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
    ax.set_xlabel("number of modes n"); ax.set_ylabel("Pauli weight"); ax.set_title("Longest bond (mode 0 to n-1)")
    fig.tight_layout(); out(fig, "fig4_cost")


def fig5():
    d = R("e4_noise")
    p = np.array(d["p2"])
    fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.3))
    for e in ENC:
        c = d["curves"]["line"][e]
        ax[0].plot(p, 1 - np.array(c["fidelity_raw"]), ls="--", marker=MARK[e], mfc="white", color=COL[e])
        ax[0].plot(p, 1 - np.array(c["fidelity_postselected"]), marker=MARK[e], color=COL[e], label=NAME[e])
        ax[1].plot(p, 1 - np.array(c["kept_fraction"]), marker=MARK[e], color=COL[e], label=NAME[e])
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("two-qubit depolarising rate $p_2$")
    ax[0].set_ylabel("infidelity $1-F$"); ax[0].legend(fontsize=7)
    ax[0].set_title("(a) Dashed: raw.  Solid: parity post-selected")
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("two-qubit depolarising rate $p_2$")
    ax[1].set_ylabel("fraction of runs discarded"); ax[1].set_title("(b) What the parity check catches")
    sc = d["steps_scan"]
    for e in ENC:
        c = sc["curves"][e]
        ax[2].plot(sc["steps"], c["fidelity_raw"], ls="--", marker=MARK[e], mfc="white", color=COL[e])
        ax[2].plot(sc["steps"], c["fidelity_postselected"], marker=MARK[e], color=COL[e], label=NAME[e])
    ax[2].set_xlabel("Trotter steps"); ax[2].set_ylabel("fidelity"); ax[2].set_title(f"(c) Depth scan at $p_2$ = {sc['p2']}")
    fig.tight_layout(); out(fig, "fig5_noise")


if __name__ == "__main__":
    for f in (fig1, fig2, fig2b, fig3, fig4, fig5):
        f()
