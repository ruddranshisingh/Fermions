#!/usr/bin/env python3
"""Run experiments E1-E4 and write results/*.json.

    python scripts/run_experiments.py            # everything (~10 min on a laptop)
    python scripts/run_experiments.py e2 e4      # a subset
    python scripts/run_experiments.py --quick    # smaller sizes, ~1 min
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from strings_attached import *  # noqa: E402
from strings_attached import cost, errors, noise  # noqa: E402
from strings_attached.encodings import ENCODINGS  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "..", "results")
QUICK = "--quick" in sys.argv
ENC = list(ENCODINGS)


def save(name, data):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name + ".json"), "w") as f:
        json.dump(data, f, indent=1)
    print(f"  -> results/{name}.json")


# =============================================================== E1 validation
def e1():
    t, delta = 1.0, 1.0
    n = 8 if QUICK else 12
    mus = [round(float(m), 4) for m in np.linspace(0.05, 4.0, 20 if QUICK else 40)]
    sweep = {"n": n, "t": t, "delta": delta, "mu": mus, "bdg": {"e0": [], "splitting": [], "gap": [], "edge": []},
             "ed": {e: {"e0": [], "splitting": [], "gap": [], "edge": []} for e in ENC}}
    encs = {e: get_encoding(e, n) for e in ENC}
    edge_ops = {e: (1j * (encs[e].gamma(0) @ encs[e].gamma(2 * n - 1))).simplify().to_matrix(sparse=True) for e in ENC}
    par_ops = {e: encs[e].parity_operator().to_matrix(sparse=True) for e in ENC}

    def edge_magnitude(e, w, v):
        """|<i g_0 g_{2N-1}>| in the parity-resolved ground state.  Deep in the
        topological phase the two lowest states are degenerate to machine
        precision, so the eigensolver returns an arbitrary mixture; we first
        rotate the pair into fermion-parity eigenstates (the correlator has
        the same magnitude and opposite sign in the two)."""
        sub = v[:, :2]
        _, u = np.linalg.eigh(sub.conj().T @ (par_ops[e] @ sub))
        states = sub @ u
        energies = [w[0] * abs(u[0, k]) ** 2 + w[1] * abs(u[1, k]) ** 2 for k in range(2)]
        psi = states[:, int(np.argmin(energies))]
        return float(abs((psi.conj() @ (edge_ops[e] @ psi)).real))
    for mu in mus:
        amat = majorana_matrix(n, chain_edges(n), t, delta, mu)
        lv = bdg_levels(amat)
        sweep["bdg"]["e0"].append(bdg_ground_energy(amat))
        sweep["bdg"]["splitting"].append(float(lv[0]))
        sweep["bdg"]["gap"].append(float(lv[1]))
        sweep["bdg"]["edge"].append(abs(bdg_edge_correlator(amat)))
        for e in ENC:
            w, v = ground_states(hamiltonian(encs[e], amat), 3)
            d = sweep["ed"][e]
            d["e0"].append(float(w[0]))
            d["splitting"].append(float(w[1] - w[0]))
            d["gap"].append(float(w[2] - w[0]))
            d["edge"].append(edge_magnitude(e, w, v))
    dev = {e: {k: float(np.max(np.abs(np.array(sweep["ed"][e][k]) - np.array(sweep["bdg"][k]))))
               for k in ("e0", "splitting", "gap", "edge")} for e in ENC}
    sweep["max_abs_deviation_from_bdg"] = dev

    # full-spectrum agreement between encodings, on three geometries, Delta != t
    full = []
    for gname, (m, edges) in {"chain(8)": (8, chain_edges(8)), "ladder(2x4)": ladder_edges(4),
                              "T-junction(7)": tjunction_edges(2)}.items():
        amat = majorana_matrix(m, edges, 1.0, 0.6, 0.9)
        spectra = {e: np.linalg.eigvalsh(hamiltonian(get_encoding(e, m), amat).to_matrix()) for e in ENC}
        ref = spectra["jordan_wigner"]
        # many-body spectrum rebuilt from the free-fermion levels
        eps = bdg_levels(amat)
        mb = np.sort([np.sum(eps * (np.array(list(np.binary_repr(k, m)), dtype=float) - 0.5)) for k in range(2 ** m)])
        full.append({"geometry": gname, "levels": int(2 ** m),
                     "max_dev_between_encodings": {e: float(np.max(np.abs(spectra[e] - ref))) for e in ENC},
                     "max_dev_from_free_fermion": float(np.max(np.abs(ref - mb)))})

    # system-size scaling of ED vs BdG at mu = 1
    scal = {"n": [], "bdg_splitting": [], "ed_splitting": [], "bdg_e0": [], "ed_e0": []}
    for m in range(3, (11 if QUICK else 15)):
        amat = majorana_matrix(m, chain_edges(m), 1.0, 1.0, 1.0)
        lv = lowest_levels(hamiltonian(jordan_wigner(m), amat), 2)
        scal["n"].append(m)
        scal["bdg_splitting"].append(float(bdg_levels(amat)[0]))
        scal["ed_splitting"].append(float(lv[1] - lv[0]))
        scal["bdg_e0"].append(bdg_ground_energy(amat))
        scal["ed_e0"].append(float(lv[0]))
    save("e1_validation", {"sweep": sweep, "full_spectrum": full, "scaling": scal})


# ========================================================= E2 error dictionary
def e2():
    t, mu, eps = 1.0, 1.0, 0.02
    nmax = 10 if QUICK else 14
    ns = list(range(4, nmax + 1))
    # (key, encoding, pauli, which qubits)
    classes = []
    for p in "XYZ":
        classes.append((f"JW {p} (all qubits)", "jordan_wigner", p, "all"))
    for p in "XYZ":
        classes.append((f"Parity {p} (qubits 0..N-2)", "parity", p, "bulk"))
        classes.append((f"Parity {p} (qubit N-1)", "parity", p, "last"))
    for p in "XYZ":
        classes.append((f"BK {p} (all qubits)", "bravyi_kitaev", p, "all"))

    def sites(which, n):
        return {"all": range(n), "bulk": range(n - 1), "last": [n - 1]}[which]

    out = {"t": t, "mu": mu, "eps": eps, "n": ns, "classes": [c[0] for c in classes], "runs": {}}
    for delta in (1.0, 0.6):
        res = {c[0]: [] for c in classes}
        res["unperturbed"] = []
        res["gap_unperturbed"] = []
        for n in ns:
            amat = majorana_matrix(n, chain_edges(n), t, delta, mu)
            res["unperturbed"].append(float(bdg_levels(amat)[0]))
            res["gap_unperturbed"].append(float(bdg_levels(amat)[1]))
            hs = {e: hamiltonian(get_encoding(e, n), amat) for e in ("jordan_wigner", "parity", "bravyi_kitaev")}
            for key, e, p, which in classes:
                v = errors.perturbation(n, p, sites(which, n))
                res[key].append(errors.splitting(hs[e], v, eps)["splitting"])
        out["runs"][f"delta={delta}"] = res

    # exact predictions at Delta = t (see docs/physics.md):
    #  JW Z_j  -> shifts mu by 2 eps                      (still free fermions)
    #  JW Y_j  -> rotates the transverse field: mu_eff/2 = sqrt((mu/2)^2 + eps^2)
    pred = {"JW Z (all qubits)": [], "JW Y (all qubits)": []}
    for n in ns:
        pred["JW Z (all qubits)"].append(float(bdg_levels(majorana_matrix(n, chain_edges(n), t, t, mu + 2 * eps))[0]))
        mu_eff = 2 * np.sqrt((mu / 2) ** 2 + eps ** 2)
        pred["JW Y (all qubits)"].append(float(bdg_levels(majorana_matrix(n, chain_edges(n), t, t, mu_eff))[0]))
    out["predictions"] = pred

    # strength scan at fixed N: first order (slope 1) vs higher order
    n = 8 if QUICK else 10
    amat = majorana_matrix(n, chain_edges(n), t, t, mu)
    eps_grid = [float(x) for x in np.logspace(-3, -1, 9)]
    scan = {"n": n, "eps": eps_grid, "unperturbed": float(bdg_levels(amat)[0]), "curves": {}}
    for key, e, p, which in classes:
        h0 = hamiltonian(get_encoding(e, n), amat)
        v = errors.perturbation(n, p, sites(which, n))
        scan["curves"][key] = [errors.splitting(h0, v, x)["splitting"] for x in eps_grid]
    out["eps_scan"] = scan

    # random-sign disorder: h_j uniform in [-1, 1]
    rng = np.random.default_rng(7)
    nd = list(range(4, (9 if QUICK else 13)))
    dis = {"n": nd, "seeds": 4 if QUICK else 12, "median": {c[0]: [] for c in classes}}
    for n in nd:
        amat = majorana_matrix(n, chain_edges(n), t, t, mu)
        hs = {e: hamiltonian(get_encoding(e, n), amat) for e in ("jordan_wigner", "parity", "bravyi_kitaev")}
        hj = rng.uniform(-1, 1, size=(dis["seeds"], n))
        for key, e, p, which in classes:
            s = list(sites(which, n))
            vals = [errors.splitting(hs[e], errors.perturbation(n, p, s, hj[k][s]), eps)["splitting"]
                    for k in range(dis["seeds"])]
            dis["median"][key].append(float(np.median(vals)))
    out["disorder"] = dis

    # qubit-resolved dictionary at N = 8
    n = 8
    amat = majorana_matrix(n, chain_edges(n), t, t, mu)
    table = {}
    for e in ENC:
        enc = get_encoding(e, n)
        fo = errors.first_order_matrix(enc, hamiltonian(enc, amat))
        rows = []
        for p in "XYZ":
            for q in range(n):
                d = errors.describe_pauli(enc, errors.single_qubit_label(n, q, p))
                d.update({"pauli": p, "qubit": q, "first_order": fo[p][q]})
                rows.append(d)
        table[e] = {"rows": rows, "majorana_strings": [enc.label(a) for a in range(2 * n)],
                    "parity_operator": enc.parity_operator().paulis[0].to_label()[::-1]}
    out["dictionary"] = {"n": n, "mu": mu, "encodings": table}
    save("e2_error_dictionary", out)


# ====================================================================== E3 cost
def e3():
    t, delta, mu = 1.0, 0.6, 0.8
    geoms = {
        "chain (12)": {"natural": (12, chain_edges(12))},
        "ladder (2x6)": {"rung order": ladder_edges(6, "rung"), "leg order": ladder_edges(6, "leg")},
        "T-junction (10)": {"bar first": tjunction_edges(3, "bar_first"), "leg inline": tjunction_edges(3, "leg_inline")},
    }
    if QUICK:
        geoms = {"chain (8)": {"natural": (8, chain_edges(8))},
                 "T-junction (7)": {"bar first": tjunction_edges(2, "bar_first"), "leg inline": tjunction_edges(2, "leg_inline")}}
    conns = ["line", "heavy_hex", "grid", "all_to_all"]
    seeds = list(range(2 if QUICK else 5))
    rows = []
    for gname, orderings in geoms.items():
        for oname, (n, edges) in orderings.items():
            amat = majorana_matrix(n, edges, t, delta, mu)
            for e in ENC:
                h = hamiltonian(get_encoding(e, n), amat)
                w = cost.pauli_weights(h)
                for c in conns:
                    runs = [cost.routed_cost(h, c, s) for s in seeds]
                    row = {"geometry": gname, "ordering": oname, "n": n, "encoding": e, "connectivity": c, **w}
                    for k in ("cx", "depth", "cx_depth"):
                        vals = [r[k] for r in runs]
                        row[k + "_mean"] = float(np.mean(vals))
                        row[k + "_min"] = int(min(vals))
                        row[k + "_max"] = int(max(vals))
                    rows.append(row)
                print(f"   {gname:16s} {oname:11s} {e:14s} done", flush=True)
    # where do the log-weight encodings start to pay off?  Pauli weight of a
    # nearest-neighbour bond and of the longest possible bond (mode 0 <-> n-1).
    sizes = [4, 8, 16, 32] if QUICK else [4, 6, 8, 12, 16, 24, 32, 48, 64]
    scaling = {"n": sizes, "nearest_mean": {e: [] for e in ENC}, "longest": {e: [] for e in ENC}}
    for n in sizes:
        for e in ENC:
            enc = get_encoding(e, n)
            wt = lambda a, b: sum(c != "I" for c in (enc.gamma(a) @ enc.gamma(b)).paulis[0].to_label())
            scaling["nearest_mean"][e].append(float(np.mean([wt(2 * j + 1, 2 * j + 2) for j in range(n - 1)])))
            scaling["longest"][e].append(int(wt(1, 2 * n - 2)))
    save("e3_cost", {"t": t, "delta": delta, "mu": mu, "seeds": len(seeds), "rows": rows, "weight_scaling": scaling})


# ===================================================================== E4 noise
def e4():
    n, t, mu, dt, steps = 6, 1.0, 1.0, 0.2, 4
    amat = majorana_matrix(n, chain_edges(n), t, t, mu)
    p2s = [float(x) for x in np.logspace(-4, -1.5, 6 if QUICK else 11)]
    out = {"n": n, "mu": mu, "dt": dt, "steps": steps, "p2": p2s, "curves": {}, "circuits": {}, "steps_scan": {}}
    for conn, line in (("line", True), ("all_to_all", False)):
        out["curves"][conn] = {}
        out["circuits"][conn] = {}
        for e in ENC:
            enc = get_encoding(e, n)
            h = hamiltonian(enc, amat)
            psi0 = noise.even_parity_ground_state(enc, h)
            tqc = noise.build_circuit(h, steps, dt, line=line)
            ops = tqc.count_ops()
            out["circuits"][conn][e] = {"cx": int(ops.get("cx", 0)), "sx_x": int(ops.get("sx", 0) + ops.get("x", 0)),
                                        "depth": int(tqc.depth()),
                                        "parity_check_weight": sum(c != "I" for c in enc.parity_operator().paulis[0].to_label())}
            runs = [noise.run(enc, h, tqc, psi0, p) for p in p2s]
            out["curves"][conn][e] = {k: [r[k] for r in runs] for k in runs[0]}
    # depth scan at fixed noise on the line
    p_fixed = 3e-3
    out["steps_scan"] = {"p2": p_fixed, "steps": list(range(1, 5 if QUICK else 9)), "curves": {}}
    for e in ENC:
        enc = get_encoding(e, n)
        h = hamiltonian(enc, amat)
        psi0 = noise.even_parity_ground_state(enc, h)
        runs = [noise.run(enc, h, noise.build_circuit(h, s, dt), psi0, p_fixed) for s in out["steps_scan"]["steps"]]
        out["steps_scan"]["curves"][e] = {k: [r[k] for r in runs] for k in runs[0]}
    save("e4_noise", out)


if __name__ == "__main__":
    todo = [a for a in sys.argv[1:] if not a.startswith("--")] or ["e1", "e2", "e3", "e4"]
    for name in todo:
        t0 = time.time()
        print(f"[{name.upper()}] running ...", flush=True)
        globals()[name]()
        print(f"[{name.upper()}] done in {time.time() - t0:.1f}s", flush=True)
