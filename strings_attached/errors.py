"""The error dictionary: what does a qubit error mean for the fermions?

Three tools:
  describe_pauli      -- algebra: decode a qubit Pauli into Majorana language
  first_order_matrix  -- perturbation theory: how strongly does a single-qubit
                         Pauli split the two ground states at first order
  splitting           -- exact diagonalisation: the full splitting under a
                         static perturbation  eps * sum_j h_j P_j
"""
from __future__ import annotations

from typing import Dict, Iterable, Optional, Sequence

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from .encodings import Encoding
from .kitaev import ground_states, lowest_levels


def single_qubit_label(n: int, q: int, p: str) -> str:
    return "".join(p if k == q else "I" for k in range(n))


def describe_pauli(enc: Encoding, label: str) -> Dict:
    """Fermionic meaning of a Pauli string (qubit 0 on the left)."""
    idx = enc.decode(label)
    n = enc.n
    even = len(idx) % 2 == 0
    span = (max(idx) - min(idx) + 1) if idx else 0
    return {
        "majoranas": idx,
        "majorana_weight": len(idx),
        "parity": "even" if even else "odd",
        "span": span,                      # extent along the 2n-site Majorana chain
        "local": span <= 4,                # fits inside two neighbouring sites
        "touches_left_edge": 0 in idx,
        "touches_right_edge": (2 * n - 1) in idx,
    }


def perturbation(n: int, pauli: str, sites: Iterable[int], weights: Optional[Sequence[float]] = None) -> SparsePauliOp:
    sites = list(sites)
    weights = [1.0] * len(sites) if weights is None else list(weights)
    return SparsePauliOp.from_sparse_list([(pauli, [q], w) for q, w in zip(sites, weights)], n)


def splitting(h0: SparsePauliOp, v: Optional[SparsePauliOp] = None, eps: float = 0.0, k: int = 4) -> Dict:
    """Lowest levels of h0 + eps*v.  'splitting' = E1 - E0, 'gap' = E2 - E0."""
    h = h0 if v is None else (h0 + eps * v).simplify()
    lv = lowest_levels(h, k=k)
    return {"splitting": float(lv[1] - lv[0]), "gap": float(lv[2] - lv[0]), "levels": [float(x) for x in lv]}


def first_order_matrix(enc: Encoding, h0: SparsePauliOp) -> Dict[str, list]:
    """First-order splitting per unit strength for every single-qubit Pauli.

    Degenerate perturbation theory: project P onto the two lowest states of
    h0 and take the spread of the 2x2 matrix's eigenvalues.  ~2 means "this
    error acts as a logical operator on the Majorana qubit"; ~0 means the
    ground space does not see it at first order.
    """
    _, vecs = ground_states(h0, 2)
    out = {}
    for p in "XYZ":
        row = []
        for q in range(enc.n):
            op = SparsePauliOp(single_qubit_label(enc.n, q, p)[::-1]).to_matrix(sparse=True)
            m = vecs.conj().T @ (op @ vecs)
            ev = np.linalg.eigvalsh(m)
            row.append(float(ev[1] - ev[0]))
        out[p] = row
    return out
