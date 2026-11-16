"""Kitaev-model Hamiltonians on arbitrary graphs, in Majorana form.

    H = sum_<i,j> [ -t (a_i^dag a_j + h.c.) + Delta (a_i a_j + h.c.) ]
        - mu sum_j (n_j - 1/2)
      = (i/4) sum_{a,b} A_ab gamma_a gamma_b ,      A real antisymmetric.

With c_j = gamma_{2j}, d_j = gamma_{2j+1}:
    bond (i -> j):  (i/2) [ (Delta - t) c_i d_j + (Delta + t) d_i c_j ]
    site j       :  -(mu/2) i c_j d_j

The same matrix A feeds (1) the qubit Hamiltonian of any encoding and
(2) the exact free-fermion (Bogoliubov-de Gennes) solution, so the two
calculations share no code beyond A itself.
"""
from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np
import scipy.sparse.linalg as spla
from qiskit.quantum_info import SparsePauliOp

from .encodings import Encoding

Edge = Tuple[int, int]


# ------------------------------------------------------------------ geometries
def chain_edges(n: int) -> List[Edge]:
    return [(j, j + 1) for j in range(n - 1)]


def ladder_edges(length: int, ordering: str = "rung") -> Tuple[int, List[Edge]]:
    """Two-leg ladder with `length` rungs.

    ordering='rung': site (leg, x) -> 2x + leg   (zig-zag, short strings on rungs)
    ordering='leg' : site (leg, x) -> leg*L + x  (leg by leg, long strings on rungs)
    """
    idx = (lambda leg, x: 2 * x + leg) if ordering == "rung" else (lambda leg, x: leg * length + x)
    edges = []
    for x in range(length):
        edges.append((idx(0, x), idx(1, x)))
        if x + 1 < length:
            edges.append((idx(0, x), idx(0, x + 1)))
            edges.append((idx(1, x), idx(1, x + 1)))
    return 2 * length, edges


def tjunction_edges(arm: int, ordering: str = "bar_first") -> Tuple[int, List[Edge]]:
    """T-junction: a horizontal bar of 2*arm+1 sites and a vertical leg of
    `arm` sites hanging off the centre.  This is the smallest geometry on
    which Majoranas can be braided (Alicea et al. 2011).

    ordering='bar_first': number the bar left to right, then the leg.
    ordering='leg_inline': left arm, leg (far end first), centre, right arm --
        keeps both left-centre and leg-centre bonds short.
    """
    nbar = 2 * arm + 1
    n = nbar + arm
    bar = list(range(nbar))
    leg = list(range(nbar, n))  # leg[0] touches the centre
    geo_edges = [(bar[k], bar[k + 1]) for k in range(nbar - 1)]
    geo_edges += [(bar[arm], leg[0])] + [(leg[k], leg[k + 1]) for k in range(arm - 1)]
    if ordering == "bar_first":
        perm = {s: s for s in range(n)}
    else:
        order = bar[:arm] + leg[::-1] + bar[arm:]
        perm = {s: k for k, s in enumerate(order)}
    return n, [(perm[i], perm[j]) for i, j in geo_edges]


# ------------------------------------------------------------- model matrices
def majorana_matrix(n: int, edges: Sequence[Edge], t: float = 1.0, delta: float = 1.0,
                    mu: float = 0.0) -> np.ndarray:
    a = np.zeros((2 * n, 2 * n))
    for i, j in edges:
        a[2 * i, 2 * j + 1] += delta - t
        a[2 * i + 1, 2 * j] += delta + t
    for j in range(n):
        a[2 * j, 2 * j + 1] += -mu
    return a - a.T


def hamiltonian(enc: Encoding, amat: np.ndarray) -> SparsePauliOp:
    """Qubit Hamiltonian  (i/2) sum_{a<b} A_ab gamma_a gamma_b  in encoding `enc`."""
    terms = []
    for a in range(2 * enc.n):
        for b in range(a + 1, 2 * enc.n):
            if abs(amat[a, b]) > 1e-14:
                terms.append(0.5j * amat[a, b] * (enc.gamma(a) @ enc.gamma(b)))
    h = SparsePauliOp.sum(terms).simplify(atol=1e-12)
    assert np.allclose(h.coeffs.imag, 0, atol=1e-10), "Hamiltonian must be Hermitian"
    return SparsePauliOp(h.paulis, h.coeffs.real)


# ------------------------------------------------- exact free-fermion solution
def bdg_levels(amat: np.ndarray) -> np.ndarray:
    """Single-particle energies eps_k >= 0;  H = sum_k eps_k (f_k^dag f_k - 1/2)."""
    ev = np.linalg.eigvalsh(1j * amat)
    return np.sort(ev[ev.size // 2:])


def bdg_ground_energy(amat: np.ndarray) -> float:
    return -0.5 * float(np.sum(bdg_levels(amat)))


def bdg_edge_correlator(amat: np.ndarray, a: int = 0, b: int = -1) -> float:
    """<i gamma_a gamma_b> in the free-fermion ground state.

    <gamma_a gamma_b> = 2 (P_+)_ab with P_+ the projector on the positive
    eigenspace of the Hermitian matrix iA.

    The lowest pair of levels (+-eps_0) is treated separately: deep in the
    topological phase eps_0 underflows to zero and a generic eigensolver
    returns an arbitrary basis of that 2D space.  The space is closed under
    complex conjugation, so it has a real orthonormal basis (r1, r2) -- the
    two Majorana zero modes -- and the fermion-parity eigenstates always
    occupy (r1 +- i r2)/sqrt(2), whatever eps_0 is.
    """
    m = 1j * amat
    w, v = np.linalg.eigh(m)
    order = np.argsort(np.abs(w))
    zero, rest = order[:2], order[2:]
    vp = v[:, rest][:, w[rest] > 0]
    ps = (v[:, zero] @ v[:, zero].conj().T).real
    _, r = np.linalg.eigh(ps)
    r1, r2 = r[:, -1], r[:, -2]
    f = (r1 + 1j * r2) / np.sqrt(2)
    if (f.conj() @ m @ f).real < 0:      # pick the occupied-energy-lowering sign
        f = f.conj()
    g = 2 * (vp @ vp.conj().T + np.outer(f, f.conj()))
    return float((1j * g[a, b]).real)


# ---------------------------------------------------- exact diagonalisation
def lowest_levels(h: SparsePauliOp, k: int = 4) -> np.ndarray:
    n = h.num_qubits
    if n <= 10:
        return np.linalg.eigvalsh(h.to_matrix())[:k]
    m = h.to_matrix(sparse=True)
    vals = spla.eigsh(m, k=k, which="SA", tol=1e-13, ncv=max(40, 4 * k), return_eigenvectors=False)
    return np.sort(vals)


def ground_states(h: SparsePauliOp, k: int = 2):
    n = h.num_qubits
    if n <= 10:
        w, v = np.linalg.eigh(h.to_matrix())
        return w[:k], v[:, :k]
    m = h.to_matrix(sparse=True)
    w, v = spla.eigsh(m, k=k, which="SA", tol=0, ncv=max(60, 6 * k))
    order = np.argsort(w)
    # Lanczos vectors of a numerically degenerate pair are accurate
    # eigenvectors but need not be orthogonal to each other: fix that.
    q, _ = np.linalg.qr(v[:, order])
    return w[order], q
