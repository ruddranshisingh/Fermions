"""Fermion-to-qubit encodings, written as lists of Majorana operators.

An encoding of n fermionic modes is fully specified by 2n Pauli strings
gamma_0 ... gamma_{2n-1} that are Hermitian, square to the identity and
mutually anticommute.  Mode j is  a_j = (gamma_{2j} + i gamma_{2j+1}) / 2.

Linear encodings (Jordan-Wigner, parity, Bravyi-Kitaev) store the qubit
bit-string  b = A n (mod 2)  for an invertible binary matrix A, and the
Majoranas follow from three index sets (Seeley, Richard & Love 2012):

    U(j) : qubits that flip when n_j flips       -> column j of A
    P(j) : qubits whose parity gives sum_{k<j} n_k -> row of (prefix . A^-1)
    F(j) : qubits whose parity gives n_j           -> row j of A^-1

    gamma_{2j}   =   X_U Z_P
    gamma_{2j+1} = i X_U Z_P Z_F

The ternary-tree encoding (Jiang et al. 2020) is not linear; its Majoranas
are the root-to-leaf paths of a ternary tree.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List

import numpy as np
from qiskit.quantum_info import Pauli, SparsePauliOp


# ----------------------------------------------------------------- GF(2) utils
def gf2_inv(a: np.ndarray) -> np.ndarray:
    """Inverse of a square binary matrix over GF(2) (Gauss-Jordan)."""
    n = a.shape[0]
    m = np.concatenate([a.astype(np.uint8) % 2, np.eye(n, dtype=np.uint8)], axis=1)
    for col in range(n):
        piv = next((r for r in range(col, n) if m[r, col]), None)
        if piv is None:
            raise ValueError("matrix is singular over GF(2)")
        m[[col, piv]] = m[[piv, col]]
        for r in range(n):
            if r != col and m[r, col]:
                m[r] ^= m[col]
    return m[:, n:]


def gf2_solve(m: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Solve m c = v over GF(2) for square invertible m."""
    return (gf2_inv(m).astype(int) @ v.astype(int)) % 2


# ------------------------------------------------------------------- Encoding
@dataclass
class Encoding:
    name: str
    n: int
    majoranas: List[SparsePauliOp]  # 2n single-term operators (phase in coeff)

    # -- basic objects ------------------------------------------------------
    def gamma(self, a: int) -> SparsePauliOp:
        return self.majoranas[a]

    def label(self, a: int) -> str:
        """Pauli string of gamma_a with qubit 0 on the LEFT (physics order)."""
        return self.majoranas[a].paulis[0].to_label()[::-1]

    def weight(self, a: int) -> int:
        return sum(ch != "I" for ch in self.label(a))

    def annihilation(self, j: int) -> SparsePauliOp:
        return (0.5 * (self.gamma(2 * j) + 1j * self.gamma(2 * j + 1))).simplify()

    def number(self, j: int) -> SparsePauliOp:
        # n_j = (1 + i gamma_{2j} gamma_{2j+1}) / 2
        eye = SparsePauliOp("I" * self.n)
        return (0.5 * (eye + 1j * (self.gamma(2 * j) @ self.gamma(2 * j + 1)))).simplify()

    def parity_operator(self) -> SparsePauliOp:
        """Total fermion parity (-1)^N_f = prod_j (-i gamma_{2j} gamma_{2j+1})."""
        op = SparsePauliOp("I" * self.n)
        for j in range(self.n):
            op = (op @ (-1j * (self.gamma(2 * j) @ self.gamma(2 * j + 1)))).simplify()
        return op

    # -- consistency --------------------------------------------------------
    def check(self) -> None:
        """Raise if the Majoranas do not satisfy {g_a, g_b} = 2 delta_ab."""
        eye = SparsePauliOp("I" * self.n)
        for a in range(2 * self.n):
            for b in range(a, 2 * self.n):
                anti = (self.gamma(a) @ self.gamma(b) + self.gamma(b) @ self.gamma(a)).simplify(atol=1e-12)
                target = 2 * eye if a == b else 0 * eye
                if not np.allclose((anti - target).simplify(atol=1e-12).coeffs, 0):
                    raise AssertionError(f"{self.name}: bad anticommutator ({a},{b})")

    # -- decoding a qubit Pauli back into fermion language --------------------
    def symplectic(self) -> np.ndarray:
        """(2n x 2n) binary matrix; column a is the [x|z] vector of gamma_a."""
        cols = []
        for g in self.majoranas:
            p = g.paulis[0]
            cols.append(np.concatenate([p.x, p.z]).astype(np.uint8))
        return np.array(cols).T

    def decode(self, pauli_label: str) -> List[int]:
        """Majorana indices whose product equals the Pauli string (up to phase).

        `pauli_label` has qubit 0 on the left, e.g. 'IXII' is X on qubit 1.
        Every Pauli string is +-(i^k) times a unique product of Majoranas,
        because the 2n Majoranas generate the whole Pauli group.
        """
        p = Pauli(pauli_label[::-1])
        v = np.concatenate([p.x, p.z]).astype(np.uint8)
        c = gf2_solve(self.symplectic(), v)
        return [int(a) for a in np.nonzero(c)[0]]


def _from_linear(name: str, amat: np.ndarray) -> Encoding:
    n = amat.shape[0]
    ainv = gf2_inv(amat)
    majs: List[SparsePauliOp] = []
    for j in range(n):
        u = amat[:, j].astype(bool)
        prefix = np.zeros(n, dtype=int)
        prefix[:j] = 1
        p = ((prefix @ ainv.astype(int)) % 2).astype(bool)
        f = ainv[j].astype(bool)
        zero = np.zeros(n, dtype=bool)
        xu = SparsePauliOp(Pauli((zero, u)))
        zp = SparsePauliOp(Pauli((p, zero)))
        zf = SparsePauliOp(Pauli((f, zero)))
        majs.append((xu @ zp).simplify())
        majs.append((1j * (xu @ zp @ zf)).simplify())
    return Encoding(name, n, majs)


def jordan_wigner(n: int) -> Encoding:
    """b_j = n_j.  Local occupation, non-local parity (Z strings)."""
    return _from_linear("jordan_wigner", np.eye(n, dtype=np.uint8))


def parity(n: int) -> Encoding:
    """b_j = n_0 + ... + n_j.  Local parity, non-local occupation updates."""
    return _from_linear("parity", np.tril(np.ones((n, n), dtype=np.uint8)))


def bravyi_kitaev_matrix(n: int) -> np.ndarray:
    size = 1
    a = np.ones((1, 1), dtype=np.uint8)
    while size < n:
        new = np.zeros((2 * size, 2 * size), dtype=np.uint8)
        new[:size, :size] = a
        new[size:, size:] = a
        new[-1, :size] = 1  # last qubit of the block stores the block parity
        a, size = new, 2 * size
    return a[:n, :n]


def bravyi_kitaev(n: int) -> Encoding:
    """Fenwick-tree compromise: both occupation and parity are O(log n)."""
    return _from_linear("bravyi_kitaev", bravyi_kitaev_matrix(n))


def ternary_tree(n: int) -> Encoding:
    """Majoranas = root-to-leaf paths of a complete ternary tree on n qubits.

    Qubit k has child slots 3k+1 (X branch), 3k+2 (Y branch), 3k+3 (Z branch).
    There are 2n+1 paths; any two differ first at one shared qubit where they
    carry different Paulis, so all of them anticommute.  We drop the last path
    and pair the remaining 2n in order.
    """
    paths: List[str] = []

    def walk(node: int, prefix: Dict[int, str]) -> None:
        for b, ch in enumerate("XYZ"):
            here = dict(prefix)
            here[node] = ch
            child = 3 * node + 1 + b
            if child < n:
                walk(child, here)
            else:
                paths.append("".join(here.get(q, "I") for q in range(n)))

    walk(0, {})
    assert len(paths) == 2 * n + 1
    majs = [SparsePauliOp(lbl[::-1]) for lbl in paths[: 2 * n]]
    return Encoding("ternary_tree", n, majs)


ENCODINGS: Dict[str, Callable[[int], Encoding]] = {
    "jordan_wigner": jordan_wigner,
    "parity": parity,
    "bravyi_kitaev": bravyi_kitaev,
    "ternary_tree": ternary_tree,
}

SHORT = {"jordan_wigner": "JW", "parity": "Parity", "bravyi_kitaev": "BK", "ternary_tree": "Ternary tree"}


def get_encoding(name: str, n: int) -> Encoding:
    return ENCODINGS[name](n)
