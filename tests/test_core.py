import numpy as np
import pytest
from qiskit.quantum_info import SparsePauliOp

from strings_attached import *
from strings_attached.encodings import ENCODINGS


@pytest.mark.parametrize("name", list(ENCODINGS))
@pytest.mark.parametrize("n", [1, 2, 5, 8])
def test_majoranas_anticommute(name, n):
    get_encoding(name, n).check()


def test_jw_matches_textbook_strings():
    enc = jordan_wigner(4)
    assert enc.label(0) == "XIII"
    assert enc.label(7) == "ZZZY"
    assert all(abs(g.coeffs[0] - 1) < 1e-12 for g in enc.majoranas)


def test_jw_hamiltonian_is_ising_at_delta_equal_t():
    n, t, mu = 5, 1.0, 0.7
    h = hamiltonian(jordan_wigner(n), majorana_matrix(n, chain_edges(n), t, t, mu))
    ref = []
    for j in range(n - 1):
        ref.append(("XX", [j, j + 1], -t))
    for j in range(n):
        ref.append(("Z", [j], mu / 2))
    ref = SparsePauliOp.from_sparse_list(ref, n)
    assert np.allclose((h - ref).simplify(atol=1e-12).coeffs, 0)


def test_parity_hamiltonian_is_kramers_wannier_dual():
    n, t, mu = 5, 1.0, 0.7
    h = hamiltonian(parity(n), majorana_matrix(n, chain_edges(n), t, t, mu))
    ref = [("X", [j], -t) for j in range(n - 1)] + [("Z", [0], mu / 2)]
    ref += [("ZZ", [j - 1, j], mu / 2) for j in range(1, n)]
    ref = SparsePauliOp.from_sparse_list(ref, n)
    assert np.allclose((h - ref).simplify(atol=1e-12).coeffs, 0)


def _explicit_jw_matrix(n, edges, t, delta, mu):
    """Independent construction from a_j matrices built with Kronecker products."""
    I, Z = np.eye(2), np.diag([1.0, -1.0])
    lower = np.array([[0, 1], [0, 0]], dtype=complex)  # |0><1| : removes a particle
    def a(j):
        m = np.array([[1.0 + 0j]])
        for q in range(n):
            m = np.kron(m, Z if q < j else (lower if q == j else I))
        return m
    ops = [a(j) for j in range(n)]
    h = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for i, j in edges:
        hop = ops[i].conj().T @ ops[j]
        pair = ops[i] @ ops[j]
        h += -t * (hop + hop.conj().T) + delta * (pair + pair.conj().T)
    for j in range(n):
        h += -mu * (ops[j].conj().T @ ops[j] - 0.5 * np.eye(2 ** n))
    return h


@pytest.mark.parametrize("geom", ["chain", "ladder", "tjunction"])
def test_majorana_form_matches_second_quantised_definition(geom):
    if geom == "chain":
        n, edges = 6, chain_edges(6)
    elif geom == "ladder":
        n, edges = ladder_edges(3)
    else:
        n, edges = tjunction_edges(2)
    t, d, mu = 1.0, 0.6, 0.9
    ref = np.linalg.eigvalsh(_explicit_jw_matrix(n, edges, t, d, mu))
    amat = majorana_matrix(n, edges, t, d, mu)
    for name in ENCODINGS:
        got = np.linalg.eigvalsh(hamiltonian(get_encoding(name, n), amat).to_matrix())
        assert np.allclose(got, ref, atol=1e-10), name


def test_ed_matches_bdg():
    n = 8
    for mu in (0.3, 1.0, 2.5):
        amat = majorana_matrix(n, chain_edges(n), 1.0, 0.8, mu)
        lv = lowest_levels(hamiltonian(bravyi_kitaev(n), amat), k=2)
        assert abs(lv[0] - bdg_ground_energy(amat)) < 1e-10
        assert abs((lv[1] - lv[0]) - bdg_levels(amat)[0]) < 1e-10


def test_edge_correlator_ed_vs_bdg():
    n, mu = 6, 0.8
    amat = majorana_matrix(n, chain_edges(n), 1.0, 1.0, mu)
    for name in ENCODINGS:
        enc = get_encoding(name, n)
        _, v = ground_states(hamiltonian(enc, amat), 1)
        op = (1j * (enc.gamma(0) @ enc.gamma(2 * n - 1))).simplify().to_matrix()
        got = (v[:, 0].conj() @ op @ v[:, 0]).real
        assert abs(got - bdg_edge_correlator(amat)) < 1e-9, name


def test_parity_operator():
    n = 6
    assert jordan_wigner(n).parity_operator().paulis[0].to_label() == "Z" * n
    assert parity(n).parity_operator().paulis[0].to_label()[::-1] == "IIIIIZ"
    for name in ENCODINGS:
        enc = get_encoding(name, n)
        p = enc.parity_operator()
        h = hamiltonian(enc, majorana_matrix(n, chain_edges(n), 1, 0.5, 0.4))
        assert np.allclose(((h @ p) - (p @ h)).simplify(atol=1e-12).coeffs, 0)


def test_decode_roundtrip():
    enc = jordan_wigner(5)
    assert enc.decode("IIZII") == [4, 5]          # Z_2 = -i g4 g5 : local, even
    assert enc.decode("XIIII") == [0]             # the left edge Majorana
    assert len(enc.decode("IIXII")) % 2 == 1      # parity-odd, carries a string
    par = parity(5)
    assert par.decode("IXIII") == [3, 4]          # the Kitaev bond itself
