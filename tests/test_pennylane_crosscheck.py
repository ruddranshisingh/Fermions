"""Independent check: PennyLane's own JW / parity / BK transforms must give
the same spectrum as our encodings.  Skipped if PennyLane is not installed."""
import numpy as np
import pytest

qml = pytest.importorskip("pennylane")

from strings_attached import chain_edges, get_encoding, hamiltonian, majorana_matrix, tjunction_edges


def _fermi_hamiltonian(n, edges, t, delta, mu):
    from pennylane.fermi import FermiWord, FermiSentence
    terms = {}
    def add(word, c):
        terms[word] = terms.get(word, 0) + c
    for i, j in edges:
        add(FermiWord({(0, i): "+", (1, j): "-"}), -t)
        add(FermiWord({(0, j): "+", (1, i): "-"}), -t)
        add(FermiWord({(0, i): "-", (1, j): "-"}), delta)
        add(FermiWord({(0, j): "+", (1, i): "+"}), delta)
    for j in range(n):
        add(FermiWord({(0, j): "+", (1, j): "-"}), -mu)
    add(FermiWord({}), mu * n / 2)
    return FermiSentence(terms)


@pytest.mark.parametrize("geom", ["chain", "tjunction"])
def test_spectrum_matches_pennylane(geom):
    n, edges = (6, chain_edges(6)) if geom == "chain" else tjunction_edges(2)
    t, delta, mu = 1.0, 0.6, 0.9
    fs = _fermi_hamiltonian(n, edges, t, delta, mu)
    amat = majorana_matrix(n, edges, t, delta, mu)
    transforms = {
        "jordan_wigner": lambda: qml.jordan_wigner(fs),
        "parity": lambda: qml.parity_transform(fs, n),
        "bravyi_kitaev": lambda: qml.bravyi_kitaev(fs, n),
    }
    for name, make in transforms.items():
        ref = np.linalg.eigvalsh(qml.matrix(make(), wire_order=range(n)))
        got = np.linalg.eigvalsh(hamiltonian(get_encoding(name, n), amat).to_matrix())
        assert np.allclose(got, ref, atol=1e-9), name
