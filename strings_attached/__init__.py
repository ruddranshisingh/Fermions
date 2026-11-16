"""Strings Attached: fermion-to-qubit encodings vs. simulated Majorana modes.

SC4055 group project code. See README.md for how to reproduce every figure.
"""
from .encodings import Encoding, jordan_wigner, parity, bravyi_kitaev, ternary_tree, get_encoding, ENCODINGS
from .kitaev import (chain_edges, ladder_edges, tjunction_edges, majorana_matrix,
                     hamiltonian, bdg_levels, bdg_ground_energy, bdg_edge_correlator,
                     lowest_levels, ground_states)
