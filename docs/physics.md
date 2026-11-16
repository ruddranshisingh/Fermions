# Physics notes

## 1. Kitaev chain

H = sum_j [ -t (a_j^dag a_{j+1} + h.c.) + Delta (a_j a_{j+1} + h.c.) ] - mu sum_j (n_j - 1/2)

Split each fermion into two Majoranas, a_j = (g_{2j} + i g_{2j+1}) / 2, with {g_a, g_b} = 2 delta_ab. Then

* n_j - 1/2 = (i/2) g_{2j} g_{2j+1}
* bond j -> j+1 = (i/2) [ (Delta - t) g_{2j} g_{2j+3} + (Delta + t) g_{2j+1} g_{2j+2} ]

At Delta = t, mu = 0 only the pairs (g_{2j+1}, g_{2j+2}) are coupled. g_0 and g_{2N-1} drop out of H: two zero modes, one non-local fermion f = (g_0 + i g_{2N-1})/2, two degenerate ground states of opposite fermion parity. For 0 < |mu| < 2|t| the modes are dressed and decay into the bulk as (mu/2t)^x; their overlap gives a splitting ~ (mu/2t)^N.

Any quadratic H is (i/4) sum A_ab g_a g_b with A real antisymmetric. The eigenvalues of iA are +-eps_k and H = sum_k eps_k (f_k^dag f_k - 1/2) (Bogoliubov-de Gennes). `kitaev.majorana_matrix` builds A; both the qubit Hamiltonian and the exact solution are derived from it.

## 2. Encodings

An encoding is a choice of 2N anticommuting Pauli strings. For linear encodings the qubits store b = A n (mod 2):

| | A | stores | g_{2j} | g_{2j+1} |
|---|---|---|---|---|
| Jordan-Wigner | identity | occupations | Z_0..Z_{j-1} X_j | Z_0..Z_{j-1} Y_j |
| Parity | lower-triangular ones | running parities | Z_{j-1} X_j X_{j+1}..X_{N-1} | Y_j X_{j+1}..X_{N-1} |
| Bravyi-Kitaev | Fenwick tree | partial sums | weight O(log N) | weight O(log N) |

General rule (Seeley-Richard-Love): g_{2j} = X_U Z_P, g_{2j+1} = i X_U Z_P Z_F with U = column j of A (qubits that flip), P = qubits giving the parity of modes < j, F = qubits giving n_j. The ternary tree uses the 2N+1 root-to-leaf paths of a ternary tree as Majoranas (weight ~ log_3 2N).

All encodings are unitarily equivalent (the Clifford algebra with 2N generators has a unique irreducible representation of dimension 2^N), so spectra agree. E1 tests this.

## 3. Same Hamiltonian, two faces (Delta = t)

* Jordan-Wigner: H = -t sum X_j X_{j+1} + (mu/2) sum Z_j, the transverse-field Ising model. Topological phase = ferromagnet. Ground states at mu = 0: (|+..+> +- |-..->)/sqrt 2, GHZ states.
* Parity: H = -t sum_{j<N-1} X_j + (mu/2) sum Z_{j-1} Z_j (Z_{-1} = 1), the Kramers-Wannier dual. Ground states at mu = 0: |+..+>|0> and |+..+>|1>, product states. Qubit N-1 does not appear in H: it *is* the fermion parity.

## 4. Why errors differ

Topological protection covers perturbations that are (i) local in the fermions and (ii) fermion-parity-even. A qubit Pauli error is some product of Majoranas (`Encoding.decode`), and whether it satisfies (i) and (ii) depends on the encoding.

* JW Z_j = -i g_{2j} g_{2j+1}: local, even. Protected.
* JW X_j = Z-string x g_{2j}: odd. In Ising language a longitudinal field: splits the ferromagnet by 2 eps N m.
* JW Y_j: odd, but harmless. The rotation U = prod_j exp(-i theta_j X_j / 2) commutes with XX and maps (mu/2) Z_j + eps Y_j to sqrt((mu/2)^2 + eps^2) Z_j. More generally H is real (class BDI) and Y is imaginary, so <psi|Y|psi> = 0 in any real eigenstate.
* Parity X_j (j < N-1) = -i g_{2j+1} g_{2j+2}: the bond itself. Protected.
* Parity Z_{N-1} = total parity: energy difference 2 eps between the two ground states. Y_{N-1} = g_{2N-1}, the bare right edge Majorana, with overlap sqrt(1 - (mu/2t)^2) on the dressed zero mode.
* Parity Z_j = (total parity) x prod_{k>j} (1 - 2 n_k): even but non-local. First-order effect 2 eps (-mu/2t)^(N-1-j).

## 5. Fermion parity as an error detector

P = prod_j (-i g_{2j} g_{2j+1}) commutes with H. Under JW it is Z on every qubit, so any single X or Y error anticommutes with it and is detected by post-selection. Under the parity encoding P = Z_{N-1} and only errors on that one qubit are detected. E4 measures the consequence.

## References

* A. Kitaev, Unpaired Majorana fermions in quantum wires, Phys.-Usp. 44, 131 (2001).
* J. T. Seeley, M. J. Richard, P. J. Love, The Bravyi-Kitaev transformation for quantum computation of electronic structure, J. Chem. Phys. 137, 224109 (2012).
* S. Bravyi, A. Kitaev, Fermionic quantum computation, Ann. Phys. 298, 210 (2002).
* Z. Jiang et al., Optimal fermion-to-qubit mapping via ternary trees, Quantum 4, 276 (2020).
* M. Greiter, V. Schnells, R. Thomale, The 1D Ising model and the topological phase of the Kitaev chain, Ann. Phys. 351, 1026 (2014).
* J. Alicea et al., Non-Abelian statistics and topological quantum information processing in 1D wire networks, Nat. Phys. 7, 412 (2011).
* X. Mi et al., Noise-resilient edge modes on a chain of superconducting qubits, Science 378, 785 (2022).
* F. Zatelli et al., Majorana parity qubit in coupled minimal Kitaev chains, arXiv:2607.09511 (2026).
