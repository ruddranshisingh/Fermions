# The Physics Behind The Functions 


## 1. Kitaev chain

$$
H=\sum_j\Big[-t\,\big(a_j^\dagger a_{j+1}+\text{h.c.}\big)+\Delta\,\big(a_j a_{j+1}+\text{h.c.}\big)\Big]-\mu\sum_j\Big(n_j-\tfrac12\Big)
$$

Split each fermion into two Majoranas:

$$
a_j=\tfrac12\big(\gamma_{2j}+i\,\gamma_{2j+1}\big),\qquad \gamma_a\gamma_b+\gamma_b\gamma_a=2\,\delta_{ab}
$$

Then

$$
n_j-\tfrac12=\tfrac{i}{2}\,\gamma_{2j}\gamma_{2j+1}
$$

$$
\text{bond } j\to j+1:\quad \tfrac{i}{2}\Big[(\Delta-t)\,\gamma_{2j}\gamma_{2j+3}+(\Delta+t)\,\gamma_{2j+1}\gamma_{2j+2}\Big]
$$

At $\Delta=t$, $\mu=0$ only the pairs $(\gamma_{2j+1},\gamma_{2j+2})$ are coupled. $\gamma_0$ and $\gamma_{2N-1}$ drop out of $H$: two zero modes, one non-local fermion

$$
f=\tfrac12\big(\gamma_0+i\,\gamma_{2N-1}\big),
$$

and two degenerate ground states of opposite fermion parity. For $0<\lvert\mu\rvert<2\lvert t\rvert$ the modes are dressed and decay into the bulk as $(\mu/2t)^x$; their overlap gives a splitting

$$
E_1-E_0\sim\Big(\frac{\mu}{2t}\Big)^{N}.
$$

Any quadratic Hamiltonian can be written

$$
H=\frac{i}{4}\sum_{a,b}A_{ab}\,\gamma_a\gamma_b,\qquad A \text{ real antisymmetric}.
$$

The eigenvalues of $iA$ are $\pm\varepsilon_k$ and

$$
H=\sum_k\varepsilon_k\Big(f_k^\dagger f_k-\tfrac12\Big)
$$

(Bogoliubov-de Gennes). `kitaev.majorana_matrix` builds $A$; both the qubit Hamiltonian and the exact solution are derived from it.

## 2. Encodings

An encoding is a choice of $2N$ anticommuting Pauli strings. For linear encodings the qubits store $b=An \pmod 2$:

| | $A$ | stores | $\gamma_{2j}$ | $\gamma_{2j+1}$ |
|---|---|---|---|---|
| Jordan-Wigner | identity | occupations | $Z_0\cdots Z_{j-1}X_j$ | $Z_0\cdots Z_{j-1}Y_j$ |
| Parity | lower-triangular ones | running parities | $Z_{j-1}X_jX_{j+1}\cdots X_{N-1}$ | $Y_jX_{j+1}\cdots X_{N-1}$ |
| Bravyi-Kitaev | Fenwick tree | partial sums | weight $O(\log N)$ | weight $O(\log N)$ |

General rule (Seeley-Richard-Love):

$$
\gamma_{2j}=X_{U(j)}\,Z_{P(j)},\qquad \gamma_{2j+1}=i\,X_{U(j)}\,Z_{P(j)}\,Z_{F(j)}
$$

with $U(j)$ = column $j$ of $A$ (qubits that flip), $P(j)$ = qubits giving the parity of modes $<j$, and $F(j)$ = qubits giving $n_j$. The ternary tree uses the $2N+1$ root-to-leaf paths of a ternary tree as Majoranas (weight $\sim\log_3 2N$).

All encodings are unitarily equivalent (the Clifford algebra with $2N$ generators has a unique irreducible representation of dimension $2^N$), so spectra agree. E1 tests this.

## 3. Same Hamiltonian, two faces ($\Delta=t$)

**Jordan-Wigner** gives the transverse-field Ising model:

$$
H_{\text{JW}}=-t\sum_{j}X_jX_{j+1}+\frac{\mu}{2}\sum_j Z_j
$$

The topological phase is the ferromagnet. Ground states at $\mu=0$ are GHZ states:

$$
\frac{1}{\sqrt2}\Big(\lvert{+}\cdots{+}\rangle\pm\lvert{-}\cdots{-}\rangle\Big)
$$

**Parity** gives its Kramers-Wannier dual (with $Z_{-1}\equiv 1$):

$$
H_{\text{par}}=-t\sum_{j=0}^{N-2}X_j+\frac{\mu}{2}\sum_{j=0}^{N-1}Z_{j-1}Z_j
$$

Ground states at $\mu=0$ are product states, $\lvert{+}\cdots{+}\rangle\lvert0\rangle$ and $\lvert{+}\cdots{+}\rangle\lvert1\rangle$. Qubit $N-1$ does not appear in $H$: it *is* the fermion parity.

## 4. Why errors differ

Topological protection covers perturbations that are (i) local in the fermions and (ii) fermion-parity-even. A qubit Pauli error is some product of Majoranas (`Encoding.decode`), and whether it satisfies (i) and (ii) depends on the encoding.

**Jordan-Wigner $Z_j$**: local and even, so protected.

$$
Z_j=-i\,\gamma_{2j}\gamma_{2j+1}
$$

**Jordan-Wigner $X_j$**: a $Z$-string times $\gamma_{2j}$, so odd. In Ising language it is a longitudinal field, which splits the ferromagnet by

$$
E_1-E_0\approx 2\,\varepsilon\,N\,m .
$$

**Jordan-Wigner $Y_j$**: odd, but harmless. The rotation

$$
U=\prod_j\exp\Big(-\tfrac{i}{2}\,\theta_j X_j\Big)
$$

commutes with every $X_jX_{j+1}$ and maps

$$
\frac{\mu}{2}Z_j+\varepsilon\,Y_j\;\longmapsto\;\sqrt{\Big(\frac{\mu}{2}\Big)^2+\varepsilon^2}\;Z_j ,
$$

so a $Y$ error is exactly a slightly larger chemical potential. More generally $H$ is real (class BDI) and $Y$ is imaginary, so $\langle\psi\rvert Y\lvert\psi\rangle=0$ in any real eigenstate.

**Parity $X_j$ with $j<N-1$**: the Kitaev bond itself, so protected.

$$
X_j=-i\,\gamma_{2j+1}\gamma_{2j+2}
$$

**Parity, qubit $N-1$**: $Z_{N-1}$ is the total parity, giving an energy difference $2\varepsilon$ between the two ground states. $Y_{N-1}=\gamma_{2N-1}$ is the bare right edge Majorana, with overlap $\sqrt{1-(\mu/2t)^2}$ on the dressed zero mode.

**Parity $Z_j$ with $j<N-1$**: even but non-local.

$$
Z_j=P\prod_{k>j}\big(1-2n_k\big),\qquad \text{first-order effect}=2\,\varepsilon\Big(-\frac{\mu}{2t}\Big)^{N-1-j}
$$

## 5. Fermion parity as an error detector

$$
P=\prod_j\big(-i\,\gamma_{2j}\gamma_{2j+1}\big)
$$

commutes with $H$. Under Jordan-Wigner $P=Z_0Z_1\cdots Z_{N-1}$, so any single $X$ or $Y$ error anticommutes with it and is detected by post-selection. Under the parity encoding $P=Z_{N-1}$ and only errors on that one qubit are detected. E4 measures the consequence.

## References

* A. Kitaev, Unpaired Majorana fermions in quantum wires, Phys.-Usp. 44, 131 (2001).
* J. T. Seeley, M. J. Richard, P. J. Love, The Bravyi-Kitaev transformation for quantum computation of electronic structure, J. Chem. Phys. 137, 224109 (2012).
* S. Bravyi, A. Kitaev, Fermionic quantum computation, Ann. Phys. 298, 210 (2002).
* Z. Jiang et al., Optimal fermion-to-qubit mapping via ternary trees, Quantum 4, 276 (2020).
* M. Greiter, V. Schnells, R. Thomale, The 1D Ising model and the topological phase of the Kitaev chain, Ann. Phys. 351, 1026 (2014).
* J. Alicea et al., Non-Abelian statistics and topological quantum information processing in 1D wire networks, Nat. Phys. 7, 412 (2011).
* X. Mi et al., Noise-resilient edge modes on a chain of superconducting qubits, Science 378, 785 (2022).
* F. Zatelli et al., Majorana parity qubit in coupled minimal Kitaev chains, arXiv:2607.09511 (2026).
