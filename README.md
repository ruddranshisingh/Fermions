# Strings Attached

**How the fermion-to-qubit encoding decides the cost, and the topological protection, of simulated Majorana modes.**

SC4055 Introduction to Quantum Computing (NTU, 26S1) group project.

A Kitaev chain stores one bit in two Majorana zero modes at its ends, and that bit is protected against *local, fermion-parity-preserving* noise. On a qubit processor the chain has to be written in qubit language first (Jordan-Wigner, parity, Bravyi-Kitaev, ternary tree). The physics is identical in every encoding, but the *hardware errors* are not: a single-qubit Pauli error can be a harmless local term in one encoding and a logical error in another. This repository measures that, and what each encoding costs in gates.

| | Question | Method | Headline result |
|---|---|---|---|
| **E1** | Do the encodings preserve the physics? | Exact diagonalisation vs exact free-fermion (BdG) solution, N = 12 | All four agree to 3 x 10^-13 |
| **E2** | Which qubit errors does topology protect against? | Ground-state splitting vs N under static Pauli fields; automatic Pauli-to-Majorana decoder | Jordan-Wigner: X errors split the states ~2 eps N, Y and Z are protected. Parity: everything is protected except errors near the last qubit |
| **E3** | What does a Trotter step cost on hardware? | Qiskit transpiler, 4 connectivities, 3 geometries, 5 seeds | 1D chain: Jordan-Wigner 22 CNOTs, Bravyi-Kitaev 74, ternary tree 144 |
| **E4** | Does fermion-parity post-selection rescue a noisy run? | Density-matrix simulation with depolarising noise (Qiskit Aer) | The ranking flips when post-selection is on: Jordan-Wigner catches 58% of the damage, parity encoding 12% |

Full numbers and discussion: [`RESULTS.md`](RESULTS.md). Physics background: [`docs/physics.md`](docs/physics.md).

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

pytest -q                            # 28 checks, ~2 s
python scripts/run_experiments.py    # E1-E4 -> results/*.json   (~3 min)
python scripts/make_figures.py       # -> figures/*.png
python scripts/build_demo.py         # -> demo/index.html
open demo/index.html                 # interactive demo, no server needed
```

`python scripts/run_experiments.py --quick` runs smaller sizes in under a minute; `python scripts/run_experiments.py e2` runs one experiment.

## Layout

```
strings_attached/
  encodings.py   Majorana operators of each encoding; Pauli -> Majorana decoder
  kitaev.py      Kitaev model on any graph (chain, ladder, T-junction); BdG; exact diagonalisation
  errors.py      error dictionary: decoder, first-order map, splitting under perturbations
  cost.py        Trotter step + routing cost on line / heavy-hex / grid / all-to-all
  noise.py       noisy density-matrix evolution + parity post-selection
scripts/         run_experiments.py, make_figures.py, build_demo.py
tests/           unit tests, including an independent cross-check against PennyLane
results/         JSON output of every experiment (committed, so figures are reproducible)
figures/         report figures
demo/            interactive demo page (index.html is self-contained)
docs/            physics notes
```

## How the code is checked

* Every encoding's 2N Majorana operators are verified to anticommute.
* The Jordan-Wigner Hamiltonian is compared term by term with the transverse-field Ising model, and the parity-encoded one with its Kramers-Wannier dual.
* The Majorana-form Hamiltonian is compared against an independent construction from explicit a, a-dagger matrices on a chain, a ladder and a T-junction.
* Spectra are cross-checked against PennyLane's own `jordan_wigner`, `parity_transform` and `bravyi_kitaev`.
* Exact diagonalisation is compared with the free-fermion solution (energies, splitting, edge correlator).

## Tools

Qiskit 2.x (operators, transpiler), Qiskit Aer (density-matrix noise simulation), NumPy/SciPy, Matplotlib. PennyLane is optional and used only as an independent check. No real quantum hardware is needed.

## Limitations

* A single Kitaev chain cannot hold a parity superposition; a real qubit needs two chains (a Majorana box qubit).
* Sizes are limited to N <= 14 (exact diagonalisation) and N = 6 (density matrices).
* Noise is idealised: static Pauli fields in E2, gate depolarising noise without idle or readout errors in E4.
* E3 uses Qiskit's default Pauli-evolution synthesis; fermionic-swap networks or encoding-aware compilers would change absolute numbers.
