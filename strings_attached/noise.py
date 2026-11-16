"""Noisy Trotter evolution of the encoded ground state, with fermion-parity
post-selection (density-matrix simulation, Qiskit Aer)."""
from __future__ import annotations

from typing import Dict

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, depolarizing_error

from .cost import BASIS, trotter_step
from .encodings import Encoding


def even_parity_ground_state(enc: Encoding, h: SparsePauliOp) -> np.ndarray:
    """Exact ground state inside the even fermion-parity sector."""
    hm = h.to_matrix()
    pm = enc.parity_operator().to_matrix()
    proj = (np.eye(hm.shape[0]) + pm) / 2
    w, v = np.linalg.eigh(proj @ hm @ proj - 100.0 * proj)  # shift selects the sector
    return v[:, 0]


def noise_model(p2: float, p1_ratio: float = 0.1) -> NoiseModel:
    nm = NoiseModel(basis_gates=BASIS)
    nm.add_all_qubit_quantum_error(depolarizing_error(p2, 2), ["cx"])
    nm.add_all_qubit_quantum_error(depolarizing_error(p2 * p1_ratio, 1), ["sx", "x"])
    return nm  # rz is a virtual (frame-change) gate: error-free


def build_circuit(h: SparsePauliOp, steps: int, dt: float, line: bool = True, seed: int = 0) -> QuantumCircuit:
    n = h.num_qubits
    cm = CouplingMap.from_line(n) if line else None
    return transpile(trotter_step(h, dt * steps, reps=steps), coupling_map=cm, basis_gates=BASIS,
                     optimization_level=2, seed_transpiler=seed, initial_layout=list(range(n)))


def run(enc: Encoding, h: SparsePauliOp, tqc: QuantumCircuit, psi0: np.ndarray, p2: float) -> Dict:
    """Returns raw / post-selected fidelity (w.r.t. the noiseless circuit
    output), energy error, and the fraction of runs the parity check keeps."""
    n = enc.n

    def simulate(noise):
        qc = QuantumCircuit(n)
        qc.set_density_matrix(np.outer(psi0, psi0.conj()))
        qc.compose(tqc, inplace=True)
        qc.save_density_matrix()
        sim = AerSimulator(method="density_matrix", noise_model=noise)
        return np.asarray(sim.run(qc, shots=1).result().data(0)["density_matrix"])

    rho_ideal = simulate(None)
    rho = simulate(noise_model(p2)) if p2 > 0 else rho_ideal
    # operators expressed on the physical qubits at the END of the routed circuit
    par = enc.parity_operator().apply_layout(tqc.layout).to_matrix()
    ham = h.apply_layout(tqc.layout).to_matrix()
    proj = (np.eye(2 ** n) + par) / 2
    keep = float(np.trace(proj @ rho).real)
    rho_ps = proj @ rho @ proj / keep
    fid = float(np.trace(rho_ideal @ rho).real)          # rho_ideal is pure
    e_ideal = float(np.trace(rho_ideal @ ham).real)
    return {
        "fidelity_raw": fid,
        "fidelity_postselected": float(np.trace(rho_ideal @ rho_ps).real),
        "kept_fraction": keep,
        "energy_error_raw": float(np.trace(rho @ ham).real) - e_ideal,
        "energy_error_postselected": float(np.trace(rho_ps @ ham).real) - e_ideal,
        "ideal_parity": float(np.trace(rho_ideal @ par).real),
    }
