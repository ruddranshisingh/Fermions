"""Circuit cost of one first-order Trotter step, after routing to hardware."""
from __future__ import annotations

from typing import Dict, Optional

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import SparsePauliOp
from qiskit.synthesis import LieTrotter
from qiskit.transpiler import CouplingMap

BASIS = ["cx", "rz", "sx", "x"]


def coupling_map(kind: str, n: int) -> Optional[CouplingMap]:
    """Smallest device of each family that fits n qubits."""
    if kind == "all_to_all":
        return None
    if kind == "line":
        return CouplingMap.from_line(n)
    if kind == "grid":
        rows = int(np.floor(np.sqrt(n)))
        cols = int(np.ceil(n / rows))
        return CouplingMap.from_grid(rows, cols)
    if kind == "heavy_hex":
        d = 3
        while CouplingMap.from_heavy_hex(d).size() < n:
            d += 2
        return CouplingMap.from_heavy_hex(d)
    raise ValueError(kind)


def trotter_step(h: SparsePauliOp, dt: float = 0.1, reps: int = 1) -> QuantumCircuit:
    qc = QuantumCircuit(h.num_qubits)
    qc.append(PauliEvolutionGate(h, time=dt, synthesis=LieTrotter(reps=reps)), range(h.num_qubits))
    return qc


def pauli_weights(h: SparsePauliOp) -> Dict:
    w = [sum(ch != "I" for ch in p.to_label()) for p in h.paulis]
    return {"terms": len(w), "max_weight": int(max(w)), "mean_weight": float(np.mean(w)),
            "total_weight": int(sum(w))}


def routed_cost(h: SparsePauliOp, kind: str, seed: int = 0, optimization_level: int = 2) -> Dict:
    tqc = transpile(trotter_step(h), coupling_map=coupling_map(kind, h.num_qubits), basis_gates=BASIS,
                    optimization_level=optimization_level, seed_transpiler=seed)
    ops = tqc.count_ops()
    return {"cx": int(ops.get("cx", 0)), "depth": int(tqc.depth()),
            "cx_depth": int(tqc.depth(lambda inst: inst.operation.num_qubits == 2)),
            "one_qubit": int(sum(v for k, v in ops.items() if k in ("sx", "x", "rz")))}
