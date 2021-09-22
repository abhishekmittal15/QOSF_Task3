# from qiskit import QuantumCircuit, Aer , execute 

# file=input("Enter the file name : ")
# file+=".qasm"

# qc=QuantumCircuit.from_qasm_file(file)

# print(qc)

from cirq.contrib.qasm_import import circuit_from_qasm 

circuit=circuit_from_qasm("""
OPENQASM 2.0;
include "qelib1.inc";
qreg q[1];
creg c[1];
x q;
measure q -> c;
""")

print(circuit)