import cirq 
from cirq import Simulator
sim=Simulator()
import numpy as np

a,b=cirq.LineQubit.range(2)

circuit=cirq.Circuit(cirq.X(a),cirq.H(b),cirq.CNOT(a,b))

print(circuit)
unitary=cirq.unitary(circuit)

hermitian=unitary.conjugate()
hermitian=hermitian.transpose()

print("Hermitian circuit")

class mygate(cirq.Gate):
    def __init__(self):
        super(mygate,self)

    def _num_qubits_(self):
        return 2

    def _unitary_(self):
        return hermitian

    def _circuit_diagram_info_(self,args):
        return "G1" , "G2"

my_gate=mygate()

# x,y=cirq.LineQubit.range(2)

hermitian_circuit=cirq.Circuit(my_gate.on(*cirq.LineQid.for_gate(my_gate)))

print(hermitian_circuit)

result=sim.simulate(hermitian_circuit)

print(np.around(result.final_state_vector,3))
