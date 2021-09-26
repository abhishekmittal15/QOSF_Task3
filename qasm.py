import functools
from ply import lex, yacc
import numpy as np 
from typing import Optional
from cirq import ops, Circuit, NamedQubit, CX
from typing import Any, Callable, cast, Dict, Iterable, List, Optional, Sequence, Union
import re
from cirq.circuits.qasm_output import QasmUGate
from cirq.contrib.qasm_import._lexer import QasmLexer
from cirq.contrib.qasm_import.exception import QasmException

class Qasm:
    """Qasm stores the final result of the Qasm parsing."""

    def __init__(
        self, supported_format: bool, stdgates: bool, qubits: dict, bits: dict, c: Circuit
    ):
        # defines whether the Quantum Experience standard header
        # is present or not
        self.stdgates = stdgates
        # defines if it has a supported format or not
        self.supportedFormat = supported_format
        # circuit
        self.qubits = qubits
        self.bits = bits
        self.circuit = c

class QasmGateStatement:
    """Specifies how to convert a call to an OpenQASM gate
    to a list of `cirq.GateOperation`s.
    Has the responsibility to validate the arguments
    and parameters of the call and to generate a list of corresponding
    `cirq.GateOperation`s in the `on` method.
    """

    def __init__(
        self,
        qasm_gate: str,
        cirq_gate: Union[ops.Gate, Callable[[List[float]], ops.Gate]],
        num_params: int,
        num_args: int,
    ):
        """Initializes a Qasm gate statement.
        Args:
            qasm_gate: the symbol of the QASM gate
            cirq_gate: the gate class on the cirq side
            num_args: the number of qubits (used in validation) this
                        gate takes
        """
        self.qasm_gate = qasm_gate
        self.cirq_gate = cirq_gate
        self.num_params = num_params

        # at least one quantum argument is mandatory for gates to act on
        assert num_args >= 1
        self.num_args = num_args

    # pylint: enable=missing-param-doc
    def _validate_args(self, args: List[List[ops.Qid]], lineno: int):
        if len(args) != self.num_args:
            raise QasmException(
                "{} only takes {} arg(s) (qubits and/or registers), "
                "got: {}, at line {}".format(self.qasm_gate, self.num_args, len(args), lineno)
            )

    def _validate_params(self, params: List[float], lineno: int):
        if len(params) != self.num_params:
            raise QasmException(
                "{} takes {} parameter(s), got: {}, at line {}".format(
                    self.qasm_gate, self.num_params, len(params), lineno
                )
            )

    def on(
        self, params: List[float], args: List[List[ops.Qid]], lineno: int
    ) -> Iterable[ops.Operation]:
        self._validate_args(args, lineno)
        self._validate_params(params, lineno)

        reg_sizes = np.unique([len(reg) for reg in args])
        if len(reg_sizes) > 2 or (len(reg_sizes) > 1 and reg_sizes[0] != 1):
            raise QasmException(
                f"Non matching quantum registers of length {reg_sizes} at line {lineno}"
            )

        # the actual gate we'll apply the arguments to might be a parameterized
        # or non-parameterized gate
        final_gate: ops.Gate = (
            self.cirq_gate if isinstance(self.cirq_gate, ops.Gate) else self.cirq_gate(params)
        )
        # OpenQASM gates can be applied on single qubits and qubit registers.
        # We represent single qubits as registers of size 1.
        # Based on the OpenQASM spec single qubit arguments can be mixed with qubit registers.
        # Given quantum registers of length reg_size and single qubits are both
        # used as arguments, we generate reg_size GateOperations via iterating
        # through each qubit of the registers 0 to n-1 and use the same one
        # qubit from the "single-qubit registers" for each operation.
        op_qubits = cast(Sequence[Sequence[ops.Qid]], functools.reduce(np.broadcast, args))
        for qubits in op_qubits:
            if isinstance(qubits, ops.Qid):
                yield final_gate.on(qubits)
            elif len(np.unique(qubits)) < len(qubits):
                raise QasmException(f"Overlapping qubits in arguments at line {lineno}")
            else:
                yield final_gate.on(*qubits)