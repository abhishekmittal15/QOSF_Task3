from ply import lex, yacc
import numpy as np 
import cirq
import functools
from typing import Optional
from cirq import ops, Circuit, NamedQubit, CX
from typing import Any, Callable, cast, Dict, Iterable, List, Optional, Sequence, Union
import re
from cirq.circuits.qasm_output import QasmUGate
from cirq.contrib.qasm_import.exception import QasmException

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

# Lexer describing OpenQASM 3
class OPENQASMLexer:
    def __init__(self):
        self.lex = lex.lex(object=self,debug=False)

    literals = ['{','}','[',']','(',')',';',',','+','/','*','-','^','=']

    reserved = {
        'qubit':'QUBIT',
        'bit':'BIT',
        'measure':'MEASURE',
        '=':'ASSIGN'
    }

    tokens = [
        'FORMAT_SPEC',
        'NUMBER',
        'NATURAL_NUMBER',
        'STDGATES',
        'ID',
        'PI',
    ]   + list(reserved.values())

    t_ignore = " \t"

    def t_COMMENT(self,t):
        r"""//.*"""
        print(t.value)
        pass

    def t_FORMAT_SPEC(self,t):
        r"""OPENQASM(\s+)([^\s\t\;]*);"""
        match = re.match(r"""OPENQASM(\s+)([^\s\t;]*);""", t.value)
        t.value = match.groups()[1]
        return t
    
    def t_NUMBER(self, t):
        r"""(
        (
        [0-9]+\.?|
        [0-9]?\.[0-9]+
        )
        [eE][+-]?[0-9]+
        )|
        (
        ([0-9]+)?\.[0-9]+|
        [0-9]+\.)"""
        t.value = float(t.value)
        return t

    def t_NATURAL_NUMBER(self, t):
        r"""\d+"""
        t.value = int(t.value)
        return t

    def t_newline(self,t):
        r'\n+'
        t.lexer.lineno += len(t.value)
        
    def t_STDGATES(self,t):
        r"""include(\s+)"stdgates.inc";"""
        print("\ninclude file encountered")
        return t
    
    def t_QUBIT(self,t):
        r"""qubit"""
        return t
    
    def t_BIT(self,t):
        r"""bit"""
        return t

    def t_ID(self,t):
        r"""[a-zA-Z][a-zA-Z\d_]*"""
        return t

    def t_MEASURE(self, t):
        r"""measure"""
        return t

    def t_ASSIGN(self, t):
        """="""
        return t
    
    def t_PI(self, t):
        r"""pi"""
        t.value = np.pi
        return 

    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)

    def token(self) -> Optional[lex.Token]:
        return self.lex.token()

# This is the main parser function, which converts
# OpenQASM 3.0 code to a cirq Circuit()
class OPENQASMParser():
    def __init__(self,file):
        # file containing OpenQASM code
        self.file = file
        self.lexer = OPENQASMLexer()
        self.parser = yacc.yacc(module=self, debug=True, write_tables=False)
        # parameters to be passed to the QASM class
        self.supported_format = False
        self.stdgatesinc = False
        self.qubits : Dict[str:int] = {}
        self.bits : Dict[str:int] = {}
        self.circuit = Circuit()
        self.variables : Dict[str:float] = {}
        # qubits are modelled as Qid
        self.qids: Dict[str, ops.Qid] = {}
        # gates described are X, Y, Z, RX, RY, RZ, H, S, S†, T, T†, CX, CCX, SWAP 
        # & CSWAP, as required
        self.gates = {
            'x': QasmGateStatement(qasm_gate='x', num_params=0, num_args=1, cirq_gate=ops.X),
            'y': QasmGateStatement(qasm_gate='y', num_params=0, num_args=1, cirq_gate=ops.Y),
            'z': QasmGateStatement(qasm_gate='z', num_params=0, num_args=1, cirq_gate=ops.Z),
            'h': QasmGateStatement(qasm_gate='h', num_params=0, num_args=1, cirq_gate=ops.H),
            's': QasmGateStatement(qasm_gate='s', num_params=0, num_args=1, cirq_gate=ops.S),
            't': QasmGateStatement(qasm_gate='t', num_params=0, num_args=1, cirq_gate=ops.T),  
            'cx': QasmGateStatement(qasm_gate='cx', num_params=0, num_args=2, cirq_gate=CX),
            'ccx': QasmGateStatement(qasm_gate='ccx', num_params=0, num_args=3, cirq_gate=ops.CCX),
            'swap': QasmGateStatement(qasm_gate='swap', num_params=0, num_args=2, cirq_gate=ops.SWAP),
            'cswap': QasmGateStatement(qasm_gate='cswap', num_params=0, num_args=3, cirq_gate=ops.CSWAP),
            'sdg': QasmGateStatement(qasm_gate='sdg', num_params=0, num_args=1, cirq_gate=ops.S ** -1),
            'tdg': QasmGateStatement(qasm_gate='tdg', num_params=0, num_args=1, cirq_gate=ops.T ** -1),
            'rx': QasmGateStatement(
                qasm_gate='rx', cirq_gate=(lambda params: ops.rx(params[0])), num_params=1, num_args=1
            ),
            'ry': QasmGateStatement(
                qasm_gate='ry', cirq_gate=(lambda params: ops.ry(params[0])), num_params=1, num_args=1
            ),
            'rz': QasmGateStatement(
                qasm_gate='rz', cirq_gate=(lambda params: ops.rz(params[0])), num_params=1, num_args=1
            ),
        }
        with open(file) as f:
            self.data = f.read()
            
        print("This is the code")
        print(self.data)
        print("")
        self.lexer.lex.input(self.data)
        self.parsedQasm = self.parser.parse(lexer=self.lexer)
    
    tokens = OPENQASMLexer.tokens
    start = 'start'

    # start : qasm
    def p_start(self,p):
        """start : qasm"""
        p[0] = p[1] 
    
    # qasm : qasm STDGATES
    #      | format
    #      | qasm circuit
    #      | empty
    def p_qasm_format_only(self, p):
        """qasm : format"""
        self.supported_format = True
        p[0] = Qasm(self.supported_format, self.stdgatesinc, self.qubits, self.bits, self.circuit)
        
    # format : FORMAT_SPEC
    def p_format(self, p):
        """format : FORMAT_SPEC"""
        if p[1] != "3":
            print("Incorrect format specified, ... overriding")
        else:
            print("GG! Correct format specified")

    def p_include(self,p):
        """qasm : qasm STDGATES"""
        self.stdgatesinc=True
        p[0] = Qasm(self.supported_format, self.stdgatesinc, self.qubits, self.bits, self.circuit)

    def p_qasm_circuit(self,p):
        """qasm : qasm circuit"""
        p[0] = Qasm(self.supported_format,self.stdgatesinc,self.qubits,self.bits,p[2])

    # circuit : new_reg circuit
    #         | gate_op circuit
    #         | measurement circuit
    #         | variable_assign circuit
    #         | empty
    def p_circuit_reg(self,p):
        """circuit : new_reg circuit"""
        p[0] = self.circuit
    
    def p_circuit_gate_or_measurement(self, p):
        """circuit :  circuit gate_op
        | circuit measurement"""
        self.circuit.append(p[2])
        p[0] = self.circuit

    def p_new_reg(self,p):
        """new_reg : QUBIT '[' NATURAL_NUMBER ']' ID ';'
        | BIT '[' NATURAL_NUMBER ']' ID ';'""" 
        name, length = p[5], p[3]
        if p[1] == "qubit":
            print(f'Qubit(s) with name: "{name}" and length "{length}"')
            self.qubits[name] = length
        else:
            print(f'Bit(s) with name: "{name}" and length "{length}"')
            self.bits[name] = length
        p[0] = (name,length)

    # def p_qasm_variable(self, p):
    #     """circuit : variable_assign circuit"""
    #     p[0] =  Qasm(self.supported_format,self.stdgatesinc,self.qubits,self.bits,p[2])
    
    # def p_variable_assign(self,p):
    #     """variable_assign : ID ASSIGN NUMBER ';'"""
    #     name, data = p[1], p[3]
    #     self.variables[name]=data
    #     print(f"{self.variables}")

    # gate_op : ID qargs
    #         | ID ( params ) qargs
    def p_gate_op_no_params(self, p):
        """gate_op :  ID qargs"""
        print(p[2], p[1])
        self.gate_operation(p[2], gate=p[1], p=p, params=[])
    
    def p_gate_op_with_params(self, p):
        """gate_op :  ID '(' params ')' qargs"""
        # print(p[3], p[1] , p[2])
        self.gate_operation(args=p[5], gate=p[1], p=p, params=p[3])
    
    def gate_operation(
        self, args: List[List[ops.Qid]], gate: str, p: Any, params: List[float]
    ):
        # assuming that stdgates has been included
        gate_set = self.gates
        if gate not in gate_set.keys():
            print("Unknown Gate: ", gate)
        p[0] = gate_set[gate].on(args=args, params=params, lineno=p.lineno(1))

    def p_measurement(self,p):
        """measurement : carg ASSIGN MEASURE qarg ';'"""
        qubits=p[4]
        bits=p[1]
        p[0]=[
            ops.MeasurementGate(num_qubits=1,key=bits[i]).on(qubits[i]) for i in range(len(qubits))
        ]

    def p_classical_arg_bit(self, p):
        """carg : ID '[' NATURAL_NUMBER ']' """
        reg = p[1]
        idx = p[3]
        arg_name = self.make_name(idx, reg)
        if reg not in self.bits.keys():
            raise QasmException(f'Undefined classical register "{reg}" at line {p.lineno(1)}')

        size = self.bits[reg]
        if idx >= size:
            raise QasmException(
                'Out of bounds bit index {} '
                'on classical register {} of size {} '
                'at line {}'.format(idx, reg, size, p.lineno(1))
            )
        p[0] = [arg_name]
    
    def p_classical_arg_register(self, p):
        """carg : ID """
        reg = p[1]
        if reg not in self.bits.keys():
            raise QasmException(f'Undefined classical register "{reg}" at line {p.lineno(1)}')

        p[0] = [self.make_name(idx, reg) for idx in range(self.bits[reg])]


    # qargs : qarg ',' qargs
    #      | qarg ';'
    def p_args_multiple(self, p):
        """qargs : qarg ',' qargs"""
        p[3].insert(0, p[1])
        p[0] = p[3]

    def p_args_single(self, p):
        """qargs : qarg ';'"""
        p[0] = [p[1]]

    def make_name(self, idx, reg):
        return str(reg) + "_" + str(idx)

    # qarg : ID
    #     | ID '[' NATURAL_NUMBER ']'
    def p_quantum_arg_register(self, p):
        """qarg : ID """
        reg = p[1]
        if reg not in self.qubits.keys():
            raise QasmException(f'Undefined quantum register "{reg}" at line {p.lineno(1)}')
        qids = []
        for idx in range(self.qubits[reg]):
            arg_name = self.make_name(idx, reg)
            if arg_name not in self.qids.keys():
                self.qids[arg_name] = NamedQubit(arg_name)
            qids.append(self.qids[arg_name])
        p[0] = qids
    
    def p_quantum_arg_bit(self, p):
        """qarg : ID '[' NATURAL_NUMBER ']' """
        reg = p[1]
        idx = p[3]
        arg_name = self.make_name(idx, reg)
        if reg not in self.qubits.keys():
            raise QasmException(f'Undefined quantum register "{reg}" at line {p.lineno(1)}')
        size = self.qubits[reg]
        if idx >= size:
            raise QasmException(
                'Out of bounds qubit index {} '
                'on register {} of size {} '
                'at line {}'.format(idx, reg, size, p.lineno(1))
            )
        if arg_name not in self.qids.keys():
            self.qids[arg_name] = NamedQubit(arg_name)
        p[0] = [self.qids[arg_name]]

    # params : param ',' params
    #        | param
    def p_params_multiple(self, p):
        """params : expr ',' params"""
        p[3].insert(0, p[1])
        p[0] = p[3]

    def p_params_single(self, p):
        """params : expr """
        p[0] = [p[1]]
    
    # expr : term
    #      | func '(' expression ')' """
    #      | binary_op
    #      | unary_op

    def p_expr_term(self, p):
        """expr : term"""
        p[0] = p[1]

    def p_expr_parens(self, p):
        """expr : '(' expr ')'"""
        p[0] = p[2]

    def p_expr_function_call(self, p):
        """expr : ID '(' expr ')'"""
        func = p[1]
        if func not in self.functions.keys():
            raise QasmException(f"Function not recognized: '{func}' at line {p.lineno(1)}")
        p[0] = self.functions[func](p[3])

    def p_expr_unary(self, p):
        """expr : '-' expr
        | '+' expr"""
        if p[1] == '-':
            p[0] = -p[2]
        else:
            p[0] = p[2]

    def p_expr_binary(self, p):
        """expr : expr '*' expr
        | expr '/' expr
        | expr '+' expr
        | expr '-' expr
        | expr '^' expr
        """
        p[0] = self.binary_operators[p[2]](p[1], p[3])

    def p_term(self, p):
        """term : NUMBER
        | NATURAL_NUMBER
        | PI"""
        p[0] = p[1]

    def p_circuit_empty(self, p):
        """circuit : empty"""
        p[0] = self.circuit

    # def p_variable_assign_empty(self, p):
    #     """variable_assign : empty"""
    #     print("Empty variable assign")

    def p_qasm_empty(self, p):
        """qasm : empty"""
        print("This is empty")

    # Acts as a sink
    def p_empty(self,p):
        """empty :"""
        pass

    # Basic error flagging
    def p_error(self, p):
        print("Error encountered")
        if p is None:
            print("p is none type")
        exit()

class mygate(cirq.Gate):
    def __init__(self,qubits_num,unitary):
        super(mygate,self)
        self.qubits_num=qubits_num
        self.unitary=unitary

    def _num_qubits_(self):
        return self.qubits_num

    def _unitary_(self):
        return self.unitary

    def _circuit_diagram_info_(self,args):
        return "G1" , "G2" , "G3"

def herm_circuit(qubits,circuit):
    unitary=cirq.unitary(circuit)
    hermitian=unitary.conjugate()
    hermitian=hermitian.transpose()

    hermitian_gate=mygate(qubits,hermitian)
    return cirq.Circuit(hermitian_gate.on(*cirq.LineQid.for_gate(hermitian_gate)))

if __name__ == "__main__":
    parser = OPENQASMParser("demo.qasm")
    circuit_properties=parser.parsedQasm
    num_qubits=list(circuit_properties.qubits.values())[0]
    circuit=circuit_properties.circuit
    
    print("\nRequired Circuit: ")
    print(circuit)
    print("\n Hermitian Circuit")
    print(herm_circuit(num_qubits,circuit))