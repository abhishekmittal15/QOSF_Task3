from ply import lex, yacc
import numpy as np 
import cirq
from typing import Optional
from cirq import ops, Circuit, NamedQubit, CX
from typing import Any, Callable, cast, Dict, Iterable, List, Optional, Sequence, Union
import re
from cirq.circuits.qasm_output import QasmUGate
from cirq.contrib.qasm_import._lexer import QasmLexer
from cirq.contrib.qasm_import.exception import QasmException
from qasm import Qasm, QasmGateStatement

# Lexer describing OpenQASM 3
class OPENQASMLexer:
    def __init__(self):
        self.lex = lex.lex(object=self,debug=False)

    literals = ['{','}','[',']','(',')',';',',','+','/','*','-','^']

    reserved = {
        'qubit':'QUBIT',
        'bit':'BIT',
        'measure':'MEASURE',
        '->':'ARROW',
        ':=' : 'ASSIGN'
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
    
    def t_ASSIGN(self, t):
        r""":="""
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

    def t_ARROW(self, t):
        """->"""
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
        """circuit :  circuit gate_op"""
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
    #     p[0] = self.circuit
    
    # def p_variable_assign(self,p):
    #     """variable_assign : ID ASSIGN NUMBER ';' variable_assign"""
    #     name, data = p[1], p[3]
    #     self.variables[name]=data
    #     print(f"{self.variables} has been added to variables")

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

    def p_variable_assign_empty(self, p):
        """variable_assign : empty"""
        print("Empty variable assign")

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
    print("\nRequired Circuit: ")
    circuit_properties=parser.parsedQasm
    num_qubits=list(circuit_properties.qubits.values())[0]
    circuit=circuit_properties.circuit
    print(circuit)

    print(herm_circuit(num_qubits,circuit))