from ply import lex,yacc
import numpy as np 
from typing import Optional
from cirq import ops,Circuit,NamedQubit,CX
from typing import Any,Callable,cast,Dict,Iterable,List,Optional,Sequence,Union
import re

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

class OPENQASMLexer:
    def __init__(self):
        self.lex=lex.lex(object=self,debug=False)

    literals=['{','}','[',']','(',')',';',',','+','/','*','-','^']

    reserved={
        'qubit':'QUBIT',
        'bit':'BIT',
        'measure':'MEASURE',
        '->':'ARROW'
    }

    tokens=[
        'FORMAT_SPEC',
        'NUMBER',
        'NATURAL_NUMBER',
        'STDGATES',
        'ID',
        'PI',
    ]   + list(reserved.values())

    t_ignore=" \t"

    def t_COMMENT(self,t):
        r"""//.*"""
        print(t.value)
        # print("Comment encountered")
        pass

    def t_FORMAT_SPEC(self,t):
        r"""OPENQASM(\s+)([^\s\t\;]*);"""
        match = re.match(r"""OPENQASM(\s+)([^\s\t;]*);""", t.value)
        t.value = match.groups()[1]
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
        print("include file encountered")
        return t
    
    def t_QUBIT(self,t):
        r"""qubit"""
        return t

    def t_ID(self,t):
        r"""[a-zA-Z][a-zA-Z\d_]*"""
        return t

    def t_error(self, t):
        print("Illegal character '%s'" % t.value[0])
        t.lexer.skip(1)

    def token(self) -> Optional[lex.Token]:
        return self.lex.token()


class OPENQASMParser():
    def __init__(self,file):
        self.file=file
        self.lexer=OPENQASMLexer()
        self.parser = yacc.yacc(module=self, debug=True, write_tables=False)
        self.supported_format = False
        self.stdgatesinc = False
        self.qubits : Dict[str:int]={}
        self.bits : Dict[str:int]={}
        self.circuit= Circuit()
        with open(file) as f:
            self.data=f.read()
        self.lexer.lex.input(self.data)
        self.parsedQasm=self.parser.parse(lexer=self.lexer)
    
    tokens=OPENQASMLexer.tokens
    start='start'

    def p_start(self,p):
        """start : qasm"""
        p[0]=p[1]
    
    def p_qasm_format_only(self, p):
        """qasm : format"""
        self.supported_format = True
        p[0] = Qasm(self.supported_format, self.stdgatesinc, self.qubits, self.bits, self.circuit)
        
    def p_format(self, p):
        """format : FORMAT_SPEC"""
        if p[1] != "3":
            print("Oh my angel")
        else:
            print("Oh my eve")

    def p_include(self,p):
        """qasm : qasm STDGATES"""
        self.stdgatesinc=True
        p[0] = Qasm(self.supported_format, self.stdgatesinc, self.qubits, self.bits, self.circuit)

    def p_qasm_circuit(self,p):
        """qasm : qasm circuit"""
        p[0]=Qasm(self.supported_format,self.stdgatesinc,self.qubits,self.bits,p[2])

    def p_circuit_reg(self,p):
        """circuit : new_reg circuit"""
        p[0]=self.circuit
    

    def p_new_reg(self,p):
        """new_reg : QUBIT '[' NATURAL_NUMBER ']' ID ';'
        | BIT '[' NATURAL_NUMBER ']' ';'""" 
        name,length=p[5],p[3]
        print(name,length)
        if p[1] == "qubit":
            self.qubits[name] = length
        else:
            self.bits[name] = length
        p[0]=(name,length)
    
    def p_circuit_empty(self, p):
        """circuit : empty"""
        p[0] = self.circuit

    def p_qasm_empty(self, p):
        """qasm : empty"""
        print("Hey i am empty")

    def p_empty(self,p):
        """empty :"""
        pass

    def p_error(self, p):
        # print(p.value)
        print("Error encountered")
        if p is None :
            print("p is none type")
        else:
            print(p)
# --------------------------------------------------------------------------------------------------------------

parser=OPENQASMParser("demo3.qasm")
value=parser.parser.parse(lexer=parser.lexer)
print(value)

