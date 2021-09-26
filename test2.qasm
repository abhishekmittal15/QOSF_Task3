OPENQASM 3;
include "stdgates.inc";
qubit[2] q;
x q[0];
h q[1];
cx q[0],q[1];

