OPENQASM 3;
include "stdgates.inc";

qubit[3] q;
bit[3] c;
h q[0];
cx q[0], q[1];
cx q[1], q[0];
sdg q[2];
t q[2];
rx(0.2) q[0];
// c = measure q;