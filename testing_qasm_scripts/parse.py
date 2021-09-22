file=input("Enter the file name")
file+=".qasm"

with open(file) as f:
    str=f.read()

str=str.split("\n")

specification_declared=False

for i in str:
    # remove trailing white spaces 
    line=i.split(" ")
    # Check for comments 
    if(line[0]=="//"):
        continue
    # Check for specification OPENQASM 3.0 
    if(i=="OPENQASM 3;"):
        specification_declared=True
    if(specification_declared==False):
        print("Specify the specification first")
        exit()
    # Check for include files 
    if(line[0]=="include"):
        file=line[1].split('"')[1]
        print(file)
    # Check for qubit 
    if(line[0].startswith("qubit")):
        print(line[0].split("["))
        qubit_name=line[1]
        print(qubit_name)
    # Check for bit 
    # Check for in built gates ( x,h,cx)

