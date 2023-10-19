import pandas as pd
import dynamics as dyn


import sys
import ast

if len(sys.argv) != 8:
    print("Usage: python script.py arg1 arg2 arg3 arg4 arg5 arg6 arg7")
else:
    namepdb = sys.argv[1]
    sequence = sys.argv[2]
    temp = sys.argv[3]
    steps = sys.argv[4]
    timestep = sys.argv[5]
    n = sys.argv[6]
    i = sys.argv[7]
print(namepdb)

a = int(temp)
t= []
t.append(a)
steps= int(steps)
timestep= float(timestep)
n= int(n)
i= int(i)
dyn.mdcalc(namepdb, sequence, t, steps, timestep, n, i)


'''      
if len(sys.argv) != 8:
    print("Usage: python script.py arg1 arg2 arg3 arg4 arg5 arg6 arg7")
else:
    namepdb = sys.argv[1]
    sequence = sys.argv[2]
    temp = sys.argv[3]
    steps = sys.argv[4]
    timestep = sys.argv[5]
    n = sys.argv[6]
    i = sys.argv[7]
print(namepdb)
dyn.mdcalc(namepdb, sequence, temp, steps, timestep, n, i)
dyn.mdcalc(arg1, arg2, arg7, arg3, arg4, arg5, arg6)
'''
