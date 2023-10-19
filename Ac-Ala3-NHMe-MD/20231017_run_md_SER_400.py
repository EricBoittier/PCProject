
import pandas as pd
import dynamics as dyn
import subprocess
#Parameters
#sequence="GLY"
sequence="SER"

#temp=[298]
#temp=[328]
#temp=[358]
temp=[400]

#number of steps
steps = 10**5

#timestep
timestep = 0.01

#number of samples
n=10 

i=0
for _ in dyn.get_list_of_pdb(sequence, temp, n):
    i += 1
    namepdb = _
    t=str(temp[0])
    subprocess.call(["python", "20231017_md.py", 
                     namepdb, 
                     sequence,
                     t,
                     str(steps),
                     str(timestep), 
                     str(n), 
                     str(i)])
                     
                     
'''                    
                     "-s", ss, "-n", n, "-t", tt, "-o", ss+"_"+str(ii+1)+"_"+tt+"K_gas"])
'''
