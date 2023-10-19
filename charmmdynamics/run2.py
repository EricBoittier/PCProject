import subprocess
s=["GLY"]
t=["328","358","400"]
n="50000"

for ss in s:
    for tt in t:
        for ii in range(10):
            subprocess.call(["python", "1FirstExample.py", "-s", ss, "-n", n, "-t", tt, "-o", ss+"_"+str(ii+1)+"_"+tt+"K_gas"])