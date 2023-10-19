import subprocess
s=["GLY", "SER"]
t=["298","328","358"]
n="50000"
wb="True"

for ss in s:
    for tt in t:
        for ii in range(10):
            subprocess.call(["python", "1FirstExample.py", "-s", ss, "-n", n, "-t", tt, "--wb", wb, "-o", ss+"_"+str(ii+1)+"_"+tt+"K_wb"])