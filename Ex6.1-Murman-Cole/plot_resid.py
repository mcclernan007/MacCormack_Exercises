import sys
import os
import matplotlib.pyplot as plt
import numpy as np
       
#    x = 
    
def readXYDataFile(filename):
    with open(filename) as f:
        lines = f.readlines()
    N = len(lines)
    iter = np.empty([N]);
    res = np.empty([N]);
     
    for line,idx in zip(lines,range(0,N)):
        splLine = line.split()
        iter[idx] = int(splLine[0])
        res[idx] = float(splLine[1])
 
    return iter,res
    
datafiles = []
for file in os.listdir():
    if ("resid" in file) and (".dat" in file):
        datafiles.append(file)
        print(file)

legendEntries = []
plt.figure()
plt.yscale("log")
for file in datafiles:
    x,Cp= readXYDataFile(file)
    plt.plot(x,Cp)
    legendEntries.append(file)
plt.legend(legendEntries)
plt.show()
#plt.savefig("resid.png")
