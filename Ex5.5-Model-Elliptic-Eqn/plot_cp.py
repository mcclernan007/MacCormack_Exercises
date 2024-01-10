import sys
import matplotlib.pyplot as plt
import numpy as np
import os

      
#    x = 
    
def readXYDataFile(filename):
    with open(filename) as f:
        lines = f.readlines()
    info = lines[0]
    infospl = info.split()
    N = int(infospl[0])
    x = np.empty([N]);
    Cp = np.empty([N]);
     
    kdx = 1 #skip info line
    for idx in range(0,N):
        curLine = lines[kdx].split()
        x[idx] = float(curLine[0])
        Cp[idx] = float(curLine[1])
            
        kdx+=1
 
    return x,Cp

datafiles = []    
if len(sys.argv)>1:
    datafiles.append(sys.argv[1])
else:
    for file in os.listdir():
        if ("cp-" in file) and (".dat" in file):
            datafiles.append(file)
            print(file)

legendEntries = []
legendLines = []
plt.figure()

for file in datafiles:
    x,Cp = readXYDataFile(file)
    plt.plot(x,Cp,"-o")
    legendEntries.append(file)
plt.legend(legendEntries)
plt.xlim([-0.5,1.5])
plt.gca().invert_yaxis()
plt.ylabel("Cp (inverted axis)")
plt.xlabel("x")
#plt.show()
plt.savefig("Cp.png")
