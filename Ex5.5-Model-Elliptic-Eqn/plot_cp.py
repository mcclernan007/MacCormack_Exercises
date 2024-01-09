import sys
import matplotlib.pyplot as plt
import numpy as np

      
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
    
xydatafile = sys.argv[1];
x,Cp= readXYDataFile(xydatafile)


plt.figure()
plt.plot(x,Cp)
plt.scatter(x,Cp)
#plt.plot(x[:,2],v[:,2],"--")
#plt.xlim([-0.5,1.5])
plt.show()
