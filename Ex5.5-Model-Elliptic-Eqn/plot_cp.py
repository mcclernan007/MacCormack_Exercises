import sys
import matplotlib.pyplot as plt
import numpy as np


def plotNodes(xdata,ydata):
    plt.figure()
    plt.grid()
    plt.scatter(xdata,ydata)
    plt.show()
def plotGrid(xNodes,yNodes):
    #pretty limited to this type of grid, but good enough for now
    #Very lazy implementation
    myColor = "#0072BD"
    myWidth = 0.5
    xRaw = np.sort(np.unique(xNodes));
    yRaw = np.sort(np.unique(yNodes));
    
    plt.figure()
    #hlines
    for idx in list(range(1,len(xRaw))):
        for jdx in list(range(0,len(yRaw))):
            plt.plot(xRaw[idx-1:idx+1],yRaw[jdx]*np.array([1,1]),color=myColor,linewidth=myWidth)
    #vlines
    for idx in list(range(0,len(xRaw))):
        for jdx in list(range(1,len(yRaw))):    
            plt.plot(xRaw[idx]*np.array([1,1]),yRaw[jdx-1:jdx+1],color=myColor,linewidth=myWidth)
    #plt.scatter(xNodes,yNodes)
    
    plt.show()        
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
