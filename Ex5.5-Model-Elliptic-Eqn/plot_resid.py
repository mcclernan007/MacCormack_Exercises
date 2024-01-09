import sys
import os
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
    if "resid.dat" in file:
        datafiles.append(file)

legendEntries = []
plt.figure()
plt.yscale("log")
for file in datafiles:
    x,Cp= readXYDataFile(file)
    plt.plot(x,Cp)
    legendEntries.append(file)
plt.legend(legendEntries)
plt.show()
