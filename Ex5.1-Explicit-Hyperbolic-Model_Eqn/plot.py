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
    
def readDataFile(filename):
    with open(filename) as f:
        lines = f.readlines()
    
    lines = lines[1:-1] #rm number points (first entry)
    print(lines)
    x = np.empty(len(lines));
    u_IC = np.empty(len(lines));
    u = np.empty(len(lines));
    for line, idx in zip(lines,list(range(0,len(lines)))) :
        curLine = line.split()
        x[idx] = float(curLine[0])
        u_IC[idx] = float(curLine[1])
        u[idx] = float(curLine[2])
    return x,u_IC,u
    
datafile = sys.argv[1];
x,u_IC,u = readDataFile(datafile)
#print([x,u_IC,u])
plt.plot(x,u_IC,color='#101010',linestyle='dashed')
plt.plot(x,u)
plt.show()
#plotGrid(xNodes,yNodes);
