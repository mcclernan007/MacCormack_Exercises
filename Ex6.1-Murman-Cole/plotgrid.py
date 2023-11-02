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
    #interior
    for idx in list(range(1,len(xRaw)-1)):
        for jdx in list(range(1,len(yRaw)-1)):
            plt.plot(xRaw[idx-1:idx+2],yRaw[jdx]*np.array([1,1,1]),color=myColor,linewidth=myWidth)
            plt.plot(xRaw[idx]*np.array([1,1,1]),yRaw[jdx-1:jdx+2],color=myColor,linewidth=myWidth)
    #edges
    for idx in list(range(1,len(xRaw))):
        plt.plot(xRaw[idx-1:idx+1],yRaw[0]*np.array([1,1]),color=myColor,linewidth=myWidth)
        plt.plot(xRaw[idx-1:idx+1],yRaw[-1]*np.array([1,1]),color=myColor,linewidth=myWidth)
    for jdx in list(range(1,len(xRaw))):
        plt.plot(xRaw[0]*np.array([1,1]),yRaw[jdx-1:jdx+1],color=myColor,linewidth=myWidth)
        plt.plot(xRaw[-1]*np.array([1,1]),yRaw[jdx-1:jdx+1],color=myColor,linewidth=myWidth)
                
    #plt.scatter(xNodes,yNodes)
    
    plt.show()        
#    x = 
    
def readXYDataFile(filename):
    with open(filename) as f:
        lines = f.readlines()
    x = np.empty(len(lines));
    y = np.empty(len(lines));
    for line, idx in zip(lines,list(range(0,len(lines)))) :
        curLine = line.split()
        x[idx] = float(curLine[0])
        y[idx] = float(curLine[1])
    return x,y
    
xydatafile = sys.argv[1];
xNodes,yNodes = readXYDataFile(xydatafile)
plotGrid(xNodes,yNodes);
