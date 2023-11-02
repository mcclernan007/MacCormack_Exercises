import sys
import matplotlib.pyplot as plt
import numpy as np


def plotNodes(xdata,ydata):
    plt.figure()
    plt.grid()
    plt.scatter(xdata,ydata)
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
xdata,ydata = readXYDataFile(xydatafile)
plotNodes(xdata,ydata);
