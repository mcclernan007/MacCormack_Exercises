import sys
import os
import matplotlib.pyplot as plt
import numpy as np
 
def readDataFile(filename):
    with open(filename) as f:
        lines = f.readlines()
    
    lines = lines[1:-1] #rm number points (first entry)
    x = np.empty(len(lines));
    u_IC = np.empty(len(lines));
    u = np.empty(len(lines));
    for line, idx in zip(lines,list(range(0,len(lines)))) :
        curLine = line.split()
        x[idx] = float(curLine[0])
        u_IC[idx] = float(curLine[1])
        u[idx] = float(curLine[2])
    return x,u_IC,u
def getExactSolution(iwhich):
    if iwhich == 1:
        x,u_IC,u = readDataFile("exact_soln1.dat")
    elif iwhich == 2: 
        x,u_IC,u = readDataFile("exact_soln2.dat")
    else:
        x,u_IC,u = readDataFile("exact_soln1.dat")
    return u
def plotResult(x,u_IC,u,u_exact,outFilePath):
    plt.figure()
    plt.plot(x,u_IC,color='#999999',linestyle='dashed')
    plt.plot(x,u_exact,color='#999999',linestyle='dashdot')
    plt.plot(x,u,marker='.')
    plt.title(outFilePath)
    #plt.show()
    plt.savefig(outFilePath)

if len(sys.argv)>1: #specified individual
    datafile = sys.argv[1];
    x,u_IC,u = readDataFile(datafile)
#print([x,u_IC,u])
    
else: #plot everything and save
    files = os.listdir();
    for file in files:
        if ".dat" in file: 
            x,u_IC,u = readDataFile(file)
            if "CFL_2" in file:
                u_exact = getExactSolution(2)
            else:
                u_exact = getExactSolution(1)
            fnames = file.split(".");
            plotResult(x,u_IC,u,u_exact,fnames[0]+".png")
            
#plotGrid(xNodes,yNodes);
