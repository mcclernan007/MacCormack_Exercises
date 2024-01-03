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
    plt.plot(x,u_IC,color='#999999',linestyle='dashed')
    plt.plot(x,u_exact,color='#999999',linestyle='dashdot')
    plt.plot(x,u,marker='.')
    plt.xlabel(r"$x$",fontsize=fs)
    plt.ylabel(r"$u$",fontsize=fs)
    plt.title(outFilePath,loc="left",fontsize=fs)

outdir = "./fig5.2/"
subplotData = [
"implicit_central.dat",
"crank_nich.dat", 
"implicit_central_CFL2.dat",
"crank_nich_CFL2.dat"]
subplotTitles= [
"a) Method 4 -- Implicit Central, CFL=0.9",
"b) Method 5 -- Crank-Nicholson, CFL=0.9",
"c) Method 4 -- Implicit Central, CFL=2.0",
"d) Method 5 -- Crank-Nicholson, CFL=2.0"]

exactX,tmp,exactU1 = readDataFile("exact_soln1.dat")
exactX,tmp,exactU2 = readDataFile("exact_soln2.dat")


#fig = plt.figure(figsize=[5.5,6])

#fig = plt.figure(figsize=[5.5,6]) #viewable on my screen
#fs = 6
fig = plt.figure(figsize=[11,6]) #printing/saving
fs = 8


for idx in range(0,len(subplotData)):
    curDataFile = subplotData[idx]
    curSubplotTitle = subplotTitles[idx]
    curX,curU_IC,curU = readDataFile(curDataFile)
    ax = plt.subplot(2,2,idx+1)
    
    exactU = exactU1#exact solutions
    if idx in [2,3]: 
        exactU = exactU2
    
    plotResult(curX,curU_IC,curU,exactU,curSubplotTitle)
    
    plt.xlim([0,2])    
    plt.xticks(np.arange(0, 2.5, step=0.5),fontsize=fs)
    plt.ylim([0.25,1.25])
    plt.yticks(np.arange(0.25, 1.50, step=0.25),fontsize=fs)  # Set label locations.
        
    #IC label
    ax.text(0.05, 0.375, "Initial\nsolution",fontsize=fs)#,,xycoords='data')
    plt.plot([0.21, 0.5],[0.55,0.75],color='#999999')
    #exact soln label
    if idx in [0,1]:
        ax.text(1.1,0.8, "Exact\nsolution",fontsize=fs)
        plt.plot([0.98,1.08],[0.95,0.95],color='#999999') 
    else:
        ax.text(1.7,0.8, "Exact\nsolution",fontsize=fs)
        plt.plot([1.54,1.68],[0.95,0.95],color='#999999') 

fig.tight_layout()        
plt.show()
#plotGrid(xNodes,yNodes);

