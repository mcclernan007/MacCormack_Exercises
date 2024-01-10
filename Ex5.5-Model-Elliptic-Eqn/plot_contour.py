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
    Nx = int(infospl[1])
    Ny = int(infospl[2])
    x = np.empty([Nx,Ny]);
    y = np.empty([Nx,Ny]);
    
    phi = np.empty([Nx,Ny]);
    
    u = np.empty([Nx,Ny]);
    v = np.empty([Nx,Ny]);
    
    kdx = 1 #skip info line
    for idx in range(0,Nx):
        for jdx in range(0,Ny):
            curLine = lines[kdx].split()
            x[idx,jdx] = float(curLine[0])
            y[idx,jdx] = float(curLine[1])
            phi[idx,jdx] = float(curLine[2])
            u[idx,jdx] = float(curLine[3])
            v[idx,jdx] = float(curLine[4])
            
            kdx+=1
 
    return x,y,phi,u,v
    
    
datafiles = []    
if len(sys.argv)>1:
    datafiles.append(sys.argv[1])
else:
    for file in os.listdir():
        if ("phiUV-" in file) and (".dat" in file):
            datafiles.append(file)
            print(file)

for file in datafiles:

    x,y, phi,u,v= readXYDataFile(file)

    fig = plt.figure(figsize=[8,5])
    ax1 = plt.subplot(1,2,1)
    plt.contourf(x,y,phi,100)
    plt.title(f"$\phi$ of flow field")
    plt.colorbar()

    ax2 = plt.subplot(1,2,2)
    plt.xlim([-0.5,1.5])
    plt.ylim([0.0,0.5])
    #which = np.empty(x.shape,dtype=bool);
    #which[:,:] = True
    #which[x<-0.5] = False
    #which[x>1.5] = False
    #which[y<0.0] = False
    #which[y>0.5] = False
    
    #print(min(phi[which]),max(phi[which]))
    pos = plt.contourf(x,y,u,200)
    fig.colorbar(pos)
    plt.title(f"$u$ near airfoil")
    
    fig.tight_layout()   
    #plt.show()
    fig.suptitle(file)
    fname = file.replace(".dat",".png")
    plt.savefig(fname)




