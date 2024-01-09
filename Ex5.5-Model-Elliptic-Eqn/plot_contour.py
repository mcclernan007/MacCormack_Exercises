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
    
xydatafile = sys.argv[1];
x,y, phi,u,v= readXYDataFile(xydatafile)

#print(x[:,1])
#plt.figure()
#plt.plot(x[:,1],v[:,1])
#plt.scatter(x[:,:],u[:,:])
#plt.plot(x[:,2],v[:,2],"--")
#plt.xlim([-0.5,1.5])
#plt.show()

#plotGrid(xNodes,yNodes);
plt.figure()
#plt.scatter(x,y,20,phi)

#plt.ylim([0,0.1])
#plt.xlim([-0.5,1.5])
plt.contourf(x,y,u,100)
#plt.quiver(x,y,u*10000,v*10000)

plt.colorbar()

plt.show()


