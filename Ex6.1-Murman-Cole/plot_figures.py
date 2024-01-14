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
    
    
    kdx = 1 #skip info line
    for idx in range(0,Nx):
        for jdx in range(0,Ny):
            curLine = lines[kdx].split()
            x[idx,jdx] = float(curLine[0])
            y[idx,jdx] = float(curLine[1])
            phi[idx,jdx] = float(curLine[2])
            
            kdx+=1
 
    return x,y,phi
    
def readXDataFile(filename):
    with open(filename) as f:
        lines = f.readlines()
    N = len(lines)
    x = np.empty([N]);
    y = np.empty([N]);
     
    for line,idx in zip(lines,range(0,N)):
        splLine = line.split()
        x[idx] = float(splLine[0])
        y[idx] = float(splLine[1])
 
    return x, y
    
####################################################################
#Figure 6.11
####################################################################
files = ["resid-C1-M0p735.dat","resid-C2-M0p908.dat"]
plt.figure()

plt.title("Figure 6.11 - Maximum Residual after 400 Iterations")
plt.yscale("log")
plt.xlabel("Steps")
plt.ylabel("Residual")
for file in files:
    iter,res = readXDataFile(file)
    plt.plot(iter,res)
plt.legend([f"$M_\infty=0.735$","$M_\infty=0.908$"])
plt.savefig("Fig-6.11-resid.png")

####################################################################
#Figure 6.7 TODO needs work, matplotlib not perfect for this
####################################################################
plt.figure()

plt.title(f"Figure 6.7 Mach Contours and streamlines of subsonic flow ($M_\infty=0.735$) ")
x, y, M = readXYDataFile("resultC3-M-M0p735.dat")
x, y, u = readXYDataFile("resultC3-u-M0p735.dat")
x, y, v = readXYDataFile("resultC3-v-M0p735.dat")

cvar = plt.contour(x,y,M,[0.7,0.75,0.8])
#plt.gca().clabel(cvar,fontsize=10) #<-not working wtf
#plt.streamplot(np.transpose(x),np.transpose(y),np.transpose(u),np.transpose(v)) requires equal spacing wtf
#plt.quiver(x,y,u,v,scale=0.02)

plt.xlim([-0.5,1.5])
plt.ylim([0,1.2])

plt.savefig("Fig-6.7-subsonic-streamlines.png")


####################################################################
#Fig 6.8 
####################################################################
plt.figure()
plt.title(f"Figure 6.8 - Subsonic Pressure Distribution($M_\infty=0.735$)")
plt.xlabel("x")
plt.ylabel("-Cp")
x,Cp = readXDataFile("resultC3-Cp-M0p735.dat")
plt.xlim([-0.5,1.5])
plt.plot(x,-Cp,"-o")
plt.savefig("Fig-6.8-subsonicCp.png")

####################################################################
#Figure 6.9 TODO needs work, matplotlib not perfect for this
####################################################################
plt.figure()

plt.title(f"Figure 6.9 Mach Contours and streamlines of transonic flow ($M_\infty=0.0.908$) ")
x, y, M = readXYDataFile("resultC4-M-M0p908.dat")
x, y, u = readXYDataFile("resultC4-u-M0p908.dat")
x, y, v = readXYDataFile("resultC4-v-M0p908.dat")

cvar = plt.contour(x,y,M,[1,1.1,1.2,1.3])
#plt.gca().clabel(cvar,fontsize=10) #<-not working wtf
#plt.streamplot(np.transpose(x),np.transpose(y),np.transpose(u),np.transpose(v)) requires equal spacing wtf
#plt.quiver(x,y,u,v,scale=0.02)

plt.xlim([-0.5,1.5])
plt.ylim([0,1.7])

plt.savefig("Fig-6.9-transonic-streamlines.png")

####################################################################
#Fig 6.10 
####################################################################
plt.figure()
plt.title(f"Figure 6.10 - Transonic Pressure Distribution($M_\infty=0.908$)")
plt.xlabel("x")
plt.ylabel("-Cp")
x,Cp = readXDataFile("resultC4-Cp-M0p908.dat")
plt.xlim([-0.5,1.5])
plt.plot(x,-Cp,"-o")
plt.savefig("Fig-6.10-transonicCp.png")

####################################################################
#Figure 6.12 TODO needs work, matplotlib not perfect for this
####################################################################
plt.figure()

plt.title(f"Figure 6.9 Mach Contours and streamlines of transonic flow ($M_\infty=0.0.908$) \n With Refined grid, 101x51")
x, y, M = readXYDataFile("resultC5-v-M0p908.dat")
x, y, u = readXYDataFile("resultC5-u-M0p908.dat")
x, y, v = readXYDataFile("resultC5-M-M0p908.dat")

cvar = plt.contour(x,y,M,[1,1.1,1.2,1.3])
#plt.gca().clabel(cvar,fontsize=10) #<-not working wtf
#plt.streamplot(np.transpose(x),np.transpose(y),np.transpose(u),np.transpose(v)) requires equal spacing wtf
#plt.quiver(x,y,u,v,scale=0.02)

plt.xlim([-0.5,1.5])
plt.ylim([0,1.7])

plt.savefig("Fig-6.12-transonic-streamlines-fine.png")

####################################################################
#Fig 6.13
####################################################################
plt.figure()
plt.title(f"Figure 6.13 - Transonic Pressure Distribution($M_\infty=0.908$)\n With Refined grid, 101x51")
plt.xlabel("x")
plt.ylabel("-Cp")
x,Cp = readXDataFile("resultC5-Cp-M0p908.dat")
plt.xlim([-0.5,1.5])
plt.plot(x,-Cp,"-o")
plt.savefig("Fig-6.13-transonicCp-fine.png")





#if len(sys.argv)>1:
#    datafiles.append(sys.argv[1])
#else:
#    for file in os.listdir():
#        if ("phi-" in file) and (".dat" in file):
#            datafiles.append(file)
#            print(file)
#
#for file in datafiles:
#
#    x,y, phi= readXYDataFile(file)
#
#    fig = plt.figure(figsize=[8,5])
#    ax1 = plt.subplot(1,2,1)
#    plt.contourf(x,y,phi,100)
#    plt.title(f"$\phi$ of flow field")
#    plt.colorbar()
#
#    ax2 = plt.subplot(1,2,2)
#    plt.xlim([-0.5,1.5])
#    plt.ylim([0.0,0.5])
#    #which = np.empty(x.shape,dtype=bool);
#    #which[:,:] = True
#    #which[x<-0.5] = False
#    #which[x>1.5] = False
#    #which[y<0.0] = False
#    #which[y>0.5] = False
#    
#    #print(min(phi[which]),max(phi[which]))
#    pos = plt.contourf(x,y,phi,200)
#    fig.colorbar(pos)
#    plt.title(f"$u$ near airfoil")
#    
#    fig.tight_layout()   
#    #plt.show()
#    fig.suptitle(file)
#    fname = file.replace(".dat",".png")
#    plt.savefig(fname)




