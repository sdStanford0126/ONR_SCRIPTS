import numpy as np
import matplotlib.pyplot as plt
from Universal_Subroutines import setPlotpref
import sys
import os
setPlotpref()
"""
the goal of this script is to utilize the axial profile probes to do the following:
1. plot the average profiles (vel, TKE, p, T) of the jet as it evolves
2. plot the spectrum as various locations to examine the turbulence development
    - at the centerline
    - at the lip lines
"""
def readPos(posFname:str):
    #read data tile 
    #posFname: full directory reference to the position file
    Pos = np.loadtxt(posFname,skiprows=1)
    ind=Pos[:,0]
    ind = ind.astype(int)
    Xr  =Pos[:,1]
    Yr  =Pos[:,2]
    Zr  =Pos[:,3]
    Npts = np.size(ind)
    X  =np.zeros(Npts)
    Y  =np.zeros(Npts)
    Z  =np.zeros(Npts)
    X[ind]  =Xr
    Y[ind]  =Yr
    Z[ind]  =Zr
    
    return ind,X,Y,Z,Npts

def readData(ind,DataName_fmt,tids,Nx,Ny,Nz):
    #inputs: 
    #ind: index to align data to positioning
    #DataName_fmt: format for the data
    #tids: array of time steps to be read
    #N{x,y,z}: number of points perlocation
    #outputs:
    #Var_t: time record of Var
    #Var list from file is u,v,w,p,T
    Nt = tids.size
    Npts = Nx*Ny*Nz
    u_t = np.zeros((Nt,Nx,Ny,Nz)) 
    v_t = np.zeros((Nt,Nx,Ny,Nz)) 
    w_t = np.zeros((Nt,Nx,Ny,Nz)) 
    p_t = np.zeros((Nt,Nx,Ny,Nz)) 
    T_t = np.zeros((Nt,Nx,Ny,Nz)) 

    u_r = np.zeros((Npts))
    v_r = np.zeros((Npts))
    w_r = np.zeros((Npts))
    p_r = np.zeros((Npts))
    T_r = np.zeros((Npts))

    for i,tid in enumerate(tids):
        DataName = DataName_fmt.format(tid)
        Data = np.loadtxt(DataName,skiprows=1)
        u_r[ind] = Data[:,0]    
        v_r[ind] = Data[:,1]    
        w_r[ind] = Data[:,2]    
        p_r[ind] = Data[:,3]    
        T_r[ind] = Data[:,4]
        u_t[i,:,:,:] = reformData(u_r,Nx,Ny,Nz)    
        v_t[i,:,:,:] = reformData(v_r,Nx,Ny,Nz)    
        w_t[i,:,:,:] = reformData(w_r,Nx,Ny,Nz)    
        p_t[i,:,:,:] = reformData(p_r,Nx,Ny,Nz)    
        T_t[i,:,:,:] = reformData(T_r,Nx,Ny,Nz)    
    
    return u_t, v_t,w_t,p_t,T_t

def reformData(Data,Nx,Ny,Nz):
    #restucture the geometry data
    #currently it is y first, x second, z last
    
    #validate data size
    if Data.size != Nx*Ny*Nz:
        raise ValueError("Data size does not match size from indicies")
    
    Data_f = np.zeros((Nx,Ny,Nz))
    for k in range(Nz):
        for i in range(Nx):
            ind_str = i*Ny+k*(Ny*Nx)
            ind_end = ind_str + Ny 
            Data_f[i,:,k] = Data[ind_str:ind_end]
        
    return Data_f

def plotSpec():
    pass

def main():
    #define input parameters 
    data_dir = "/anvil/scratch/x-sdai/AR2_base_151M_str/pcprobes"
    out_dir  = "/anvil/scratch/x-sdai/AR2_base_151M_post/ax_prof_post/"
    posFname = os.path.join(data_dir,"axial.pxyz")
    ind,X,Y,Z,Npts = readPos(posFname)
    print("total number of points is: ", Npts)
    Xu = np.unique(X)
    Yu = np.unique(Y)
    Zu = np.unique(Z)

    Nx = Xu.size
    Ny = Yu.size
    Nz = Zu.size
    print("Nx is %d, Ny is %d, Nz is %d" %(Nx,Ny,Nz)) 

    Xf = reformData(X,Nx,Ny,Nz)
    Yf = reformData(Y,Nx,Ny,Nz)
    Zf = reformData(Z,Nx,Ny,Nz)
    tid_str = 249400
    tid_end = 249800
    dt = 50
    tids = np.arange(tid_str,tid_end+dt,dt)
    print(tids[:51])

    DataName_fmt = "axial.{:08d}.pcd"
    DataName_fmt = os.path.join(data_dir,DataName_fmt)
    print("test output ", DataName_fmt.format(tid_str))

    u_t,v_t,w_t,p_t,T_t = readData(ind,DataName_fmt,tids,Nx,Ny,Nz)

    #test plot
    plt.figure()
    plt.pcolor(Yu,Zu,u_t[0,9,:,:])  
    plt.savefig("ax_u_test.png")  

    """
    #test reformData
    plt.figure()
    plt.plot(Xf[:,0,0])
    plt.savefig("ax_Xf.png")
    plt.figure()
    plt.plot(Yf[0,:,0])
    plt.savefig("ax_Yf.png")
    plt.figure()
    plt.plot(Zf[0,0,:])
    plt.savefig("ax_Zf.png")
    """    

    """
    #testing of rolling indices
    #plt.figure()
    #plt.plot(X[:1001],label="x")
    #plt.savefig(out_dir+"ax_X.png")
    #plt.figure()
    #plt.plot(Y[:1001],label="y")
    #plt.savefig(out_dir+"ax_Y.png")
    #plt.figure()
    #plt.plot(Z[:2*Nx*Ny],label="z")
    #plt.savefig(out_dir+"ax_Z.png")
    """ 
if __name__ == "__main__":
    main()
