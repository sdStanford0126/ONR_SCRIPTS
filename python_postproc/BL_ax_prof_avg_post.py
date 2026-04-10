#THis is being used to post processing time averaged
#results for the BL_interior_probes (axial profiles)
import os 
import numpy as np
import matplotlib.pyplot as plt

L = 2.56 #baselength
LEVEL=9 #mesh refinement level
delta = L/2**LEVEL

#Nz is in fact a function of x
x_str = -1.0152632
x_end = 0.0
z_str = 0.42498737
z_end = 0.5
slope = (z_end-z_str)/(x_end-x_str)

#old versoin
xs = np.linspace(-1,0,5)
#new version

y = np.arange(-(1.0 - delta / 2.0), (1.0 - delta / 2.0) + delta * 0.5, delta)




def z_lim(x):
    if x > 0:
        raise ValueError("x must be smaller or equal to 0")
    else:
        z_lim = z_str + slope*(x - x_str)
    return z_lim

def findNz(x,delta):
    z_max = z_lim(x)
    z = np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta)
    Nz = len(z)
    return Nz

Nx = xs.size
Ny = y.size
Nzs = [findNz(x,delta) for x in xs]

print("Nz_position by count",Nzs)
print("total probe pts:", np.sum(Ny*Nzs))

#TODO: need to implement extraction of each axial plane of probes, 
#first need to get index from .pxyz and record x,y,z locatiions (no read original probe file)
#then for a given x value, find the closest plane and extract the data from 
# that plane and organize them in the y,z plane

def extract_data(x,posName,DataName_fmt,tid_str,tid_end,dt):
    """
    INPUT: x: x position
           posName: .pxyz file name
           DataName_fmt: .pcd file name format (takes in INT tid to generate DataName)
           tid_str: starting tid
           tid_end: ending tid
           dt: step size between samples

    OUPUT:  
            ALL OUTPUTS ARE 3D ARRAYS, FIRST INDEX IN TIME, SECOND INDEX IN Y, THIRD INDEX IN Z
            u: x velocity component
            v: y velocity component
            w: z velocity component
            p: pressure
            rho: density
            mf: rho*u (mass flux)
    """ 
    #grab position data    
    Pos=np.loadtxt(posName,skiprows=1)
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
    """
    #print("data size", Npts)
    #print("expected data size", np.sum(np.fromiter((Ny*Nz for Nz in Nzs),int)))
    #plt.figure()
    #plt.plot(X)
    #plt.plot(Y)
    #plt.plot(Z)
    #plt.savefig("test.png")
    # grab time series data from the given x position
    """
    X_ind = np.where(np.abs(X-x) < 1e-6)
    Nz = findNz(x,delta)
    print(X_ind[0].size)
    print(Nz*Ny)
    #now read data
    tids = np.arange(tid_str,tid_end,dt)
    Nt = np.size(tids)
    u  = np.zeros((Nt,Ny,Nz))
    rho  = np.zeros((Nt,Ny,Nz))
    for i,tid in enumerate(tids):
        DataName = DataName_fmt.format(int(tid))
        print(DataName)
        DataR = np.loadtxt(DataName,skiprows=1)
        DataRo = np.zeros(np.shape(DataR))
        DataRo[ind,:]=DataR[:,:]
        u[i,:,:] = np.reshape(DataRo[X_ind[0],0],(Ny,Nz))
        v[i,:,:] = np.reshape(DataRo[X_ind[0],1],(Ny,Nz))
        w[i,:,:] = np.reshape(DataRo[X_ind[0],2],(Ny,Nz))
        p[i,:,:] = np.reshape(DataRo[X_ind[0],3],(Ny,Nz))
        rho[i,:,:] = np.reshape(DataRo[X_ind[0],4],(Ny,Nz))
        
    
    print(u.size)    
    print(Nt*Ny*Nz)
    mf = u*rho
    return u,v,w,p,rho,mf
   
    



if __name__ == "__main__":
    data_dir = "/anvil/scratch/x-sdai/BL_test_baseline_0.01875/pcprobe_int_axprof"
    posName = os.path.join(data_dir,"int_axprof.pxyz")
    fname_fmt = "int_axprof.{:08d}.pcd"
    DataName_fmt = os.path.join(data_dir,fname_fmt)
    print(xs)
    u,v,w,p,rho,mf= extract_data(0,posName,DataName_fmt,50000,52500,50)
    z_max = z_lim(0)
    z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
    print(mf.size)
    print(mf.shape)
    print(y.size * z.size)
    print(z.size)
    print(y.size)
    mf_avg = np.squeeze(np.mean(mf,axis=0))
    mf_tot = np.trapezoid(np.trapezoid(mf,z),y)
    
    print(np.mean(mf_tot))

