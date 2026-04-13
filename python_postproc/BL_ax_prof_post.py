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
    #print("data size", Npts)
    #print("expected data size", np.sum(np.fromiter((Ny*Nz for Nz in Nzs),int)))
    #plt.figure()
    #plt.plot(X)
    #plt.plot(Y)
    #plt.plot(Z)
    #plt.savefig("test.png")
    # grab time series data from the given x position
    X_ind = np.where(np.abs(X-x) < 1e-6)
    Nz = findNz(x,delta)
    print(X_ind[0].size)
    print(Nz*Ny)
    #now read data
    tids = np.arange(tid_str,tid_end,dt)
    Nt = np.size(tids)
    u  = np.zeros((Nt,Ny,Nz))
    v  = np.zeros((Nt,Ny,Nz))
    w  = np.zeros((Nt,Ny,Nz))
    p  = np.zeros((Nt,Ny,Nz))
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
   
#TODO: plot various profile data at the 4 different lines defined with in BL_probes and assmemble these quantities
#Avarilable quantities, by order are:
#u,v,w,p,rho
#we would like to at least assemble the following:
#avg(u) profile, avg(u'v'), avg(u'w'), 
def plotTurbProf_C(u,v,w,xc,y,z, out_dir = "./",eps=delta/2.0):
    """
    takes in time history and the 1d y,z vector and produces the following profiles
    u_avg over y centerline (major axis)
    u_avg over z centerline (minor axis)
    u'v' over y centerline (major axis)
    u'w' over z centerline (minor axis)
    tke over y centerline (major axis)
    tke over z centerline (minor axis)
    """
    #find centerline 
    y_ind = np.where(np.abs(y) < eps)
    z_ind = np.where(np.abs(z) < eps)
    y_ind = y_ind[0]
    z_ind = z_ind[0]
    #debug (print loc)
    print(y[y_ind])
    print(z[z_ind])
    
    #get time averages
    u_avg = np.mean(u,axis=0)
    v_avg = np.mean(v,axis=0)
    w_avg = np.mean(w,axis=0)
    #print(u_avg.shape)

    #get flucutations
    uf   = u - u_avg
    vf   = v - v_avg
    wf   = w - w_avg

    u_avg_y = u_avg[:,z_ind]
    #print("u_avg_y_size: ",u_avg_y.size)
    #print("y: ",y.size)
    u_avg_z = u_avg[y_ind,:]     
    #print("u_avg_z_size: ",u_avg_z.size)
    #print("y: ",z.size)

    #build fluctuation
    ufvf_avg = np.mean(uf*vf,axis=0)
    ufvf_avg_y = ufvf_avg[:,z_ind]

    ufwf_avg = np.mean(uf*wf, axis=0)
    ufwf_avg_z = ufwf_avg[y_ind,:]


    PlotName_fmt = "{}_{}.png"
    PlotName_fmt = os.path.join(out_dir, PlotName_fmt)
    #u_avg_y
    fig_id = 1
    var_label = "u_avg"
    side_label = "y"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,y,u_avg_y,r"$\overline{u}$", side_label,PlotName,fig_id)

    
    #u_avg_z
    fig_id = 2
    var_label = "u_avg"
    side_label = "z"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,z,u_avg_z,r"$\overline{u}$", side_label,PlotName,fig_id)

    #ufvf_avg_y
    fig_id = 3
    var_label = "ufvf_avg"
    side_label = "y"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,y,ufvf_avg_y,r"$\overline{u}$", side_label,PlotName,fig_id)

    #ufvf_avg_z
    fig_id = 3
    var_label = "ufvf_avg"
    side_label = "y"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,y,ufvf_avg_y,r"$\overline{u}$", side_label,PlotName,fig_id)

def plotProfile(xc,pos,var,var_label:str,side_label:str,PlotName: str,fig_id):
    plt.figure(fig_id)
    plt.plot(pos,var, label = "x = %.2f" % xc)
    plt.xlabel(var_label)
    plt.ylabel(side_label)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PlotName, dpi = 300, bbox_inches="tight")

def evalMfAvg(mf,z,y):
    mf_avg = np.squeeze(np.mean(mf,axis=0))
    mf_tot = np.trapezoid(np.trapezoid(mf,z),y)
    mf_tot = np.mean(mf_tot) 
    return mf_tot

if __name__ == "__main__":
    #these are for test purposes
    data_dir = "/anvil/scratch/x-sdai/BL_test_baseline_0.025/pcprobe_int_axprof"
    posName = os.path.join(data_dir,"int_axprof.pxyz")
    fname_fmt = "int_axprof.{:08d}.pcd"
    DataName_fmt = os.path.join(data_dir,fname_fmt)
    print(xs)
    u,v,w,p,rho,mf= extract_data(0,posName,DataName_fmt,162200,163000,50)
    z_max = z_lim(0)
    #print(y)
    z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
    print(mf.size)
    print(mf.shape)
    print(y.size * z.size)
    print(z.size)
    print(y.size)
    print(evalMfAvg(mf,z,y))
    out_dir = "/anvil/scratch/x-sdai/"
    plotTurbProf_C(u,v,w,0,y,z) 

    for x in xs:
        u,v,w,p,rho,mf= extract_data(x,posName,DataName_fmt,162200,163000,50)
        z_max = z_lim(x)
        #print(y)
        z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
        print("mass flow at x=%.2f is mf = %.2f " % (x,evalMfAvg(mf,z,y)))
        plotTurbProf_C(u,v,w,x,y,z)


