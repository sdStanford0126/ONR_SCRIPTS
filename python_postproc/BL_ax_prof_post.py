#THis is being used to post processing time averaged
#results for the BL_interior_probes (axial profiles)

import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpi4py import MPI
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
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

#thermodymaics reference for Power-Law viscosity
mu_ref = 3.26788998e-6
T_ref = 1.0
n = 0.7
gamma = 1.4

def muCalc(T):
    return mu_ref*np.power((T/T_ref),n)

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


def _load_timestep(DataName_fmt, tid, ind, X_ind, Ny, Nz):
    DataName = DataName_fmt.format(int(tid))
    DataR = np.loadtxt(DataName, skiprows=1)
    DataRo = np.zeros(np.shape(DataR))
    DataRo[ind, :] = DataR[:, :]
    return (
        np.reshape(DataRo[X_ind[0], 0], (Ny, Nz)),
        np.reshape(DataRo[X_ind[0], 1], (Ny, Nz)),
        np.reshape(DataRo[X_ind[0], 2], (Ny, Nz)),
        np.reshape(DataRo[X_ind[0], 3], (Ny, Nz)),
        np.reshape(DataRo[X_ind[0], 4], (Ny, Nz)),
    )

#first need to get index from .pxyz and record x,y,z locatiions (no read original probe file)
#then for a given x value, find the closest plane and extract the data from
# that plane and organize them in the y,z plane
#TODO: modify by using a dictionary of dictionaries to get to individual axial profiles at one go (instead of reading each file again at each x location)
def extract_data_serial(x,posName,DataName_fmt,tid_str,tid_end,dt):
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
    for i,tid in tqdm(enumerate(tids), total=Nt, desc=f"x={x:.3f}"):
        DataName = DataName_fmt.format(int(tid))
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

def extract_data(x,posName,DataName_fmt,tid_str,tid_end,dt,max_workers=None):
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
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_load_timestep, DataName_fmt, tid, ind, X_ind, Ny, Nz): i
                   for i, tid in enumerate(tids)}
        for future in tqdm(as_completed(futures), total=Nt, desc=f"x={x:.3f}"):
            i = futures[future]
            u[i], v[i], w[i], p[i], rho[i] = future.result()
        
    
    print(u.size)    
    print(Nt*Ny*Nz)
    mf = u*rho
    return u,v,w,p,rho,mf
   
#TODO: plot various profile data at the 4 different lines defined with in BL_probes and assmemble these quantities
#Avarilable quantities, by order are:
#u,v,w,p,rho
#we would like to at least assemble the following:
#avg(u) profile, avg(u'v'), avg(u'w'), 
def plotTurbProf_C(u,v,w,rho,p,xc,y,z, out_dir = "./",eps=delta/2.0, fig_offset=0):
    """
    takes in time history and the 1d y,z vector and produces the following profiles at given streamwise station x = xc
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
    u_avg  = np.mean(u,axis=0)
    v_avg  = np.mean(v,axis=0)
    w_avg  = np.mean(w,axis=0)
    rho_avg= np.mean(rho,axis=0)  
    p_avg  = np.mean(p,axis=0)
    #get flucutations
    uf   = u - u_avg
    vf   = v - v_avg
    wf   = w - w_avg

    u_avg_y = np.squeeze(u_avg[:,z_ind])
    #print("u_avg_y_size: ",u_avg_y.size)
    #print("y: ",y.size)
    u_avg_z = np.squeeze(u_avg[y_ind,:])     
    #print("u_avg_z_size: ",u_avg_z.size)
    #print("y: ",z.size)

    rho_avg_y = np.squeeze(rho_avg[:,z_ind])
    rho_avg_z = np.squeeze(rho_avg[y_ind,:])

    p_avg_y = np.squeeze(p_avg[:,z_ind])
    p_avg_z = np.squeeze(p_avg[y_ind,:])

    T_avg_y = p_avg_y/rho_avg_y
    T_avg_z = p_avg_z/rho_avg_z

    mu_avg_y = muCalc(T_avg_y)
    mu_avg_z = muCalc(T_avg_z)

    mu_w_y = (mu_avg_y[0]+mu_avg_y[-1])/2.0
    mu_w_z = (mu_avg_z[0]+mu_avg_z[-1])/2.0
    dy = np.abs(y[1]-y[0])
    dz = np.abs(z[1]-z[0])
    dudy_w1 = (-3*u_avg_y[0] + 4*u_avg_y[1] - u_avg_y[2])/(2*dy)
    dudz_w1 = (-3*u_avg_z[0] + 4*u_avg_z[1] - u_avg_z[2])/(2*dz)
    #dudy_w2 = (3*u_avg_y[-1] - 4*u_avg_y[-2] + u_avg_y[-3])/(2*dy)
    #dudz_w2 = (3*u_avg_z[-1] - 4*u_avg_z[-2] + u_avg_z[-3])/(2*dz)
    #dudy_w  = (dudy_w1+dudy_w2)/2
    #dudz_w  = (dudz_w1+dudz_w2)/2
    dudy_w = dudy_w1
    dudz_w = dudz_w1

    print("dudy at wall for x_c=%.2f is %.4f" % (xc,dudy_w))
    print("dudz at wall for x_c=%.2f is %.4f" % (xc,dudz_w))


    tauy_w = mu_w_y*dudy_w
    tauz_w = mu_w_z*dudz_w

    print("tau_y at wall for x_c=%.2f is %.4f" % (xc,tauy_w))
    print("tau_z at wall for x_c=%.2f is %.4f" % (xc,tauz_w))
    rhoy_w = (rho_avg_y[0]+rho_avg_y[-1])/2.0
    rhoz_w = (rho_avg_z[0]+rho_avg_z[-1])/2.0
    print("rho_y at wall for x_c=%.2f is %.4f" % (xc,rhoy_w))
    print("rho_z at wall for x_c=%.2f is %.4f" % (xc,rhoz_w))

    u_tau_y = np.sqrt(tauy_w/rhoy_w) 
    u_tau_z = np.sqrt(tauz_w/rhoz_w) 
    print("friction vel for y is ", u_tau_y)
    print("friction vel for z is ", u_tau_z)
    #build fluctuation
    ufvf_avg = np.mean(uf*vf,axis=0)
    ufvf_avg_y = np.squeeze(ufvf_avg[:,z_ind])

    ufwf_avg = np.mean(uf*wf, axis=0)
    ufwf_avg_z = np.squeeze(ufwf_avg[y_ind,:])


    PlotName_fmt = "{}_{}.png"
    PlotName_fmt = os.path.join(out_dir, PlotName_fmt)

    fig_count = 7
    #u_avg_y
    fig_id = 0 + fig_count*fig_offset
    var_label = "u_avg"
    side_label = "y"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,y,u_avg_y,r"$\overline{u}$", side_label,PlotName,fig_id)

    
    #u_avg_z
    fig_id = 1 + fig_count*fig_offset
    var_label = "u_avg"
    side_label = "z"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,z,u_avg_z,r"$\overline{u}$", side_label,PlotName,fig_id)

    #ufvf_avg_y
    fig_id = 2+ fig_count*fig_offset
    var_label = "ufvf_avg"
    side_label = "y"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,y,ufvf_avg_y,r"$\overline{u'v'}$", side_label,PlotName,fig_id)

    #ufwf_avg_z
    fig_id = 3 + fig_count*fig_offset
    var_label = "ufwf_avg"
    side_label = "z"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,z,ufwf_avg_z,r"$\overline{u'w'}$", side_label,PlotName,fig_id)

    #ufvf_avg_y
    fig_id = 4+ fig_count*fig_offset
    var_label = "ufvf_avg_div_u_tau_y"
    side_label = "y"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,y,ufvf_avg_y/u_tau_y,r"$\frac{\overline{u'v'}}{u_{\tau,y}}$", side_label,PlotName,fig_id)

    #ufwf_avg_z
    fig_id = 5 + fig_count*fig_offset
    var_label = "ufwf_avg_div_u_tau_z"
    side_label = "z"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotProfile(xc,z,ufwf_avg_z/u_tau_z,r"$\frac{\overline{u'w'}}{u_{\tau,z}}$", side_label,PlotName,fig_id)

def plotProfile(xc,pos,var,var_label:str,side_label:str,PlotName: str,fig_id):
    plt.figure(fig_id)
    plt.plot(var,pos, label = "x = %.2f" % xc)
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


def plotAxialProf(xc,var,z,y,var_label:str,out_dir):
    PlotName_fmt = "{}_x_{}.png"
    PlotName_fmt = os.path.join(out_dir,PlotName_fmt)
    PlotName = PlotName_fmt.format(var_label,xc) 
    plt.figure()
    plt.pcolormesh(y,z,np.transpose(var),shading='nearest')
    plt.xlabel("y")
    plt.ylabel("z")
    plt.colorbar()
    plt.title(var_label)
    plt.tight_layout()
    plt.axis('equal')
    plt.savefig(PlotName,dpi=300,bbox_inches="tight")
    
def main():
        #these are for test purposes
    TBL_test_cases = ["0.0125", "0.01875", "0.025"]
    tid_str_cases  = [40000, 40000, 130000]
    tid_end_cases  = [100000, 100000, 190000]

    #TBL_test_cases = ["0.025"]
    #tid_str_cases  = [130000]
    #tid_end_cases  = [160000]
    data_dir_fmt = "/anvil/scratch/x-sdai/BL_test_baseline_{:s}/pcprobe_int_axprof"
    out_dir_fmt = "/anvil/scratch/x-sdai/BL_post_proc/BL_{:s}"
    for i, case in enumerate(TBL_test_cases):
        print("Now processing for TBL case: ", case)
        data_dir = data_dir_fmt.format(case)
        print(xs)
        z_max = z_lim(0)
        z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
        out_dir = out_dir_fmt.format(case)
        posName = os.path.join(data_dir,"int_axprof.pxyz")
        fname_fmt = "int_axprof.{:08d}.pcd"
        DataName_fmt = os.path.join(data_dir,fname_fmt)
        tid_str = tid_str_cases[i]
        tid_end = tid_end_cases[i]
        for x in xs:
            u,v,w,p,rho,mf= extract_data(x,posName,DataName_fmt,tid_str,tid_end,50)
            z_max = z_lim(x)
            #print(y)
            T = p/rho
            c = np.sqrt(gamma*T)
            z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
            print("mass flow at x=%.2f is mf = %.2f " % (x,evalMfAvg(mf,z,y)))
            plotTurbProf_C(u,v,w,rho,p,x,y,z,out_dir = out_dir)
            plotAxialProf(x,np.mean(u,axis=0),z,y,"u_avg",out_dir)
            plotAxialProf(x,np.mean(u,axis=0)/np.mean(c,axis=0),z,y,"Mach_u",out_dir)
            

    #load baseline h5 make the same plots based using the same functions
    #need to supply new x,y,z
    #grab x slices
if __name__ == "__main__":
    t_start = time.perf_counter()
    main()
    t_end = time.perf_counter()
    print(f"Total execution time: {t_end - t_start:.2f} s")



