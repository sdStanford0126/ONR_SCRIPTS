#get 3D data of the interior shock structure 
import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py
#from mpi4py import MPI
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed



def _load_timestep_grad(DataName_fmt, tid, ind):
    # columns: grad_rho_x(0) grad_rho_y(1) grad_rho_z(2) u(3) v(4) w(5) p(6) rho(7)
    DataName = DataName_fmt.format(int(tid))
    DataR = np.loadtxt(DataName, skiprows=1)
    DataRo = np.zeros(np.shape(DataR))
    DataRo[ind, :] = DataR[:, :]
    return (
        DataRo[:, 0],
        DataRo[:, 1],
        DataRo[:, 2],
        DataRo[:, 3],
        DataRo[:, 4],
        DataRo[:, 5],
        DataRo[:, 6],
        DataRo[:, 7],
    )

def buildTimeRecord(posName,dataName_fmt,tids,max_workers=None):
    """
    INPUTs:
    posName: *.pxyz (the file that contains the positional data)
    dataName_fmt: *.{%08d}.pcd  datafiles for the entire grid
    tids: timestep indices: tid_start:dt:tid_end
    
    Due to the strcture at play here, will not attempt to reshape into 3D spatial representation
    as Nz = f(x) (number of pts on Minor dir. is a func. of x)
    
    OUTPUTs:
    X,Y,Z: flat arrays of positive in X (streamwise), Y(spanwise, major width), Z(transverse, minor height)
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
    Nt = np.size(tids)
    #storage variables
    grad_rho_x = np.zeros((Nt,Npts))
    grad_rho_y = np.zeros((Nt,Npts))
    grad_rho_z = np.zeros((Nt,Npts))
    u   = np.zeros((Nt,Npts))
    v   = np.zeros((Nt,Npts))
    w   = np.zeros((Nt,Npts))
    p   = np.zeros((Nt,Npts))
    rho = np.zeros((Nt,Npts))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        print("active threads: ", executor._max_workers)
        futures = {executor.submit(_load_timestep_grad, dataName_fmt, tid, ind): i
                   for i, tid in enumerate(tids)}
        for future in tqdm(as_completed(futures), total=Nt):
            i = futures[future]
            grad_rho_x[i], grad_rho_y[i], grad_rho_z[i], u[i], v[i], w[i], p[i], rho[i] = future.result()
    #datasize validation
    print("valiate that extracted variable size aligns with expectations")
    print(u.size)
    print(Nt*Npts)
    mf = u*rho
    return grad_rho_x, grad_rho_y, grad_rho_z, u, v, w, p, rho,mf,X,Y,Z

def _calc_c(p,rho,R=1.0/1.4,gamma=1.4):
    """
    computes speed of sound (c) assuming ideal gas laws
    """
    return np.sqrt(gamma*R*p/rho)
def calc_avg(var_time,Nt):
    return var_time/Nt

def buildAvgData(posName,dataName_fmt,tids):
    """
    serial code doing the following:
    1. read a given file at a given tid
    2. extract the variables and add to sum
    3. increment the avg counter
    4. compute the the average
    5. output the average
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
    #build the sum variables
    #raw data variables
    grad_rho_x = np.zeros(Npts)
    grad_rho_y = np.zeros(Npts)
    grad_rho_z = np.zeros(Npts)
    u   = np.zeros(Npts)
    v   = np.zeros(Npts)
    w   = np.zeros(Npts)
    p   = np.zeros(Npts)
    rho = np.zeros(Npts)
    #computed variables
    mf = np.zeros(Npts)
    Mach = np.zeros(Npts)
    Nsamp = np.size(tids)
    for i,tid in tqdm(enumerate(tids),total=np.size(tids)):
        dataName = dataName_fmt.format(tid)
        DataR = np.loadtxt(dataName,skiprows=1)
        DataRo = np.zeros(np.shape(DataR))
        DataRo[ind,:] = DataR[:,:]
        #grab data at timestep
        grad_rho_x_temp = DataRo[:,0]
        grad_rho_y_temp = DataRo[:,1]
        grad_rho_z_temp = DataRo[:,2]
        u_temp          = DataRo[:,3]
        v_temp          = DataRo[:,4]
        w_temp          = DataRo[:,5]
        p_temp          = DataRo[:,6]
        rho_temp        = DataRo[:,7]

        mf_temp         = u_temp*rho_temp
        c_temp          = _calc_c(p_temp,rho_temp)
        vel_mag_temp    = np.sqrt(u_temp**2 + v_temp**2 + w_temp**2)
        Mach_temp       = vel_mag_temp/c_temp
        
        
        grad_rho_x = grad_rho_x+grad_rho_x_temp
        grad_rho_y = grad_rho_y+grad_rho_y_temp
        grad_rho_z = grad_rho_z+grad_rho_z_temp

        u = u+u_temp
        v = v+v_temp
        w = w+w_temp
        p = p+p_temp
        rho = rho+rho_temp
        mf = mf+mf_temp
        Mach = Mach+Mach_temp

        grad_rho_x_avg     =calc_avg(grad_rho_x,Nsamp)
        grad_rho_y_avg     =calc_avg(grad_rho_y,Nsamp)
        grad_rho_z_avg     =calc_avg(grad_rho_z,Nsamp)
        u_avg     =calc_avg(u,Nsamp)
        v_avg     = calc_avg(v,Nsamp)
        w_avg     = calc_avg(w,Nsamp)
        p_avg     = calc_avg(p,Nsamp)
        rho_avg   = calc_avg(rho,Nsamp)
    return grad_rho_x_avg,grad_rho_y_avg,grad_rho_z_avg,u_avg,v_avg,w_avg,p_avg,rho_avg,X,Y,Z

def _load_timestep_avg(dataName_fmt, tid, ind):
    DataName = dataName_fmt.format(int(tid))
    DataR = np.loadtxt(DataName, skiprows=1)
    DataRo = np.zeros(np.shape(DataR))
    DataRo[ind, :] = DataR[:, :]
    u_t   = DataRo[:, 3]
    v_t   = DataRo[:, 4]
    w_t   = DataRo[:, 5]
    p_t   = DataRo[:, 6]
    rho_t = DataRo[:, 7]
    mf_t  = u_t * rho_t
    Mach_t = np.sqrt(u_t**2 + v_t**2 + w_t**2) / _calc_c(p_t, rho_t)
    return DataRo[:, 0], DataRo[:, 1], DataRo[:, 2], u_t, v_t, w_t, p_t, rho_t, mf_t, Mach_t

def buildAvgData_parallel(posName, dataName_fmt, tids, max_workers=None):
    """
    Parallelized version of buildAvgData using ThreadPoolExecutor.
    Each timestep is loaded concurrently; results are accumulated in the main thread.

    Returns same variables as buildAvgData plus mf_avg, Mach_avg, X, Y, Z.
    """
    Pos = np.loadtxt(posName, skiprows=1)
    ind = Pos[:, 0].astype(int)
    Npts = np.size(ind)
    X = np.zeros(Npts); X[ind] = Pos[:, 1]
    Y = np.zeros(Npts); Y[ind] = Pos[:, 2]
    Z = np.zeros(Npts); Z[ind] = Pos[:, 3]

    Nsamp = np.size(tids)
    grad_rho_x = np.zeros(Npts)
    grad_rho_y = np.zeros(Npts)
    grad_rho_z = np.zeros(Npts)
    u   = np.zeros(Npts)
    v   = np.zeros(Npts)
    w   = np.zeros(Npts)
    p   = np.zeros(Npts)
    rho = np.zeros(Npts)
    mf  = np.zeros(Npts)
    Mach = np.zeros(Npts)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        print("active threads: ", executor._max_workers)
        futures = {executor.submit(_load_timestep_avg, dataName_fmt, tid, ind): i
                   for i, tid in enumerate(tids)}
        for future in tqdm(as_completed(futures), total=Nsamp):
            grx, gry, grz, u_t, v_t, w_t, p_t, rho_t, mf_t, Mach_t = future.result()
            grad_rho_x += grx;  grad_rho_y += gry;  grad_rho_z += grz
            u += u_t;  v += v_t;  w += w_t
            p += p_t;  rho += rho_t;  mf += mf_t;  Mach += Mach_t

    return (grad_rho_x/Nsamp, grad_rho_y/Nsamp, grad_rho_z/Nsamp,
            u/Nsamp, v/Nsamp, w/Nsamp, p/Nsamp, rho/Nsamp, mf/Nsamp, Mach/Nsamp,
            X, Y, Z)


def write_time_history_h5(grad_rho_x, grad_rho_y, grad_rho_z,
                                u, v, w, p, rho, mf,
                                X,Y,Z, tids,
                                caseName: str, out_dir="./", max_workers=None):
    """
    Write the full time history (grad_rho_* + flow vars) for a single axial
    plane to HDF5.  Direct counterpart to write_time_history_h5 for data
    produced by extract_data_grad.

    Parameters
    ----------
    grad_rho_x, grad_rho_y, grad_rho_z : ndarray, shape (Nt, Npts)
    u, v, w, p, rho, mf                : ndarray, shape (Nt, Npts)
    X,Y,Z: position data of the gird
    tid_str, tid_end : int — timestep range (same as passed to extract_data_grad)
    dt         : int   — timestep stride
    caseName   : str   — prepended to the filename
    out_dir    : str   — output directory
    max_workers: int | None — thread count for parallel dataset writes

    Returns
    -------
    str  — absolute path of the written HDF5 file
    """
    fname = os.path.join(out_dir, f"{caseName}_noz_int_grad_time_history_full.h5")

    flow_vars = {
        "grad_rho_x": grad_rho_x,
        "grad_rho_y": grad_rho_y,
        "grad_rho_z": grad_rho_z,
        "u":   u,
        "v":   v,
        "w":   w,
        "p":   p,
        "rho": rho,
        "mf":  mf,
    }

    with h5py.File(fname, "w") as hf:
        hf.create_dataset("tids", data=tids)
        hf.create_dataset("X",    data=X)
        hf.create_dataset("Y",    data=Y)
        hf.create_dataset("Z",    data=Z)

        for name, arr in flow_vars.items():
            hf.create_dataset(name, shape=arr.shape, dtype=arr.dtype,
                              compression="gzip", compression_opts=4)

        def _write_var(name):
            hf[name][...] = flow_vars[name]
            return name

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_write_var, name): name for name in flow_vars}
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc=f"writing h5 grad"):
                fut.result()

    print(f"  wrote: {fname}")
    return fname

def write_avg_h5(grad_rho_x_avg, grad_rho_y_avg, grad_rho_z_avg,
                      u_avg, v_avg, w_avg, p_avg, rho_avg, mf_avg, Mach_avg,
                      X, Y, Z,T,
                      caseName: str, out_dir="./"):
    """
    Write the time-averaged flow fields produced by buildAvgData_parallel to HDF5.

    Parameters
    ----------
    grad_rho_x_avg, ..., Mach_avg : ndarray, shape (Npts,)  — averaged variables
    X, Y, Z                       : ndarray, shape (Npts,)  — probe coordinates
    caseName : str   — prepended to the filename
    out_dir  : str   — output directory
    T        :       - time (in ATU) that the result has been averaged over
    Returns
    -------
    str  — absolute path of the written HDF5 file
    """
    fname = os.path.join(out_dir, f"{caseName}_noz_int_avg_T_{T}.h5")

    avg_vars = {
        "grad_rho_x_avg": grad_rho_x_avg,
        "grad_rho_y_avg": grad_rho_y_avg,
        "grad_rho_z_avg": grad_rho_z_avg,
        "u_avg":    u_avg,
        "v_avg":    v_avg,
        "w_avg":    w_avg,
        "p_avg":    p_avg,
        "rho_avg":  rho_avg,
        "mf_avg":   mf_avg,
        "Mach_avg": Mach_avg,
    }

    with h5py.File(fname, "w") as hf:
        hf.create_dataset("X", data=X)
        hf.create_dataset("Y", data=Y)
        hf.create_dataset("Z", data=Z)
        for name, arr in avg_vars.items():
            hf.create_dataset(name, data=arr, compression="gzip", compression_opts=4)

    print(f"  wrote: {fname}")
    return fname

def read_avg_h5(fname: str):
    """
    Read the time-averaged HDF5 file written by write_avg_h5.

    Parameters
    ----------
    fname : str — path to the HDF5 file

    Returns
    -------
    grad_rho_x_avg, grad_rho_y_avg, grad_rho_z_avg,
    u_avg, v_avg, w_avg, p_avg, rho_avg, mf_avg, Mach_avg,
    X, Y, Z  — all ndarray, shape (Npts,)
    """
    with h5py.File(fname, "r") as hf:
        X             = hf["X"][:]
        Y             = hf["Y"][:]
        Z             = hf["Z"][:]
        grad_rho_x_avg = hf["grad_rho_x_avg"][:]
        grad_rho_y_avg = hf["grad_rho_y_avg"][:]
        grad_rho_z_avg = hf["grad_rho_z_avg"][:]
        u_avg          = hf["u_avg"][:]
        v_avg          = hf["v_avg"][:]
        w_avg          = hf["w_avg"][:]
        p_avg          = hf["p_avg"][:]
        rho_avg        = hf["rho_avg"][:]
        mf_avg         = hf["mf_avg"][:]
        Mach_avg       = hf["Mach_avg"][:]
    return (grad_rho_x_avg, grad_rho_y_avg, grad_rho_z_avg,
            u_avg, v_avg, w_avg, p_avg, rho_avg, mf_avg, Mach_avg,
            X, Y, Z)

def plot_slices(x,y,z,y_pos,delta,data,cmap='viridis',ub = np.inf, lb = -np.inf, title='test',outdir="./",figname="test.png"):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    data_min = np.min(data)
    data_max = np.max(data)
    data_range = data_max - data_min
    data_mean = np.mean(data)
    data_center = (data_max-data_min)/2.0 + data_min
    max_section = np.maximum(data_mean-data_min,data_max-data_mean)
    alpha = np.abs(data)/(max_section)
    from matplotlib import cm,colors
    for i, y_c in enumerate(y_pos):
        condition = (y > (y_c - delta/2)) & (y < (y_c + delta/2))
        x_plot=x[condition]
        y_plot=y[condition]
        z_plot=z[condition]
        data_plot=data[condition]
        alpha_plot = alpha[condition]
        norm_c = colors.Normalize(vmin=-3, vmax=3)
        cmap = cm.viridis

        norm_a = colors.Normalize(vmin=alpha.min(), vmax=alpha.max())
        alpha_plot = norm_a(alpha_plot)
        
        #print(alpha.max)
        #print(alpha.min)
        print(np.min(alpha))
        print(np.max(alpha))
        rgba = cmap(norm_c(data_plot))
        rgba[:,3]=alpha_plot
        sp=ax.scatter(x_plot,y_plot,z_plot,c=rgba)
    # 4. Set axes labels and title
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.title(title)
    #plt.colorbar(sp)
    # 5. save plot
    fig.savefig(os.path.join(outdir,figname), dpi=300,bbox_inches="tight")
    plt.close(fig)
 

def plot_scatter_contour_3d(x, y, z, data, cmap='viridis',tgt_val= 1.0, delta=0.05,marker_size=50, alpha=0.8, title='test',outdir="./",figname="test.png"):
    """
    Plots a 3D scatter plot where the 4th dimension (data) is represented by a colormap.

    Parameters:
    x (array-like): 1D array of X coordinates.
    y (array-like): 1D array of Y coordinates.
    z (array-like): 1D array of Z coordinates.
    data (array-like): 1D array of data values corresponding to the (x, y, z) spatial coordinates.
    cmap (str): Colormap to use for the data values. Default is 'viridis'.
    marker_size (float): Size of the scatter plot markers. Default is 50.
    alpha (float): Transparency of the markers (0 to 1). Default is 0.8.
    title (str): Title of the plot.
    """
    # 1. Initialize the figure and 3D axis
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    condition = (data > (tgt_val - delta)) & (data < (tgt_val+delta)) #bounds
    x_plot=x[condition]
    y_plot=y[condition]
    z_plot=z[condition]
    data_plot=data[condition]
    # 2. Create the scatter plot mapping 'data' to the color parameter 'c'
    scatter = ax.scatter(x_plot, y_plot, z_plot, c=data_plot, cmap=cmap, s=marker_size, alpha=alpha)

    # 3. Add a colorbar to act as the contour legend
    cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
    cbar.set_label('Data Values')

    # 4. Set axes labels and title
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.set_zlabel('Z Axis')
    plt.title(title)

    # 5. save plot
    fig.savefig(os.path.join(outdir,figname), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    L = 2.56 #baselength
    LEVEL=9 #mesh refinement level
    delta = L/2**LEVEL

    #diverging nozzle grid assessment 
    #Nz is in fact a function of x
    x_str = -1.0152632
    x_end = 0.0
    z_str = 0.42498737
    z_end = 0.5
    L = 2.56 #baselength
    LEVEL=9 #mesh refinement level
    delta = L/2**LEVEL


    ##################
    #I/O specfications here
    #################

    data_dir_fmt = "/anvil/scratch/x-sdai/AR2_{:s}/pcprobes_noz_int"
    out_dir = "/anvil/scratch/x-sdai/noz_int_post"

    caseNames = ["base_151M_akhil_restart", "port_252M_V8", "port_252M_V10"]
    titles = [rf"Baseline 151M", rf"Port Injection 252M, $IPR=1.45$", rf"Port Injection 252M, $IPR=1.98"]
    tid_strs = [884550,1225500,750]
    tid_ends = [984600,1683000,500250]
    dts      = [150,750,750]
    Ts = [(tid_end - tid_str)/dt * 0.15 for tid_str,tid_end,dt in zip(tid_strs,tid_ends,dts)]
    print("total sample times are: ", Ts)
   #delta_t = [1e-3,2e-4,2e-4]
    #files naming scheme
    posName = "nov_int.pxyz"
    dataName_fmt = "nov_int.{:08d}.pcd"

    debug_flag = False
    if debug_flag:
        cpus = 20
    else:
        cpus = int(os.getenv("SLURM_NTASKS"))
    print("assigned core count: ", cpus)
    figname_fmt = "Mach_3D_scatter_{:s}.png"

    for i,caseName in enumerate(caseNames):
        data_dir = data_dir_fmt.format(caseName)
        tid_str = tid_strs[i]
        tid_end = tid_ends[i]
        dt      = dts[i]
        T       = int(Ts[i])
        tids = np.arange(tid_str,tid_end+dt,dt)
        print(np.size(tids))
        #grad_rho_x, grad_rho_y, grad_rho_z, u, v, w, p, rho,mf,X,Y,Z = buildTimeRecord(
        #    os.path.join(data_dir,posName),os.path.join(data_dir,dataName_fmt),
        #    tids,max_workers=cpus)
        #write_time_history_h5_grad(grad_rho_x, grad_rho_y, grad_rho_z,
        #                           u, v, w, p, rho, mf,
        #                           X,Y,Z, tids,caseName, 
        #                           out_dir, max_workers=cpus) 
        #u_avg   = np.mean(u,axis=0)
        #v_avg   = np.mean(v,axis=0)
        #w_avg   = np.mean(w,axis=0)
        #p_avg   = np.mean(p,axis=0)
        #rho_avg = np.mean(rho,axis=0)

        #R = 1
        #gamma = 1.4 
        #T_avg = p_avg/(rho_avg*R)

        #vel_mag_avg = np.sqrt(u_avg**2 + v_avg**2 + w_avg**2)
        #c_avg = np.sqrt(gamma*R*T_avg)
        #Mach_avg = vel_mag_avg/c_avg
        (grad_rho_x_avg,grad_rho_y_avg,grad_rho_z_avg,
         u_avg,v_avg,w_avg,p_avg,rho_avg,
         mf_avg,Mach_avg
         ,X,Y,Z)=buildAvgData_parallel(os.path.join(data_dir,posName), os.path.join(data_dir,dataName_fmt), tids, max_workers=cpus)
        plot_scatter_contour_3d(X,Y,Z,Mach_avg,tgt_val=1.0, delta=0.05,title=titles[i],outdir=out_dir,figname=figname_fmt.format(caseName))
        write_avg_h5(grad_rho_x_avg, grad_rho_y_avg, grad_rho_z_avg,
                      u_avg, v_avg, w_avg, p_avg, rho_avg, mf_avg, Mach_avg,
                      X, Y, Z,T,
                      caseName, out_dir)
        #fname = os.path.join(out_dir, f"{caseName}_noz_int_avg_T_{T}.h5")
        #grad_rho_x_avg, grad_rho_y_avg, grad_rho_z_avg,u_avg, v_avg, w_avg, p_avg, rho_avg, mf_avg, Mach_avg,X, Y, Z = read_avg_h5(fname)

        plot_scatter_contour_3d(X, Y, Z, grad_rho_x_avg, cmap='viridis',tgt_val= 10, delta=0.05,marker_size=50, alpha=0.8, title='test',outdir="./",figname="test.png") 
        plot_slices(X,Y,Z,np.linspace(0,1-delta,10),delta,grad_rho_x_avg)
if __name__ =="__main__":
    main()