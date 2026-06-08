#rewritten to take advantage of existing pcprobes_noz_int data for boundary layer processing
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
L = 2.56 #baselength
LEVEL=9 #mesh refinement level
delta = L/2**LEVEL

#Nz is in fact a function of x
x_str = -1.0152632
x_end = 0.0
z_str = 0.42498737
z_end = 0.5
slope = (z_end-z_str)/(x_end-x_str)

#old version (used for BL_test_<roughness_height>)
xs = np.linspace(-1,0,5)
#new version (used for AR2_base_<*>)
#Nx = 101 
#xs = np.linspace(-1,0,Nx) #direct defined from write_BL_interior_probes.py

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

def _load_timestep_grad(DataName_fmt, tid, ind, X_ind, Ny, Nz):
    # columns: grad_rho_x(0) grad_rho_y(1) grad_rho_z(2) u(3) v(4) w(5) p(6) rho(7)
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
        np.reshape(DataRo[X_ind[0], 5], (Ny, Nz)),
        np.reshape(DataRo[X_ind[0], 6], (Ny, Nz)),
        np.reshape(DataRo[X_ind[0], 7], (Ny, Nz)),
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
    tids = tids.astype(int) #type cast to int
    Nt = np.size(tids)
    u  = np.zeros((Nt,Ny,Nz))
    v  = np.zeros((Nt,Ny,Nz))
    w  = np.zeros((Nt,Ny,Nz))
    p  = np.zeros((Nt,Ny,Nz))
    rho  = np.zeros((Nt,Ny,Nz))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        print("active threads: ", executor._max_workers)
        futures = {executor.submit(_load_timestep, DataName_fmt, tid, ind, X_ind, Ny, Nz): i
                   for i, tid in enumerate(tids)}
        for future in tqdm(as_completed(futures), total=Nt, desc=f"x={x:.3f}"):
            i = futures[future]
            u[i], v[i], w[i], p[i], rho[i] = future.result()
        
    
    print(u.size)
    print(Nt*Ny*Nz)
    mf = u*rho
    return u,v,w,p,rho,mf

def extract_data_grad(x, posName, DataName_fmt, tid_str, tid_end, dt, max_workers=None):
    # columns in .pcd: grad_rho_x(0) grad_rho_y(1) grad_rho_z(2) u(3) v(4) w(5) p(6) rho(7)
    Pos  = np.loadtxt(posName, skiprows=1)
    ind  = Pos[:, 0].astype(int)
    Npts = np.size(ind)
    X    = np.zeros(Npts)
    X[ind] = Pos[:, 1]

    X_ind = np.where(np.abs(X - x) < 1e-6)
    Nz    = findNz(x, delta)
    print(X_ind[0].size)
    print(Nz * Ny)

    tids = np.arange(tid_str, tid_end, dt).astype(int)
    Nt   = np.size(tids)
    grad_rho_x = np.zeros((Nt, Ny, Nz))
    grad_rho_y = np.zeros((Nt, Ny, Nz))
    grad_rho_z = np.zeros((Nt, Ny, Nz))
    u          = np.zeros((Nt, Ny, Nz))
    v          = np.zeros((Nt, Ny, Nz))
    w          = np.zeros((Nt, Ny, Nz))
    p          = np.zeros((Nt, Ny, Nz))
    rho        = np.zeros((Nt, Ny, Nz))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        print("active threads: ", executor._max_workers)
        futures = {executor.submit(_load_timestep_grad, DataName_fmt, tid, ind, X_ind, Ny, Nz): i
                   for i, tid in enumerate(tids)}
        for future in tqdm(as_completed(futures), total=Nt, desc=f"x={x:.3f}"):
            i = futures[future]
            grad_rho_x[i], grad_rho_y[i], grad_rho_z[i], u[i], v[i], w[i], p[i], rho[i] = future.result()

    print(u.size)
    print(Nt * Ny * Nz)
    mf = u*rho
    return grad_rho_x, grad_rho_y, grad_rho_z, u, v, w, p, rho,mf


def write_time_history_h5(u, v, w, p, rho, mf, x, tid_str, tid_end, dt,
                           caseName:str, out_dir="./", max_workers=None):
    """
    Write the full time history of a single axial-profile plane to HDF5.

    Parameters
    ----------
    u, v, w, p, rho, mf : ndarray, shape (Nt, Ny, Nz)
        Direct output of extract_data for the plane at x.
    x : float
        Streamwise coordinate of the plane.
    tid_str, tid_end : int
        First and last (exclusive) timestep indices passed to extract_data.
    dt : int
        Timestep stride used in extract_data.
    out_dir : str
        Directory in which to create the HDF5 file.
    max_workers : int or None
        Thread count for parallel dataset writes (None = os.cpu_count()).

    Returns
    -------
    str
        Absolute path to the written HDF5 file.
    """
    tids  = np.arange(tid_str, tid_end, dt)
    z_max = z_lim(x)
    z     = np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta)

    x_tag = f"{x:.4f}".replace("-", "n").replace(".", "p")
    fname = os.path.join(out_dir, f"{caseName}_noz_int_time_history_x_{x_tag}.h5")

    flow_vars = {"u": u, "v": v, "w": w, "p": p, "rho": rho, "mf": mf}

    with h5py.File(fname, "w") as hf:
        # metadata and coordinates — written sequentially (structural ops)
        hf.attrs["x"]       = x
        hf.attrs["tid_str"] = tid_str
        hf.attrs["tid_end"] = tid_end
        hf.attrs["dt"]      = dt
        hf.create_dataset("tids", data=tids)
        hf.create_dataset("y",    data=y)
        hf.create_dataset("z",    data=z)

        # pre-allocate all flow-variable datasets before spawning threads
        for name, arr in flow_vars.items():
            hf.create_dataset(name, shape=arr.shape, dtype=arr.dtype,
                              compression="gzip", compression_opts=4)

        def _write_var(name):
            hf[name][...] = flow_vars[name]
            return name

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_write_var, name): name for name in flow_vars}
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc=f"writing h5 x={x:.3f}"):
                fut.result()

    print(f"  wrote: {fname}")
    return fname


def write_time_history_h5_grad(grad_rho_x, grad_rho_y, grad_rho_z,
                                u, v, w, p, rho, mf,
                                x, tid_str, tid_end, dt,
                                caseName: str, out_dir="./", max_workers=None):
    """
    Write the full time history (grad_rho_* + flow vars) for a single axial
    plane to HDF5.  Direct counterpart to write_time_history_h5 for data
    produced by extract_data_grad.

    Parameters
    ----------
    grad_rho_x, grad_rho_y, grad_rho_z : ndarray, shape (Nt, Ny, Nz)
    u, v, w, p, rho, mf               : ndarray, shape (Nt, Ny, Nz)
    x          : float  — streamwise coordinate of the plane
    tid_str, tid_end : int — timestep range (same as passed to extract_data_grad)
    dt         : int   — timestep stride
    caseName   : str   — prepended to the filename
    out_dir    : str   — output directory
    max_workers: int | None — thread count for parallel dataset writes

    Returns
    -------
    str  — absolute path of the written HDF5 file
    """
    tids  = np.arange(tid_str, tid_end, dt).astype(int)
    z_max = z_lim(x)
    z     = np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta)

    x_tag = f"{x:.4f}".replace("-", "n").replace(".", "p")
    fname = os.path.join(out_dir, f"{caseName}_noz_int_grad_time_history_x_{x_tag}.h5")

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
        hf.attrs["x"]       = x
        hf.attrs["tid_str"] = tid_str
        hf.attrs["tid_end"] = tid_end
        hf.attrs["dt"]      = dt
        hf.create_dataset("tids", data=tids)
        hf.create_dataset("y",    data=y)
        hf.create_dataset("z",    data=z)

        for name, arr in flow_vars.items():
            hf.create_dataset(name, shape=arr.shape, dtype=arr.dtype,
                              compression="gzip", compression_opts=4)

        def _write_var(name):
            hf[name][...] = flow_vars[name]
            return name

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(_write_var, name): name for name in flow_vars}
            for fut in tqdm(as_completed(futures), total=len(futures),
                            desc=f"writing h5 grad x={x:.3f}"):
                fut.result()

    print(f"  wrote: {fname}")
    return fname


def read_time_history_h5_grad(fname):
    """
    Read an HDF5 file written by write_time_history_h5_grad.

    Parameters
    ----------
    fname : str
        Path to the HDF5 file.

    Returns
    -------
    dict with keys:
        x, tid_str, tid_end, dt          : scalar metadata (from file attrs)
        tids, y, z                       : 1-D coordinate arrays
        grad_rho_x, grad_rho_y, grad_rho_z : ndarray, shape (Nt, Ny, Nz)
        u, v, w, p, rho, mf             : ndarray, shape (Nt, Ny, Nz)
    """
    with h5py.File(fname, "r") as hf:
        out = {
            "x":       hf.attrs["x"],
            "tid_str": hf.attrs["tid_str"],
            "tid_end": hf.attrs["tid_end"],
            "dt":      hf.attrs["dt"],
            "tids":    hf["tids"][:],
            "y":       hf["y"][:],
            "z":       hf["z"][:],
        }
        for var in ("grad_rho_x", "grad_rho_y", "grad_rho_z",
                    "u", "v", "w", "p", "rho", "mf"):
            out[var] = hf[var][:]
    return out


#TODO: plot various profile data at the 4 different lines defined with in BL_probes and assmemble these quantities
#Avarilable quantities, by order are:
#u,v,w,p,rho
#we would like to at least assemble the following:
#avg(u) profile, avg(u'v'), avg(u'w'), 
def plotTurbProf_C(u,v,w,rho,p,xc,y,z, out_dir = "./",eps=delta/1.9, fig_offset=0):
    """
    takes in time history and the 1d y,z vector and produces the following profiles at given streamwise station x = xc
    u_avg over y centerline (major axis)
    u_avg over z centerline (minor axis)
    u'v' over y centerline (major axis)
    u'w' over z centerline (minor axis)
    tke over y centerline (major axis)
    tke over z centerline (minor axis)
    """
    #do averging across an region (see axial line profiles)
    #definte limits
    y_lim = 0.85  #+/-y_lim
    z_lim = 0.35  #+/-z_lim
    #find bounds    
    y_ind_lo = np.where(((-y_lim - eps) < y) & (y < (-y_lim + eps)))
    y_ind_hi = np.where((( y_lim - eps) < y) & (y < ( y_lim + eps)))
    z_ind_lo = np.where(((-z_lim - eps) < z) & (z < (-z_lim + eps)))
    z_ind_hi = np.where((( z_lim - eps) < z) & (z < ( z_lim + eps)))

    
    y_ind_lo = y_ind_lo[0][0]
    y_ind_hi = y_ind_hi[0][0]
    z_ind_lo = z_ind_lo[0][0]
    z_ind_hi = z_ind_hi[0][0]

    print(y_ind_lo)
    print(y_ind_hi) 
    print(z_ind_lo)
    print(z_ind_hi) 

    y_ind = np.arange(y_ind_lo, y_ind_hi+1)
    z_ind = np.arange(z_ind_lo, z_ind_hi+1)
    #preivous centerline implementation
    #y_ind = np.where(np.abs(y) < eps)
    #z_ind = np.where(np.abs(z) < eps)
    #y_ind = y_ind[0]
    #z_ind = z_ind[0]
    
    
    Nz = np.size(z)


    #debug (print loc)
    #print(y[y_ind])
    #print(z[z_ind])
    
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

    u_avg_y = np.squeeze(np.mean(u_avg_y,axis=1))
    u_avg_z = np.squeeze(np.mean(u_avg_z,axis=0))

    rho_avg_y = np.squeeze(rho_avg[:,z_ind])
    rho_avg_z = np.squeeze(rho_avg[y_ind,:])

    rho_avg_y = np.squeeze(np.mean(rho_avg_y,axis=1))
    rho_avg_z = np.squeeze(np.mean(rho_avg_z,axis=0))

    p_avg_y = np.squeeze(p_avg[:,z_ind])
    p_avg_z = np.squeeze(p_avg[y_ind,:])

    p_avg_y = np.squeeze(np.mean(p_avg_y,axis=1))
    p_avg_z = np.squeeze(np.mean(p_avg_z,axis=0))

    print("shape of rho_avg_y: ", np.shape(rho_avg_y))
    print("shape of p_avg_y: ", np.shape(p_avg_y))
    T_avg_y = p_avg_y/rho_avg_y
    T_avg_z = p_avg_z/rho_avg_z
    print("Shape of T_avg_y: ", np.shape(T_avg_y))
    mu_avg_y = muCalc(T_avg_y)
    mu_avg_z = muCalc(T_avg_z)

    print("shape of mu_avg_y: ",np.shape(mu_avg_y))
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
    ufvf_avg_y = np.mean(ufvf_avg_y,axis=1)

    ufwf_avg = np.mean(uf*wf, axis=0)
    ufwf_avg_z = np.squeeze(ufwf_avg[y_ind,:])
    ufwf_avg_z = np.mean(ufwf_avg_z,axis=0)

    PlotName_fmt = "{}_{}.png"
    PlotName_fmt = os.path.join(out_dir, PlotName_fmt)
    LogLogName_fmt = "{}_{}_loglog.png"
    LogLogName_fmt = os.path.join(out_dir, LogLogName_fmt)
    SemiLogXName_fmt = "{}_{}_semilogx.png"
    SemiLogXName_fmt = os.path.join(out_dir, SemiLogXName_fmt)
    ## y_plus calc
    y_p = _calcYplus(y+1.0,u_tau_y,mu_w_y/rhoy_w)
    z_p = _calcYplus(z+0.5,u_tau_z,mu_w_z/rhoz_w) 

    y_p = y_p
    print("first ten elements of y_p: ", y_p[:10])
    z_p = z_p
    print("first ten elements of z_p", z_p[:10])
    u_p_y = u_avg_y / u_tau_y    
    u_p_z = u_avg_z / u_tau_z    
    ## ploting happending below
    
    fig_count = 9
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

    ##ufvf_avg_y
    #fig_id = 2+ fig_count*fig_offset
    #var_label = "ufvf_avg"
    #side_label = "y"
    #PlotName = PlotName_fmt.format(var_label,side_label)
    #plotProfile(xc,y,ufvf_avg_y,r"$\overline{u'v'}$", side_label,PlotName,fig_id)

    ##ufwf_avg_z
    #fig_id = 3 + fig_count*fig_offset
    #var_label = "ufwf_avg"
    #side_label = "z"
    #PlotName = PlotName_fmt.format(var_label,side_label)
    #plotProfile(xc,z,ufwf_avg_z,r"$\overline{u'w'}$", side_label,PlotName,fig_id)

    #ufvf_avg_y
    fig_id = 4+ fig_count*fig_offset
    var_label = "ufvf_avg_div_u_tau_y_2"
    side_label = "y+"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotSemiLogXProfile(xc,ufvf_avg_y[:Ny//2]/(u_tau_y**2),y_p[:Ny//2],side_label, r"$\frac{\overline{u'v'}}{u_{\tau,y}^2}$",PlotName,fig_id)

    #ufwf_avg_z
    fig_id = 5 + fig_count*fig_offset
    var_label = "ufwf_avg_div_u_tau_z_2"
    side_label = "z+"
    PlotName = PlotName_fmt.format(var_label,side_label)
    plotSemiLogXProfile(xc,ufwf_avg_z[:Nz//2]/(u_tau_z**2),z_p[:Nz//2],side_label, r"$\frac{\overline{u'w'}}{u_{\tau,z}}$",PlotName,fig_id)


    #u_avg_y law of the wall
    fig_id = 6 + fig_count*fig_offset
    var_label = "u+"
    side_label = "y+"
    PlotName = SemiLogXName_fmt.format(var_label,side_label)
    plotSemiLogXProfile(xc,u_p_y[:Ny//2],y_p[:Ny//2],side_label, r"$\overline{u}_+$",PlotName,fig_id)

    #u_avg_z law of the wall
    fig_id = 7 + fig_count*fig_offset
    var_label = "u+"
    side_label = "z+"
    PlotName = SemiLogXName_fmt.format(var_label,side_label)
    plotSemiLogXProfile(xc,u_p_z[:Nz//2],z_p[:Nz//2],side_label,r"$\overline{u}_+$",PlotName,fig_id)




def CalcTurbProf(u, v, w, rho, p, xc, y, z, eps=delta/1.9):
    """
    Compute boundary-layer turbulence profiles at streamwise station xc.

    Mirrors plotTurbProf_C but returns all plotted quantities as a dict
    instead of saving figures.

    Returns
    -------
    dict with keys:
        y, u_avg_y           : mean-u profile on major (y) axis
        z, u_avg_z           : mean-u profile on minor (z) axis
        y_half, ufvf_avg_y_norm  : u'v'/u_tau_y^2, wall-side half-profile
        z_half, ufwf_avg_z_norm  : u'w'/u_tau_z^2, wall-side half-profile
        y_p, u_p_y           : law-of-the-wall profile (y side, wall units)
        z_p, u_p_z           : law-of-the-wall profile (z side, wall units)
        u_tau_y, u_tau_z     : friction velocities
        tauy_w, tauz_w       : wall shear stresses
        rhoy_w, rhoz_w       : wall densities
        mu_w_y,  mu_w_z      : wall dynamic viscosities
        rho_avg_y, rho_avg_z : spanwise-averaged mean density profiles
        T_avg_y,   T_avg_z   : spanwise-averaged mean temperature profiles
    """
    y_lim_val = 0.85
    z_lim_val = 0.35

    y_ind_lo = np.where(((-y_lim_val - eps) < y) & (y < (-y_lim_val + eps)))[0][0]
    y_ind_hi = np.where((( y_lim_val - eps) < y) & (y < ( y_lim_val + eps)))[0][0]
    z_ind_lo = np.where(((-z_lim_val - eps) < z) & (z < (-z_lim_val + eps)))[0][0]
    z_ind_hi = np.where((( z_lim_val - eps) < z) & (z < ( z_lim_val + eps)))[0][0]

    y_ind = np.arange(y_ind_lo, y_ind_hi + 1)
    z_ind = np.arange(z_ind_lo, z_ind_hi + 1)

    Ny_loc = np.size(y)
    Nz_loc = np.size(z)

    u_avg   = np.mean(u,   axis=0)
    v_avg   = np.mean(v,   axis=0)
    w_avg   = np.mean(w,   axis=0)
    rho_avg = np.mean(rho, axis=0)
    p_avg   = np.mean(p,   axis=0)

    uf = u - u_avg
    vf = v - v_avg
    wf = w - w_avg

    u_avg_y = np.squeeze(u_avg[:, z_ind])
    u_avg_y = np.squeeze(np.mean(u_avg_y, axis=1))
    u_avg_z = np.squeeze(u_avg[y_ind, :])
    u_avg_z = np.squeeze(np.mean(u_avg_z, axis=0))

    rho_avg_y = np.squeeze(rho_avg[:, z_ind])
    rho_avg_y = np.squeeze(np.mean(rho_avg_y, axis=1))
    rho_avg_z = np.squeeze(rho_avg[y_ind, :])
    rho_avg_z = np.squeeze(np.mean(rho_avg_z, axis=0))

    p_avg_y = np.squeeze(p_avg[:, z_ind])
    p_avg_y = np.squeeze(np.mean(p_avg_y, axis=1))
    p_avg_z = np.squeeze(p_avg[y_ind, :])
    p_avg_z = np.squeeze(np.mean(p_avg_z, axis=0))

    T_avg_y = p_avg_y / rho_avg_y
    T_avg_z = p_avg_z / rho_avg_z
    mu_avg_y = muCalc(T_avg_y)
    mu_avg_z = muCalc(T_avg_z)

    mu_w_y = (mu_avg_y[0] + mu_avg_y[-1]) / 2.0
    mu_w_z = (mu_avg_z[0] + mu_avg_z[-1]) / 2.0

    dy = np.abs(y[1] - y[0])
    dz = np.abs(z[1] - z[0])
    dudy_w = (-3*u_avg_y[0] + 4*u_avg_y[1] - u_avg_y[2]) / (2*dy)
    dudz_w = (-3*u_avg_z[0] + 4*u_avg_z[1] - u_avg_z[2]) / (2*dz)

    tauy_w = mu_w_y * dudy_w
    tauz_w = mu_w_z * dudz_w
    rhoy_w = (rho_avg_y[0] + rho_avg_y[-1]) / 2.0
    rhoz_w = (rho_avg_z[0] + rho_avg_z[-1]) / 2.0

    u_tau_y = np.sqrt(tauy_w / rhoy_w)
    u_tau_z = np.sqrt(tauz_w / rhoz_w)

    ufvf_avg   = np.mean(uf * vf, axis=0)
    ufvf_avg_y = np.squeeze(ufvf_avg[:, z_ind])
    ufvf_avg_y = np.mean(ufvf_avg_y, axis=1)

    ufwf_avg   = np.mean(uf * wf, axis=0)
    ufwf_avg_z = np.squeeze(ufwf_avg[y_ind, :])
    ufwf_avg_z = np.mean(ufwf_avg_z, axis=0)

    y_p = _calcYplus(y + 1.0, u_tau_y, mu_w_y / rhoy_w)
    z_p = _calcYplus(z + 0.5, u_tau_z, mu_w_z / rhoz_w)
    u_p_y = u_avg_y / u_tau_y
    u_p_z = u_avg_z / u_tau_z

    return {
        "y":       y,
        "u_avg_y": u_avg_y,
        "z":       z,
        "u_avg_z": u_avg_z,
        "y_half":          y[:Ny_loc // 2],
        "ufvf_avg_y_norm": ufvf_avg_y[:Ny_loc // 2] / u_tau_y**2,
        "z_half":          z[:Nz_loc // 2],
        "ufwf_avg_z_norm": ufwf_avg_z[:Nz_loc // 2] / u_tau_z**2,
        "y_p":   y_p[:Ny_loc // 2],
        "u_p_y": u_p_y[:Ny_loc // 2],
        "z_p":   z_p[:Nz_loc // 2],
        "u_p_z": u_p_z[:Nz_loc // 2],
        "u_tau_y": u_tau_y,
        "u_tau_z": u_tau_z,
        "tauy_w":  tauy_w,
        "tauz_w":  tauz_w,
        "rhoy_w":  rhoy_w,
        "rhoz_w":  rhoz_w,
        "mu_w_y":  mu_w_y,
        "mu_w_z":  mu_w_z,
        "rho_avg_y": rho_avg_y,
        "rho_avg_z": rho_avg_z,
        "T_avg_y":   T_avg_y,
        "T_avg_z":   T_avg_z,
    }


def _calcYplus(y,u_tau,nu):
    #Return coordinate in y+ units
    #INPUT: y: coordinate
    #       u_tau: friction velocity
    #       nu: kinematic viscosity
    #OUTPUT:  y in y+ units
    # y+ def: y+ = y * u_tau/nu 
    return y*u_tau/nu

def plotProfile(xc,pos,var,var_label:str,side_label:str,PlotName: str,fig_id,titName=None,lineLabel=None):
    plt.figure(fig_id)
    if lineLabel == None:
        plt.plot(var,pos, label = "x = %.2f" % xc)
    else:
        plt.plot(var,pos,label=lineLabel)
    plt.xlabel(var_label)
    plt.ylabel(side_label)
    if titName != None:
        plt.title(titName)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PlotName, dpi = 300, bbox_inches="tight")

def plotLogLogProfile(xc,pos,var,var_label:str,side_label:str,PlotName: str,fig_id,titName=None,lineLabel=None):
    plt.figure(fig_id)
    if lineLabel == None:
        plt.loglog(var,pos, label = "x = %.2f" % xc)
    else:
        plt.loglog(var,pos,label=lineLabel)
    plt.xlabel(var_label)
    plt.ylabel(side_label)
    if titName != None:
        plt.title(titName)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PlotName, dpi = 300, bbox_inches="tight")


def plotSemiLogXProfile(xc,pos,var,var_label:str,side_label:str,PlotName: str,fig_id, titName=None,lineLabel=None):
    plt.figure(fig_id)
    if lineLabel == None:
        plt.semilogx(var,pos, label = "x = %.2f" % xc)
    else:
        plt.semilogx(var,pos,label=lineLabel)
    plt.xlabel(var_label)
    plt.ylabel(side_label)
    if titName != None:
        plt.title(titName)
    plt.legend()
    plt.tight_layout()
    plt.savefig(PlotName, dpi = 300, bbox_inches="tight")

def evalMfAvg(mf,z,y):
    mf_avg = np.squeeze(np.mean(mf,axis=0))
    mf_tot = np.trapezoid(np.trapezoid(mf,z),y)
    mf_tot = np.mean(mf_tot) 
    return mf_tot

def plotAxialProf(xc,var,z,y,var_label:str,out_dir,var_lim=None):
    PlotName_fmt = "{}_x_{}.png"
    PlotName_fmt = os.path.join(out_dir,PlotName_fmt)
    PlotName = PlotName_fmt.format(var_label,xc) 
    plt.figure()
    if var_lim ==None:
        plt.pcolormesh(y,z,np.transpose(var),shading='nearest')
    else:
        #if var_lim is defined
        plt.pcolormesh(y,z,np.transpose(var),shading='nearest',vmin=var_lim[0],vmax=var_lim[1])
    plt.xlabel("y")
    plt.ylabel("z")
    plt.colorbar()
    plt.title(var_label)
    plt.tight_layout()
    plt.axis('equal')

    plt.axvline(x= 0.85,color='r',linestyle='--')
    plt.axvline(x=-0.85,color='r',linestyle='--')
    plt.axhline(y= 0.35,color='c',linestyle='-.')
    plt.axhline(y=-0.35,color='c',linestyle='-.')

    plt.savefig(PlotName,dpi=300,bbox_inches="tight")

def calcDel99(u_avg,y):
    #TODO: calculate delta_99
    """
    INPUT:
    u_avg: the averaged (both in time and in trans)
    y: normal direction
    OUTPUT:
    """
    pass

def calcDelStr():
    #TODO: calculate displacement thickness (delta_star)
    pass

def calcTheta():
    #TODO: calculate momentum thickness (theta)
    pass

def main():
    #these are for test purposes
    #General BL test cases
    #TBL_test_cases = ["0.0125", "0.01875", "0.025"]
    #tid_str_cases  = [40000, 40000, 130000]
    #tid_end_cases  = [100000, 100000, 190000]

    #BL test debug    
    #TBL_test_cases = ["0.025"]
    #tid_str_cases  = [155000]
    #tid_end_cases  = [160000]
    #data_dir_fmt = "/anvil/scratch/x-sdai/BL_test_baseline_{:s}/pcprobe_int_axprof"
    #out_dir_fmt = "/anvil/scratch/x-sdai/BL_post_proc/BL_{:s}"
    
    
    #TBL_test_cases = ["151M_akhil_restart","TBL_0.0125_157M"]
    #tid_str_cases = [529050]
    #tid_str_cases = [549600]
    #id_end_cases = [549900]
    data_dir_fmt = "/anvil/scratch/x-sdai/AR2_base_{:s}/pcprobes_noz_int"
    out_dir = "/anvil/scratch/x-sdai/AIAA_SciTech_2027/BL_comp/"
    dt = 150
    debug_flag = False



    TBL_test_cases = ["151M_akhil_restart","TBL_0.0125_157M"]
    Plot_labels=["Baseline 151M", "BaseRough 157M"]


    if debug_flag:
        #launched on login-node
        cpus = 4
    else:
        #launched on compute node
        cpus = int(os.getenv("SLURM_NTASKS"))
    print("assigned core count: ", cpus)
    #DEBUG
    print(data_dir_fmt.format(TBL_test_cases[0]))
    print(out_dir)
    print("targeting postions are: ", xs)
    # read data and get get plots
    for x in xs: 
        #for each axial section
        #line plots of the following:
        # u_avg_y vs y        -> fig0
        x_tag = f"{x:.4f}".replace("-", "n").replace(".", "p")
        figName0 = f"u_avg_y_x_{x_tag}.png"
        titName0 = rf"$\overline{{u}}, x={x_tag}$"
        figName0 = os.path.join(out_dir,figName0)
        # u+_y vs y+          -> fig1
        figName1 = f"up_avg_yp_x_{x_tag}.png"
        titName1 = rf"$\overline{{u^+}}, x={x_tag}$"
        figName1 = os.path.join(out_dir,figName1)
        # <u'v'>/u_tau^2 vs y+-> fig2
        figName2 = "upvp_avg_y_x_{x_tag}.png"
        titName2 = rf"$\frac{{\overline{{u'v'}}}}{{u_{{\tau,y}}^2}}, x={x}$"
        figName2 = os.path.join(out_dir,figName2)
        
        # u_avg_z vs z        -> fig3
        figName3 = f"u_avg_z_x_{x_tag}.png"
        titName3 = rf"$\overline{{u}}, x={x}$"
        figName3 = os.path.join(out_dir,figName3)
        # u+_z vs z+          -> fig4
        figName4 = f"up_avg_zp_x_{x_tag}.png"
        titName4 = rf"$\overline{{u^+}}, x={x}$"
        figName4 = os.path.join(out_dir,figName4)
        # <u'w'>/u_tau^2 vs z+-> fig5
        figName5 = f"upvp_avg_z_x_{x_tag}.png"
        titName5 = rf"$\frac{{\overline{{u'w'}}}}{{u_{{\tau,z}}^2}}, x={x}$"
        figName5 = os.path.join(out_dir,figName5)

        for i,caseName in enumerate(TBL_test_cases):
            fname = os.path.join(out_dir, f"{caseName}_noz_int_grad_time_history_x_{x_tag}.h5")
            #extract data for a given case
            xcData = read_time_history_h5_grad(fname)
            #positional data
            y = xcData.get("y")
            z = xcData.get("z")
            #variable (Nt * Ny * Nz)
            grad_rho_x=xcData.get("grad_rho_x")
            u = xcData.get("u")
            v = xcData.get("v")
            w = xcData.get("w")
            rho = xcData.get("rho")
            p = xcData.get("p")

            #get turbulence data
            turbStat = CalcTurbProf(u, v, w, rho, p, x, y, z)
            u_avg_y         = turbStat.get("u_avg_y")
            y_p             = turbStat.get("y_p")
            u_p_y           = turbStat.get("u_p_y")
            ufvf_avg_y_norm = turbStat.get("ufvf_avg_y_norm")
            u_avg_z         = turbStat.get("u_avg_z")
            z_p             = turbStat.get("z_p")
            u_p_z           = turbStat.get("u_p_z")
            ufwf_avg_z_norm = turbStat.get("ufwf_avg_z_norm")

            #plot 2D axial profiles
            plotProfile(x,u_avg_y,y,r"$y$",r"$\overline{u}/c_\infty$",figName0,0,titName0,Plot_labels[i])
            plotSemiLogXProfile(x,u_p_y,y_p,r"$y+$", r"$\overline{u}_+$",figName1,1,titName1,Plot_labels[i])
            plotProfile(x,np.abs(ufvf_avg_y_norm),y_p,r"$y+$",rf"$\frac{{\overline{{u'v'}}}}{{u_{{\tau,y}}^2}}$",figName2,2,titName2,Plot_labels[i])
            plotProfile(x,u_avg_z,z,r"$z$",r"$\overline{u}/c_\infty$",figName3,3,titName3,Plot_labels[i])
            plotSemiLogXProfile(x,u_p_z,z_p,r"$z+$", r"$\overline{u}_+$",figName4,4,titName4,Plot_labels[i])
            plotProfile(x,np.abs(ufwf_avg_z_norm),z_p,r"$z+$",rf"$\frac{{\overline{{u'v'}}}}{{u_{{\tau,y}}^2}}$",figName5,5,titName5,Plot_labels[i])
                        
        plt.close('all') #close all per axial station
 
    """
    #for i, case in enumerate(TBL_test_cases):
    #    plt.close('all')
    #    print("Now processing for TBL case: ", case)
    #    data_dir = data_dir_fmt.format(case)
    #    print(xs)
    #    z_max = z_lim(0)
    #    z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
    #    #out_dir = out_dir_fmt.format(case)
    #    #posName = os.path.join(data_dir,"int_axprof.pxyz")
    #    #fname_fmt = "int_axprof.{:08d}.pcd"
    #    posName = os.path.join(data_dir,"nov_int.pxyz")
    #    fname_fmt = "nov_int.{:08d}.pcd"
    #    DataName_fmt = os.path.join(data_dir,fname_fmt)
    #    tid_str = tid_str_cases[i]
    #    tid_end = tid_end_cases[i]
    #    tids = np.arange(tid_str,tid_end+dt,dt)
    #    print(np.size(tids))
    #    for x in xs:
    #        #cpus = int(os.environ.get('SLURM_NTASK', 1))
    #        grad_rho_x,grad_rho_y,grad_rho_z,u,v,w,p,rho,mf= extract_data_grad(x,posName,DataName_fmt,tid_str,tid_end,dt,max_workers=cpus)
    #        z_max = z_lim(x)
    #        #print(y)
    #        T = p/rho
    #        c = np.sqrt(gamma*T)
    #        z =np.arange(-(z_max - delta / 2.0), (z_max - delta / 2.0) + delta * 0.5, delta) 
    #        print("mass flow at x=%.2f is mf = %.2f " % (x,evalMfAvg(mf,z,y)))
    #        #plotTurbProf_C(u,v,w,rho,p,x,y,z,out_dir = out_dir)
    #        #TODO: write h5 data to data storage
    #        write_time_history_h5_grad(grad_rho_x,grad_rho_y,grad_rho_z,u,v,w,p,rho,mf,x,tid_str,tid_end,dt,case,out_dir)
    #        plotAxialProf(x,np.mean(u,axis=0),z,y,"u_avg",out_dir)
    #        plotAxialProf(x,np.mean(p,axis=0),z,y,"p_avg",out_dir)
    #        plotAxialProf(x,np.mean(u,axis=0)/np.mean(c,axis=0),z,y,"Mach_u",out_dir,var_lim=(0,1.5))

            

    #load baseline h5 make the same plots based using the same functions
    #need to supply new x,y,z
    #grab x slices
    """
if __name__ == "__main__":
    t_start = time.perf_counter()
    main()
    t_end = time.perf_counter()
    print(f"Total execution time: {t_end - t_start:.2f} s")



