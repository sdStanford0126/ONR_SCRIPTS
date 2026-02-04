
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
import os
import sys
from Universal_Subroutines import setPlotpref
setPlotpref() 
from tqdm import tqdm
import python_general_util_func as gutils
import scipy.signal

input_dir  = "/Users/steven/Library/CloudStorage/OneDrive-Personal/Stanford/ONR project/results/baseline_new/144M/pcprobes_turb/"
output_dir = "/Users/steven/Library/CloudStorage/OneDrive-Personal/Stanford/ONR project/results/baseline_new/144M/turb_probes_post/"

def readProbes(input_dir, fname_prefix, startTime,endTime,step):
    """
    #from the input file and directory, read Charles Probe output and organized them based on the original order in the specified file
    #TODO: read in .pxyz file to get the location information and the indexes information, reorganized based on the order
    #should be a 2D array (# probes being the row, x,y,z corresponding to the columns)
    """
    pos_fname = fname_prefix+".pxyz"
    #load position file as needed
    prob_pos = np.loadtxt(input_dir+pos_fname, skiprows=1)
    #print(prob_pos[:,2]) # debug 
    numProb,_=np.shape(prob_pos)
    inds = prob_pos[:,0]
    inds = inds.astype(int)
    

    #print(inds) # debug
    xyz_pos = np.zeros(shape=(numProb,3))
    xyz_pos[inds,:] = prob_pos[:,1:]
    #debug block
    #     #print(xyz_pos[:,0])
    #print(xyz_pos[:,1])
    #print(xyz_pos[:,2])
    #TODO: read in the data files
    times = np.arange(startTime,endTime+step,step)
    print(times[-1])
    L_hist = np.size(times) #length of the time record
    data_fmt = fname_prefix + ".{:08d}.pcd"
    #debug
    print(data_fmt.format(168600))
    data_cur = np.loadtxt(input_dir + data_fmt.format(168600),skiprows=1)
    u_test = np.zeros((numProb))
    u_test[inds] = data_cur[:,0]
    #print(data_cur[:,0])
    #print(u_test)    
    #print(inds)
    """
    for now, it is assume that the data file for the pointcloud probes
    have the following structure: (by columns)
    u, v, w, p, rho, T
    """
    u_hist   = np.zeros((numProb,L_hist),dtype = np.float64)
    v_hist   = np.zeros((numProb,L_hist),dtype = np.float64)
    w_hist   = np.zeros((numProb,L_hist),dtype = np.float64)
    p_hist   = np.zeros((numProb,L_hist),dtype = np.float64)
    rho_hist = np.zeros((numProb,L_hist),dtype = np.float64)
    T_hist   = np.zeros((numProb,L_hist),dtype = np.float64)

    print(np.shape(u_hist))
    index = 0
    for time in tqdm(times):
        #load data at timestep time
        
        data_name = data_fmt.format(time)
        data_cur =  np.loadtxt(input_dir + data_name,skiprows=1) #load data at this time step
        #write to time history with the corrected order
        u_hist[inds,index]   = data_cur[:,0]
        v_hist[inds,index]   = data_cur[:,1]
        w_hist[inds,index]   = data_cur[:,2]
        p_hist[inds,index]   = data_cur[:,3]
        rho_hist[inds,index] = data_cur[:,4]
        T_hist[inds,index]   = data_cur[:,5]
        percent = int(((index+1)/L_hist)*100)
        index +=1
        #if percent%10 == 0:
        #    print(f"data read completed at {percent}%")
        #if index == 0:
        #    print(u_hist[:,index])
    #return the results
    vars_hist = {
        "u": u_hist,
        "v": v_hist,
        "w": w_hist,
        "p": p_hist,
        "rho": rho_hist,
        "T": T_hist,
    }
    return xyz_pos, vars_hist

def calcFFT_Welch(var_hist,avg_flag,df,step,overlap=0.5,dt = 1e-3):
    #TODO: implement Welch method for spectral convergence with FFT
    """
    input:
    var_hist: variable time history (mean removed!)
    avg_flag: spatial averaging (True:do averaging, False: no averaging)
    df: frequency resolution
    step: sampling timestep size (1/f_samp)
    overlap: how much overlap between neighboring blocks
    output:
    Welch Method based spectrum (freq and spectra) for the given variable!
    """
    
    #TODO: determine block size
    T_b = 1/df #block sample time length
    t_samp = step*dt; 
    Nsamp_b = T_b//t_samp
    df_act = 1/(Nsamp_b*t_samp)
    print(f"frequency resolution is {df_act}\n") 
    Nprobs, Nsamp_tot = np.shape(var_hist)
    #block_step = int((1-overlap)*Nsamp_b)
    """
    block_str_inds = []
    block_end_inds = []
    if Nsamp_tot > Nsamp_b:
        start_ind = 0
        end_ind = int(Nsamp_b -1)
        block_num = 0
        print(end_ind)
        print(block_step)
        while start_ind < Nsamp_tot:
            if end_ind > Nsamp_tot-1:
                break #truncate the timeseries for now
            block_str_inds.append(start_ind)
            block_end_inds.append(end_ind)
            start_ind += block_step
            end_ind += block_step
            block_num += 1
            print(start_ind)
            print(end_ind)
        print(f"number of blocks is {block_num}")
    else:
        print(f"time record too short for specified freq. res.")
        print(f"treating this as a single block")
        df_act = 1/(Nsamp_tot*step)
        print(f"actual frequency resolution is {df_act}\n")
        block_str_inds.append(0)
        block_end_inds.append(Nsamp_tot - 1) 
    """
    #TODO: divide to blocks and do the math
    if Nsamp_tot <= Nsamp_b:
        freq,var_spec=scipy.signal.welch(var_hist,fs = 1/t_samp,nperseg=Nsamp_tot,noverlap = 0,axis = 1 )
    else:
        freq,var_spec = scipy.signal.welch(var_hist,fs = 1/t_samp,window='hann',nperseg=Nsamp_b,noverlap=int(Nsamp_b*overlap),axis = 1)
    #TODO: average between the probes that should be symmetric
    if avg_flag:
        #we are doing averaging here
        assert(Nprobs%2 == 0) #needs to be even number
        var_spec_avg = (var_spec[:Nprobs//2,:] + var_spec[Nprobs//2:,:])/2
        var_spec = var_spec_avg
    return freq,var_spec  

def plot_spectra(St, var_spec,x_pos_probes,fig_name_fmt,varName, SpecColors = matplotlib.color_sequences['tab20']):
    """
    Input: 
        St: Strouhal number range (1D vector)
        var_spec: spectrum 2D array (number of unique probes x spectrum length), spectrum length = length(St)
        x_pos_probes: the x position of the unique probes
        fig_name_fmt: format for the figure names
        varName:  variable name
        colors: the discrete color map to be used
    
    Output:
        Plots of St vs. var_spectrum, each probe get its own plot    
    """
    Nprobs, L_spec = np.shape(var_spec)
    y_min = max(np.min(var_spec),1e-8)
    y_max = max(np.max(var_spec),1e-1)
    assert np.size(St) == L_spec, "St vector and the spectrum does NOT have compatible length"
    for i in range(Nprobs):
        plt.figure()
        plt.tight_layout()
        plt.xlabel(rf"$St$")
        plt.ylabel(rf"${varName}$")
        plt.loglog(St,var_spec[i,:],color = SpecColors[i])
        plt.ylim([y_min, y_max])
        x_pos = x_pos_probes[i]
        title_name_fmt = "{:s} PSD, vs $St$,$x/h=$ {:.02f}"
        plt.title(title_name_fmt.format(varName,x_pos))
        fig_name = output_dir + fig_name_fmt.format(x_pos)
        plt.savefig(fig_name,bbox_inches='tight', dpi = 300)
        plt.close()

def plot_spectra_combined(St, var_spec,x_pos_probes,fig1_name_fmt,varName, SpecColors = matplotlib.color_sequences['tab20'],x_step=1):
    Nprobs, L_spec = np.shape(var_spec)
    y_min = max(np.min(var_spec),1e-8)
    y_max = np.max(var_spec)
    assert np.size(St) == L_spec, "St vector and the spectrum does NOT have compatible length"
    plt.figure()
    plt.tight_layout()
    plt.xlabel(rf"$St$")
    plt.ylabel(rf"${varName}$")
    plt.ylim([y_min, y_max])
    legend_fmt = "$x/h= ${:.02f}"
    title_fmt = "{:s} PSD vs St at multiple x position"
    for i in np.arange(0,Nprobs,x_step):
        x_pos = x_pos_probes[i]
        plt.loglog(St,var_spec[i,:],color = SpecColors[i], label=legend_fmt.format(x_pos))
    plt.title(title_fmt.format(varName))
    fig_name = output_dir + fig1_name_fmt
    plt.legend()
    plt.savefig(fig_name,bbox_inches='tight', dpi = 300)
    plt.close()
def main():
    #pre-amble to the data
    print(f"starting Turbulence Probe Processing!\n")
    print(f"input directory is {input_dir}\n")
    print(f"Output_dir is: {output_dir}\n")

    sides = ["xy", "xz"]
    for side in sides:
        #side = "xy" #"xy" or "xz"
        if side == "xy":
            fname_prefix = "turb_xy"
        else:
            fname_prefix = "turb_xz"
        print(f"working on side {side}\n")
        startTime = 168600
        endTime = 292300
        step = 100

        xyz_pos,vars_hist=readProbes(input_dir, fname_prefix, startTime,endTime,step)
        u_hist = vars_hist["u"]
        u_mean = np.mean(u_hist,axis=1)
        #print(xyz_pos[:,0])
        #print(u_mean)

        up_hist = np.transpose(np.transpose(u_hist) - u_mean)
        #print(up_hist[:,0])
        df = 1/75
        avg_flag = True
        overlap = 0.75
        dt = 1e-3
        freq, up_spec = calcFFT_Welch(up_hist,avg_flag,df,step,overlap,dt)
        NPR = 3
        Ma = gutils.computeMa(NPR)
        h = 1
        L = 2
        Deq_s = gutils.computeDeq(h,L)
        St_spec= gutils.calcSt_sim(freq,Ma,Deq_s)
        #print(St_spec)
        plt.figure()
        SpecColors = matplotlib.color_sequences['tab20b']
        if side == "xy":
            fig_name_fmt = "up_spec_144M_xy_x{:.02f}_.pdf"
            fig1_name_fmt = "up_spec_144M_xy_x_all_.pdf"
            x_step = 3
        else:
            fig_name_fmt = "up_spec_144M_xz_x{:.02f}_.pdf"
            fig1_name_fmt = "up_spec_144M_xz_x_all_.pdf"
            x_step = 2
        varName = "u'"
        x_pos_probes = xyz_pos[:,0]
        plot_spectra(St_spec, up_spec,x_pos_probes,fig_name_fmt,varName, SpecColors)
    
        plot_spectra_combined(St_spec, up_spec,x_pos_probes,fig1_name_fmt,varName, SpecColors,x_step)
    

        p_hist = vars_hist["p"]
        p_mean = np.mean(p_hist, axis = 1)
        pp_hist = np.transpose(np.transpose(p_hist) - p_mean)
        freq, pp_spec = calcFFT_Welch(pp_hist,avg_flag, df, step, overlap, dt)
        St_spec = gutils.calcSt_sim(freq,Ma,Deq_s)
        if side == "xy":
            fig_name_fmt = "pp_spec_144M_xy_x{:.02f}_.pdf"
            fig1_name_fmt = "pp_spec_144M_xy_x_all_.pdf"
            x_step = 3
        else:
            fig_name_fmt = "pp_spec_144M_xz_x{:.02f}_.pdf"
            fig1_name_fmt = "pp_spec_144M_xz_x_all_.pdf"
            x_step = 2
        varName = "p'"
        plot_spectra(St_spec, pp_spec,x_pos_probes,fig_name_fmt,varName, SpecColors)
        plot_spectra_combined(St_spec, pp_spec,x_pos_probes,fig1_name_fmt,varName, SpecColors, x_step)
        """
        Nprobs_unique,_ = np.shape(up_spec) 
        # test print section 
        for i in range(Nprobs_unique):
            plt.semilogx(St_spec,up_spec[i,:])
        plt.xlabel("St_spec")
        plt.ylabel("Spectra")
        plt.show()
        #print(u_hist[0,:])
        """
if __name__ == "__main__":
    main()