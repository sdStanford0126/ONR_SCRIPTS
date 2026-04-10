import os
import numpy as np
import matplotlib.pyplot as plt
from Universal_Subroutines import setPlotpref
setPlotpref()

Ny = 400
Nz = 200

def readprobes(pos_name):
    """
    Read probe position file (.pxyz) and split into 4 lines:
    1) y-line: x=0, z=0
    2) y-line: x=0, z=0.45
    3) z-line: x=0, y=0
    4) z-line: x=0, y=0.95
    Returns:
        xyz_pos: (N,3) probe positions ordered by probe index
        lines: list of 4 arrays, each (Ni,3) sorted along the varying axis
        line_inds: list of 4 arrays, indices into xyz_pos for each line
    """
    pos_data = np.loadtxt(pos_name, skiprows=1)
    if pos_data.ndim != 2 or pos_data.shape[1] < 4:
        raise ValueError(f"Position file {pos_name} must have at least 4 columns: index,x,y,z")

    inds = pos_data[:, 0].astype(int)
    num_prob = pos_data.shape[0]
    xyz_pos = np.zeros((num_prob, 3), dtype=float)
    xyz_pos[inds, :] = pos_data[:, 1:4]

    x = xyz_pos[:, 0]
    y = xyz_pos[:, 1]
    z = xyz_pos[:, 2]
    tol = 1e-6

    def _select(mask, sort_axis):
        idx = np.where(mask)[0]
        if idx.size == 0:
            return idx, np.empty((0, 3))
        coords = xyz_pos[idx]
        order = np.argsort(coords[:, sort_axis])
        return idx[order], coords[order]

    line1_idx, line1_xyz = _select(np.isclose(x, 0.0, atol=tol) & np.isclose(z, 0.0, atol=tol), 1)
    line2_idx, line2_xyz = _select(np.isclose(x, 0.0, atol=tol) & np.isclose(z, 0.45, atol=tol), 1)
    line3_idx, line3_xyz = _select(np.isclose(x, 0.0, atol=tol) & np.isclose(y, 0.0, atol=tol), 2)
    line4_idx, line4_xyz = _select(np.isclose(x, 0.0, atol=tol) & np.isclose(y, 0.95, atol=tol), 2)


    lines = [line1_xyz, line2_xyz, line3_xyz, line4_xyz]
    print(line1_xyz)
    line_inds = [line1_idx, line2_idx, line3_idx, line4_idx]

    return inds, xyz_pos, lines, line_inds

def readData(fname,line_inds,inds):
    #reading data from file
    data = np.loadtxt(fname,skiprows=1)
    num_prob = data.shape[0]
    num_var  = data.shape[1]
    data_org = np.zeros((num_prob,num_var),dtype=float)
    data_org[inds,:] = data
    #the data are
    #avg(p), avg(rho), comp(avg(u),0), comp(avg(u),1), comp(avg(u),2), 
    #comp(rms(u),0), comp(rms(u),1), comp(rms(u),2)
    #grab only the streamwise velocity and rms data for now
    u_avg_x = data_org[:,2]
    u_rms_x = data_org[:,5]
    u_rms_y = data_org[:,6]
    u_rms_z = data_org[:,7]

    u_avg_x1 = u_avg_x[line_inds[0]]
    u_avg_x2 = u_avg_x[line_inds[1]]
    u_avg_x3 = u_avg_x[line_inds[2]]
    u_avg_x4 = u_avg_x[line_inds[3]]
    """
    print(line_inds[0])
    plt.figure()
    plt.plot(u_avg_x1)
    plt.show()
    """
    u_rms_x1 = u_rms_x[line_inds[0]]
    u_rms_x2 = u_rms_x[line_inds[1]]
    u_rms_x3 = u_rms_x[line_inds[2]]
    u_rms_x4 = u_rms_x[line_inds[3]]

    u_rms_y1 = u_rms_y[line_inds[0]]
    u_rms_y2 = u_rms_y[line_inds[1]]
    u_rms_y3 = u_rms_y[line_inds[2]]
    u_rms_y4 = u_rms_y[line_inds[3]]

    u_rms_z1 = u_rms_z[line_inds[0]]
    u_rms_z2 = u_rms_z[line_inds[1]]
    u_rms_z3 = u_rms_z[line_inds[2]]
    u_rms_z4 = u_rms_z[line_inds[3]]
    
    u_avg_x_lines =[u_avg_x1,u_avg_x2,u_avg_x3,u_avg_x4]
    u_rms_x_lines =[u_rms_x1,u_rms_x2,u_rms_x3,u_rms_x4]
    u_rms_y_lines =[u_rms_y1,u_rms_y2,u_rms_y3,u_rms_y4]
    u_rms_z_lines =[u_rms_z1,u_rms_z2,u_rms_z3,u_rms_z4]

    return u_avg_x_lines,u_rms_x_lines,u_rms_y_lines,u_rms_z_lines

def readData1(fname,line_inds,inds):
    #reading data from file
    data = np.loadtxt(fname,skiprows=1)
    num_prob = data.shape[0]
    num_var  = data.shape[1]
    data_org = np.zeros((num_prob,num_var),dtype=float)
    data_org[inds,:] = data
    #the data are
    #avg(p), avg(rho), comp(avg(u),0), comp(avg(u),1), comp(avg(u),2), 
    #comp(rms(u),0), comp(rms(u),1), comp(rms(u),2)
    #grab only the streamwise velocity and rms data for now
    u_avg_x = data_org[:,0]
    u_rms_x = data_org[:,1]
    u_rms_y = data_org[:,2]
    u_rms_z = data_org[:,3]

    u_avg_x1 = u_avg_x[line_inds[0]]
    u_avg_x2 = u_avg_x[line_inds[1]]
    u_avg_x3 = u_avg_x[line_inds[2]]
    u_avg_x4 = u_avg_x[line_inds[3]]
    """
    print(line_inds[0])
    plt.figure()
    plt.plot(u_avg_x1)
    plt.show()
    """
    u_rms_x1 = u_rms_x[line_inds[0]]
    u_rms_x2 = u_rms_x[line_inds[1]]
    u_rms_x3 = u_rms_x[line_inds[2]]
    u_rms_x4 = u_rms_x[line_inds[3]]

    u_rms_y1 = u_rms_y[line_inds[0]]
    u_rms_y2 = u_rms_y[line_inds[1]]
    u_rms_y3 = u_rms_y[line_inds[2]]
    u_rms_y4 = u_rms_y[line_inds[3]]

    u_rms_z1 = u_rms_z[line_inds[0]]
    u_rms_z2 = u_rms_z[line_inds[1]]
    u_rms_z3 = u_rms_z[line_inds[2]]
    u_rms_z4 = u_rms_z[line_inds[3]]
    
    u_avg_x_lines =[u_avg_x1,u_avg_x2,u_avg_x3,u_avg_x4]
    u_rms_x_lines =[u_rms_x1,u_rms_x2,u_rms_x3,u_rms_x4]
    u_rms_y_lines =[u_rms_y1,u_rms_y2,u_rms_y3,u_rms_y4]
    u_rms_z_lines =[u_rms_z1,u_rms_z2,u_rms_z3,u_rms_z4]

    return u_avg_x_lines,u_rms_x_lines,u_rms_y_lines,u_rms_z_lines


def plotQuantities(lines,var_lines,var_name,label="",prev_fig=0):
    #plot for all four lines the given quantity
    #for line 1/2
    ind_offset = prev_fig * len(lines)
    print("offset level", ind_offset)
    for ind, line_xyz in enumerate(lines):
        x = line_xyz[0,0]
        print(x)
        var_line = var_lines[ind]
        if ind <2: #first two
            z = line_xyz[0,2]
            y = line_xyz[:,1]
            plt.figure(ind + ind_offset)
            plt.tight_layout()
            title_name_fmt = "{:s}, at x/h={:.02f}, z/h={:.02f}"
            title_name = title_name_fmt.format(var_name,x,z)
            print("title_name: ", title_name)
            fig_name_fmt = "x_{:.02f}_line{:d}_{:s}.pdf"
            fig_name = fig_name_fmt.format(x,ind,var_name)
            print("fig_name: ", fig_name)
            plt.title(title_name)
            plt.plot(var_line,y,"-o",label=label)
            plt.xlabel(var_name)
            plt.ylabel("y")
            plt.legend()
            plt.savefig(os.path.join(outdir,fig_name),bbox_inches='tight', dpi = 300)
        else:
            z = line_xyz[:,2]
            y = line_xyz[0,1]
            plt.figure(ind+ind_offset)
            plt.tight_layout()
            title_name_fmt = "{:s}, at $x/h={:.02f}, y/h=${:.02f}"
            title_name = title_name_fmt.format(var_name,x,y)
            fig_name_fmt = "x_{:.02f}_line{:d}_{:s}.pdf"
            fig_name = fig_name_fmt.format(x,ind,var_name)
            plt.title(title_name)
            plt.plot(var_line,z,"-o",label=label)
            plt.xlabel(var_name)
            plt.ylabel("z")
            plt.legend()
            plt.savefig(os.path.join(outdir,fig_name),bbox_inches='tight', dpi = 300)


def main():
    
    
    #I/O settings
    indir_fmt = "/pcprobes_{:s}_avg"
    probename = "BL"
    fname_fmt = "{:s}.{:08d}.pcd"
    pos_name_fmt = "{:s}.pxyz"
    
    casedir1 = "/Users/steven/OneDrive/Stanford/ONR project/results/BL_test/smooth"
    casedir2 = "/Users/steven/OneDrive/Stanford/ONR project/results/BL_test/0_025"
    casedir3 = "/Users/steven/OneDrive/Stanford/ONR project/results/BL_test/0_0125"
    casedir4 = "/Users/steven/OneDrive/Stanford/ONR project/results/BL_test/0_01875"
    
    indir = indir_fmt.format(probename)
    inputdir1 = casedir1+indir
    inputdir2 = casedir2+indir
    inputdir3 = casedir3+indir
    inputdir4 = casedir4+indir
    print(inputdir1)
    outdir = "/Users/steven/OneDrive/Stanford/ONR project/results/BL_test/"
    
    tid1 = 168500
    fname1 = fname_fmt.format(probename,tid1)
    pos_name = pos_name_fmt.format(probename)
    
    
    tid2 = 109001
    fname2 = fname_fmt.format(probename,tid2)
    
    fname1 = os.path.join(inputdir1,fname1)
    pos_name1 = os.path.join(inputdir1,pos_name)
    print(fname1)
    print(pos_name1)
    
    fname2 = os.path.join(inputdir2,fname2)
    pos_name2 = os.path.join(inputdir2,pos_name)
    print(fname2)
    print(pos_name2)
    
    tid3 = 54501
    fname3 = fname_fmt.format(probename,tid3)
    fname3 = os.path.join(inputdir3,fname3)
    pos_name3 = os.path.join(inputdir3,pos_name)
    print(fname3)
    print(pos_name3)
    
    tid4 = 52501
    fname4 = fname_fmt.format(probename,tid4)
    fname4 = os.path.join(inputdir4,fname4)
    pos_name4 = os.path.join(inputdir4,pos_name)
    print(fname4)
    print(pos_name4)
    #parameters from data
        
        
    inds1,_,lines1,line_inds1 = readprobes(pos_name1)

    print(lines1)
    u_avg_x_lines1,u_rms_x_lines1,u_rms_y_lines1,u_rms_z_lines1 = readData(fname1,line_inds1,inds1)
    #plotQuantities(lines1,u_avg_x_lines1,"u_avg_x","smooth")
    plotQuantities(lines1,u_avg_x_lines1,"u_avg_x","smooth")
    plotQuantities(lines1,u_rms_x_lines1,"u_rms_x","smooth",prev_fig=1)
    inds2,_,lines2,line_inds2 = readprobes(pos_name2)
    u_avg_x_lines2,u_rms_x_lines2,u_rms_y_lines2,u_rms_z_lines2 = readData(fname2,line_inds2,inds2)
    #plotQuantities(lines2,u_avg_x_lines2,"u_avg_x","rough 0.025h")
    #plotQuantities(lines2,u_avg_x_lines2,"u_avg_x","rough 0.025h")
    #plotQuantities(lines2,u_rms_x_lines2,"u_rms_x","rough 0.025h")

    inds3,_,lines3,line_inds3 = readprobes(pos_name3)
    u_avg_x_lines3,u_rms_x_lines3,u_rms_y_lines3,u_rms_z_lines3 = readData(fname3,line_inds3,inds3)
    plotQuantities(lines3,u_avg_x_lines3,"u_avg_x","rough 0.0125h")
    plotQuantities(lines3,u_rms_x_lines3,"u_rms_x","rough 0.0125h",prev_fig=1)

    inds4,_,lines4,line_inds4 = readprobes(pos_name4)
    u_avg_x_lines4,u_rms_x_lines4,u_rms_y_lines4,u_rms_z_lines4 = readData(fname4,line_inds4,inds4)
    plotQuantities(lines4,u_avg_x_lines4,"u_avg_x","rough 0.01875h")
    plotQuantities(lines4,u_rms_x_lines4,"u_rms_x","rough 0.01875h",prev_fig=1)
if __name__ == "__main__":
    main()
