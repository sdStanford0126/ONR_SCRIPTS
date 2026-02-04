import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import python_general_util_func as gutil
from Universal_Subroutines import setPlotpref
setPlotpref()

outdir = "../"
def main():
    fname_fmt = "nearfield_probes_{:s}.txt"
    planes = ["xy", "xz"]
    #physical dimensions
    h = 0.475 #in    
    x_st = -0.305 #in
    x_sep = np.array([0.89,0.853,0.794,0.715,0.835,0.725,0.765]) #in
    offset_p = 3.07 #in
    x_pos_p = np.array([x_st])
    c_val = x_st
    for dx in x_sep:
        c_val += dx
        print(c_val)
        x_pos_p = np.append(x_pos_p,c_val)
    print("phsyical x positions" + str(x_pos_p))
    print("separation is %.2f" % (x_pos_p[-1] - x_st))
    print("number of mic is %d" % np.size(x_pos_p))

    #simulation dims
    x_pos = x_pos_p /h
    offset = offset_p/h
    N_mic = np.size(x_pos) # double (also do negative for better convergence)
    print("sim x pos" + str(x_pos))
    print("sim_offset %.4f" % offset)
    print(N_mic)
    header = 'x,y,z #near field probes'
    for plane in planes:
        fname =outdir + fname_fmt.format(plane)
        print(fname)
        xyz = np.zeros((N_mic*2,3))
        xyz[:,0] = np.append(x_pos,x_pos)
        if plane == "xy":
            #sampling the major plane (on xy plane, z = 0)
            xyz[:(N_mic),1] = offset
            xyz[(N_mic):,1] = -offset
        else:
            #sampling the minor plane (on xz plane, y = 0)
            xyz[:(N_mic),2] = offset
            xyz[(N_mic):,2] = -offset 
        np.savetxt(fname,xyz,delimiter=" ",header=header)
        
if __name__ == "__main__":
    main()
    