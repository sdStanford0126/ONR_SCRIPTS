import numpy as np
import scipy 
import matplotlib.pyplot as pyplot
import python_general_util_func as gutil



def main():
    #jet physical configuration

    AR = 2.0
    NPR = 3.0
    h = 0.475 / 39.37           # [m]
    L = AR*h                     # [m]
    Deq = gutil.computeDeq(h,L) # [m]
    T_inf = 297.3960            # [K]
    R = 287.0528                # [J/kg*K]
    gamma = 1.4
    c_inf = np.sqrt(gamma*R*T_inf)
    Ma = gutil.computeMa(NPR,gamma)
    Uj = Ma*c_inf
    
    #query limits
    fAcq = 204800
    minfact = 8.0 #this is the factor to modify the minimum requirements based on realization count
    fmin = 300*minfact
    fmax = fAcq/2
    St_min = gutil.calcSt_phys(fmin,Uj,Deq)
    St_max = gutil.calcSt_phys(fmax,Uj,Deq)
    print(f"St_min is {St_min}\n")
    print(f"St_max is {St_max}\n")
    #simulation frequency limits
    h_s = 1
    L_s = h_s*AR
    Deq_s = gutil.computeDeq(h_s, L_s)
    f_ND_min = gutil.calcSt2fND(St_min,Ma,Deq_s)
    f_ND_max = gutil.calcSt2fND(St_max, Ma, Deq_s)

    print(f"f_ND_min is {f_ND_min} \n")
    print(f"f_ND_max is {f_ND_max} \n")

        
    f_samp = f_ND_max * 3
    print(f"T_samp is {1/f_samp}")
    St_samp = gutil.calcSt_sim(f_samp,Ma,Deq_s)
    print(f"St_samp: %.2f \n" % St_samp)
    dt = 0.001
    f_fwh = 1/(50*dt)
    St_fwh = gutil.calcSt_sim(f_fwh,Ma,Deq_s)
    print(f"St_fwh: %.2f \n" % St_fwh)
    step_samp = (1/f_samp)/dt
    T_samp = 1/f_ND_min * 10

    print(f"step_samp is {step_samp}\n")
    print(f"T_samp is {T_samp}\n")

    St_min = gutil.calcSt_sim(1/T_samp,Ma,Deq_s) 
    print(f"St_min is : {St_min}\n")

    f_med = gutil.calcSt2fND(0.2,Ma,Deq_s)
    N_f_med = T_samp//(1/f_med) 
    print(f"number of realization of St0.2 given smapling time {N_f_med}")
if __name__ == "__main__":
    main()
