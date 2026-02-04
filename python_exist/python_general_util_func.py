import numpy as np
import matplotlib.pyplot as pyplot
import scipy

def computeDeq(h,L):
    '''
    computes equivalent diameter for the rectangular nozzle based on nozzle height and length
    save def. as Juen et al., 2022, AIAA Journal
    '''
    Deq = 1.3*(L*h)**(0.625)/(L+h)**(0.25)
    return Deq

def computeMj(NPR,gamma=1.4):
    """
    compute the Jet Mach number of the jet flow, this is defined as
    Mj = Uj/c_j
    Uj is the perfectly expanded jet velocity using isentropic rules
    c_j is the speed of jet 
    """
    Mj = np.sqrt(2/(gamma-1)*(NPR**((gamma-1)/gamma)-1)) #
    return Mj

def computeMa(NPR,gamma=1.4):
    """
    compute the acoustic Mach number of the jet flow
    this is based on the assumption that T_tot of the inflow is the same as the T_inf (TPR = 1, unheated jet)
    """
    Mj = computeMj(NPR,gamma)
    TjoT_inf = 1/(1+ (gamma-1)/2*Mj**2)
    cjoc_inf = np.sqrt(TjoT_inf)
    Ma = Mj*cjoc_inf
    return Ma


def calcSt_phys(freq,Uj,Deq):
    '''
    takes in physical frequency freq (Hz), physical jet velocity Uj (m/s), and Equivalent Diameter (m)
    return corresponding Strouhal number 
    '''
    St = freq * Deq/Uj
    return St

def calcSt_sim(f_ND,Ma,Deq_s):
    St = f_ND*Deq_s/Ma 
    return St

def calcSt2fND(St,Ma,Deq_s):
    '''
    takes in Strouhal number, Acoustic Mach number ( = Uj_ND), Deq_s equivalent diameter in nondimensional units
    return nondimension frequency f_ND corresponding to the St number
    '''
    f_ND = St * Ma/Deq_s
    return f_ND

