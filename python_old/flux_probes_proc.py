import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import python_general_util_func as gutil
from Universal_Subroutines import setPlotpref
setPlotpref()




#what needs to be done 
# a) read in the port files
# b) comprehend the result

#global variables
indir  = "/Users/steven/OneDrive/Stanford/ONR project/results/port_inject/131M"
outdir = "/Users/steven/OneDrive/Stanford/ONR project/results/port_inject/131M"

port_fmt = "/fluxprobes_{st}/{st}{ind}.fp"

def readPortprobes(fname):
    #INPUT fname: this is the full path to the file to be read (will be used directly)
    #OUTPUT portData: directionary containing mass_flux, x-momentum flux, y-momentum flux, z_momentum flux
    '''
    the format for the port probes as of May 23, 2025 is defined as follow:
    first two rows are of headers
    the columns as specified in the header are:
    1:name, 2:step, 3:time, 4:xp, 5:yp, 6:zp, 7:proj_area, 
    8:mass_flux, 9:mass_flux(comp(u,0)), 10:mass_flux(comp(u,1)) 11:mass_flux(comp(u,2)) 
    
    We will only be reading 8-11 in the above description
    '''
    # load raw data
    Data = np.loadtxt(fname, skiprows=2,usecols=(7,8,9,10))
    
    mf  = Data[:,0] #mass flux 
    xmf = Data[:,1] #x-momentum flux
    ymf = Data[:,2] #y-momentum flux
    zmf = Data[:,3] #z-momentum flux

    portData={
        "mf" : mf,
        "xmf" : xmf,
        "ymf" : ymf,
        "zmf" : zmf,
    }
    return portData
def readOutprobes(fname):
    #INPUT fname: this is the full path to the file to be read (will be used directly)
    #OUTPUT portData: directionary containing mass_flux, x-momentum flux, y-momentum flux, z_momentum flux
    '''
    the format for the outlet probes as of May 23, 2025 is defined as follow:
    first two rows are of headers
    the columns as specified in the header are:
    1:name, 2:step, 3:time, 4:xp, 5:yp, 6:zp, 7:proj_area, 
    8:mass_flux, 9:mass_flux(comp(u,0))
    
    We will only be reading 8-9 in the above description
    '''
    # load raw data
    Data = np.loadtxt(fname, skiprows=2,usecols=(7,8))
    
    mf  = Data[:,0] #mass flux 
    xmf = Data[:,1] #x-momentum flux
    

    outData={
        "mf" : mf,
        "xmf" : xmf,
    }
    return outData
    


def main():
    ## debug testing
    #port_bot_1= port_fmt.format(st="bot",ind = 1)
    #fname = indir + port_bot_1
    #print(fname)
    #probData = readPortprobes(fname)
    #print(probData["xmf"])
    
    #aggregated data
    mf_Ports_avg = 0 #this will be the combined mass flow of all ports averaged over samples
    mf_Ports_var = 0
    #setup reading of probes
    for i in range(6):
        i1 = i+1
        #get the file names
        fname_top = indir + port_fmt.format(st="top", ind=i1)
        fname_bot = indir + port_fmt.format(st="bot", ind=i1)
        #get the data
        portData_top = readPortprobes(fname_top)
        portData_bot = readPortprobes(fname_bot)

        #take mass flow rate and average
        mf_top_avg = np.average(portData_top["mf"]) 
        mf_top_var = np.var(portData_top["mf"]) 
        mf_bot_avg = np.average(portData_bot["mf"])
        mf_bot_var = np.var(portData_bot["mf"]) 
        print("mass flow top %.5f, port number %d" %(mf_top_avg,i1))
        print("mass flow bot %.5f, port number %d" %(mf_bot_avg,i1))
        mf_Ports_avg += (mf_top_avg + mf_bot_avg) #accumulate
        mf_Ports_var += (mf_top_var + mf_bot_var)
    print(mf_Ports_avg)

    #for outlet
    fname_out = indir + "/fluxprobes_out/outlet.fp"
    portData_out = readOutprobes(fname_out)
    mf_out_avg = np.average(portData_out["mf"])    
    mf_out_var = np.var(portData_out["mf"]) 
    print("mass flow out %.5f" %(mf_out_avg))

    MFR = mf_Ports_avg/mf_out_avg
    MFR_max = (mf_Ports_avg+np.sqrt(mf_Ports_var))/(mf_out_avg-np.sqrt(mf_out_var))
    MFR_min = (mf_Ports_avg-np.sqrt(mf_Ports_var))/(mf_out_avg+np.sqrt(mf_out_var))
    mesh_size = "131M"
    print("MFR using %s mesh is %.2f percent " %(mesh_size, MFR*100))
    print("MFR_max using %s mesh is %.2f percent " %(mesh_size, MFR_max*100))
    print("MFR_min using %s mesh is %.2f percent " %(mesh_size, MFR_min*100))

if __name__ == "__main__":
    main()