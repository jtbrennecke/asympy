from pylab import *
import mdtraj as md
from optparse import OptionParser

def Rot(index,B,chains):
    """Change the sorting of chains. It is changed
    by as many chains as index tells.
    Performes an in place rotation.

    Parameters
    ----------
    index:  (int)  The number of chains it is rotated by.
    B:      (array) array of the positions of trajectory.
    chains: (list of lists) a list containing a list of 
            the first and last atom of any chain.
            
    Returns
    -------
    array of the rotated positions.
    
    """
    coordinates = array([B[ch[0]:ch[1]] for ch in chains],copy=True)
    for i,ch in enumerate(chains):
        B[ch[0]:ch[1]] = coordinates[(i+index)%len(chains)]
    return B

def prepare_files(traj,dat,Cn,chains = None, Rotfunc=Rot):
    """Prepares the input files for pls. It duplicates the input
    trajectory as well as the data input. The input xtc will
    be appended in the rotated versions to ensure symmetry
    of the pls vector.
    
    Parameters
    ----------
    traj:
        mdtraj trajectory file
    dat:
        1d or 2d numpy array
    Cn:
        symmetry group which will later be used for pls
    chains:
        Can specify the chains of a protein by giving
        a list of start and end atomindices of the 
        chains. This is only usefull if the sorting 
        of the protein is not in a single direction 
        (clockwise or anti-clockwise but is sorted 
        differently). For more complex sortings see 
        Rotfunc option.
    Rotfunc:
        If for the rotation it is not enough to cycle
        the protein but to have a more complex symmetry
        a user defined rotation function can be 
        implemented.
    """

    traj_tot = traj[:]
    dat_tot = dat
    
    atNCn = traj.n_atoms / Cn
    if chains == None:
        chains = []
        for i in range(Cn):
            chains.append([atNCn * i,atNCn * (i+1)])
    
    for n in range(1,Cn):
        B = traj[:]
        for t in range(traj.n_frames):
            B.xyz[t] = Rotfunc(n,B.xyz[t],chains)
        traj_tot = traj_tot.join(B)
        dat_tot = append(dat_tot,dat,axis=0)
        
    traj_tot.superpose(traj_tot,parallel=False)
    traj_tot.superpose(traj_tot,parallel=False)
    return traj_tot, dat_tot
    
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f","--file", dest="filename", help="xtc file to process")
    parser.add_option("-s","--structure", dest="structure", help="pdb file")
    parser.add_option("-d","--data", dest="data", help="data file 1 or 2d input")
    parser.add_option("-n","--subunits", dest="subunits", help="number of subunits")
    parser.add_option("-x","--xtcout", dest="xtcout", default="xtcout.xtc", help="xtc output of correctly rotated trajectory")
    parser.add_option("-p","--datout", dest="datout", default="datout.xtc", help="dat correct data file matching xtcout")
    parser.add_option("-c","--chains", dest="chains", default=None, help="list of chains if sorting is non standard usage: -c '[[0,10],[20,30],[10,20]]'")
    (options,args) = parser.parse_args()
    
    print "loading trajectory"
    traj = md.load(options.filename, top=options.structure)

    try:
        sub = int(options.subunits)
    except:
        raise IOError("-n option has to be an integer")
    
    print "loading data file"
    try:
        dat = loadtxt(options.data)
    except:
        raise IOError("-d file could not be loaded")
    
    try:
        chains = eval(options.chains)
    except:
        raise IOError("Could not convert chains to list")
    
    print "processing files..."
    res = prepare_files(traj,dat,sub,chains = None, Rotfunc=Rot)
    
    print "saving xtc to " + options.xtcout
    res[0].save_xtc(options.xtcout)
    
    print "saving data file to " + options.datout
    savetxt(options.datout,res[1])
    
    print "finished successfully"
    
    
    
    
    
