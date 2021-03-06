from pylab import *
import mdtraj as md
from scipy.optimize import fsolve
from optparse import OptionParser

try:
    from ipy_progressbar import ProgressBar
    pbb = True
except:
    pbb = False

def RMatrix_Axis_Angle(axis, angle, normed=False, deg=False, ):
    """ 
    Returns the rotation matrix to rotate around a axis (vector) centered at
    the origin by an angle 
    Parameters
    ----------
    axis : ndarray (T,3)
        vector around wich to rotate for each time T
    angle : ndarray (T)
        angle in rad to rotate by for each time T
    normed : Bool (False)
        if True will expect a normalized axis (of length 1 in cartesian space)
    deg : Bool (False)
        if True indicates that the angle is given in degrees instead of radians

    Returns
    -------
    narray ( T,3,3 )
    Rotation matrix for every time T
    """
    ### input parsing
    if normed:
        ax = axis
    else:
        ax = axis / linalg.norm(axis, axis=1)[:, None]

    if deg:
        ang = deg2rad(angle)
    else:
        ang = angle

    ### calculations 
    cp = cos(ang)
    a_sp = (sin(ang)[:, None] * ax) 
    RMat = einsum("ta,tb->tab", ax, ax) * (1 - cp)[:, None, None]
    #diagonal elements
    RMat[:,0,0] += cp
    RMat[:,1,1] += cp
    RMat[:,2,2] += cp
    #off diagonal elements
    RMat[:,1,2] += a_sp[:,0]
    RMat[:,2,1] -= a_sp[:,0]
    RMat[:,2,0] += a_sp[:,1]
    RMat[:,0,2] -= a_sp[:,1]
    RMat[:,0,1] += a_sp[:,2]
    RMat[:,1,0] -= a_sp[:,2]

    return RMat

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

def calc_rot_axis(B,Cn):
    """This function calculates the optimal rotation axis following
    Pinsky et al. 2008.
    
    Parameters
    ----------
    B:  structure to calculate the rotation on
    Cn: Symmetry
    
    Returns
    -------
    rot_axis: 3D vector of the optimal rotation axis
    
    """
    # Calculate the A Matrix following Pinsky et al. 2008
    Am = zeros((3,3))
    for rot in range(1,Cn):
        prefactor = (1 - cos( rot * 2 * pi / Cn ))
    if prefactor == 0:
        return None
    Atmp = zeros((3,3))
    for k in range(B.n_atoms):
        Atmp += outer(B.xyz[0][k],B.xyz[rot][k])
        Atmp += outer(B.xyz[rot][k],B.xyz[0][k])
    Am += prefactor * Atmp
    
    
    # Calculate the eigenvector and eigenvalues of A
    eig_val,eig_vec = linalg.eigh(Am)
    eig_max = argmax(eig_val)
    
    
    # If we have a Cn2 symmetry we are done.

    if Cn == 2:
        # The eigenvector corresponding to the largest eigenvalue
        # gives the rotation axis.
        rot_axis = eig_vec[eig_max]
        
    # If we have more subunits than two we need some more 
    # to calculate the rotation axis.        
    else: 
        # Calculate B Matrix (Pinsky et al. 2008)
        Bm = zeros(3)
        for rot in range(Cn):
            prefactor = sin( rot * 2 * pi / Cn )
            if prefactor == 0:
                continue
            Btmp = zeros(3)
            for k in range(B.n_atoms):
                Btmp += cross( B.xyz[0][k] , B.xyz[rot][k] )
            Bm += prefactor * Btmp

        # Get the lambda_max
        f = lambda lambda_max: sum([(dot(eig_vec[i],Bm)/(eig_val[i]-lambda_max))**2 for i in range(3)])-1
        # by solving the f function.
        # The solution is not unique. 
        # However, by starting by from the largest eigenvalue we get the solution we are looking for.
        # If we start exactly at the eigenvalue we run into a singularity (hence eig_val + 1).
        lambda_max = fsolve(f, eig_val[eig_max] +1)
        m_max = sum([eig_vec[i] * (dot(eig_vec[i],Bm)/(eig_val[i]-lambda_max)) for i in range(3)], axis=0)
        # And finally we get the optimal rotation axis m_max
        rot_axis = m_max
       
    # END ELSE
    
    return rot_axis
    
def create_permutations(Cn):
    """Calculate all possible permutations of the chains.
    
    Parameters
    ----------
    Cn:  Number of residues
    
    Returns
    -------
    perm_list: list of residues which have to be permutated.
    """
    perm_list = []
    for i in range(2**(Cn-1)-1):
        perm_list.append(list(str(int("".join(list(bin(i+1))[2:]))).zfill(Cn)))
    return perm_list

def CSM(Sobject, Cn,chains = None, Rotfunc=Rot,perm_opt = False):
    """Calculates the continuous symmetry measure as
    defined by Zabrodsky et al. (Continuous Symmetry
    Measures JAmChemSoc 1992, 114, 7843-7851) for the
    calculation of the optimal rotation axes the 
    analytical solution was used from with Pinsky et al. 
    (Analytical Methods for Calculating Continuous
    Symmetry Measures and the Chirality 
    JComputChem 29: 2712-2721, 2008)
    
    Parameters
    ----------
    Sobject:
        Model or SelObj of the ngmx class. In the 
        current implementation this has to be a single 
        frame!
    Cn: (int)
        Cn is the symmetry group. Currently only Cn 
        symmetry is supported.
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

    Returns
    -------
    asymmetry_measure: (float)
        A value between 0. and 1. where 1. means perfect 
        symmetry according to Cn and 0. no symmetry.
    symm_struct: (SelObj)
        Symmetric structure corresponding to trajectory
    """
    
    Smean = Sobject[:]
    # Create the structure which will be used
    # to calculate the symmetric structure.
    B = Sobject[:Cn]
        
    Sn = []
    
    if pbb:
        # Use progessbar to show the progress
        pb = ProgressBar(Sobject.n_frames, title='Progress')
        pb.start()
        
    atNCn = Sobject.n_atoms / Cn
    if chains == None:
        chains = []
        for i in range(Cn):
            chains.append([atNCn * i,atNCn * (i+1)]) 
    # Do the optimization if activated:
    if perm_opt:       
        # dict of all ambiguous residues with the matching residues
        ambiguous_atoms = {"VAL":[["CG1","CG2"],["CG2","CG1"]],
                          "TYR": [["CD1","CD2","CE1","CE2"],["CD2","CD1","CE2","CE1"]],
                          "ASP": [["OD1","OD2"],["OD2","OD1"]],
                          "GLU": [["OE1","OE2"],["OE2","OE1"]]}

        table, bonds = Sobject.topology.to_dataframe()
        # find all residues which are of the type defined in keys
        keys = ambiguous_atoms.keys()
        ambiguous_residues = []
        for key in keys:
                ambiguous_residues.append(list(unique(table[table.resName == key].resSeq)))
        ambiguous_residues = [item for sublist in ambiguous_residues for item in sublist]


        # Create the lists which atoms have to be exchanged from the original sorting 
        # for any of the permutations.
        from_id_list = []
        to_id_list = []
        indices_list = []

        perms = create_permutations(Cn)
        for perm in perms:
            indices_list.append([i for i, x in enumerate(perm) if x == '1'])
        # Get the list of chains (indices) which chains have to be relabeled:
        for indices in indices_list:
            from_id = []
            to_id = []

            # And find the atoms to be relabeled for any of the chains:
            for chain in indices:
                for sel_res in ambiguous_residues:

                    #First identify the residues:
                    dfs = table[(table.resSeq == sel_res) & (table.chainID == chain)]
                    # find the id of the atoms we need to permute from
                    for at in ambiguous_atoms[unique(dfs.resName)[0]][0]:
                        from_id.append( dfs[dfs.name == at].index.tolist()[0] )

                    # and to permute to
                    for at in ambiguous_atoms[unique(dfs.resName)[0]][1]:
                        to_id.append( dfs[dfs.name == at].index.tolist()[0] )
                        
            from_id_list.append(from_id)
            to_id_list.append(to_id)
    

    # END OPTIMIZATION
    for t in range(Sobject.n_frames):

        # Center COG at (0,0,0)
        Sobject.xyz[t] = Sobject[t].xyz - Sobject[t].xyz.mean(axis=1)

        # rotate to have all possible rotations
        # and relabel them
        for rot in range(Cn):
            B.xyz[rot] = Sobject.xyz[t]
            B.xyz[rot] = Rotfunc(rot,B.xyz[rot],chains)
            

        rot_axis = calc_rot_axis(B,Cn)
        
        if rot_axis is None:
            continue
        
        # We rotate around the axis and average over the rotations
        # to get the closest symmetric structure.
        for rot in range(1,Cn):
            # Calculate rotation matrix around rot_axis by an angle of rot*360/Cn
            RMat = RMatrix_Axis_Angle([rot_axis],[ rot * 360./Cn ],deg=True)
            # Do the rotation using this matrix
            B.xyz[rot] = einsum("tnc,tcp->tnp",B[rot].xyz,RMat,casting='same_kind')
        # And average over it
        Smean.xyz[t] = B.xyz.mean(axis=0)
        
        
        # optimize if activated:
        if perm_opt:
            rmsd_list = []
            B_opt = B[0]

            for i, sel_res in enumerate(ambiguous_residues):
                indlist = table[(table.resSeq == sel_res)].index.tolist()
                rmsd_list.append(sum((Smean.atom_slice(indlist).xyz[t] - B.atom_slice(indlist).xyz[0])**2))
                
            for from_id,to_id in zip(from_id_list,to_id_list):
                # Do the permutation
                chains_to_perm = []
                
                # reset B
                for rot in range(Cn):
                    B.xyz[rot] = Sobject.xyz[t]
                    B.xyz[rot][to_id] = B[rot].xyz[0][from_id]
                    # And do the rotation again
                    B.xyz[rot] = Rotfunc(rot,B.xyz[rot],chains)
                    
                for rot in range(1,Cn):
                    # Calculate rotation matrix around rot_axis by an angle of rot*360/Cn
                    RMat = RMatrix_Axis_Angle([rot_axis],[ rot * 360./Cn ],deg=True)
                    # Do the rotation using this matrix
                    B.xyz[rot] = einsum("tnc,tcp->tnp",B[rot].xyz,RMat,casting='same_kind')
                    
                # Calculate the new Mean
                Smean.xyz[t] = B.xyz.mean(axis=0)
                
                # Calculate the distance for all relevant residues
                for i, sel_res in enumerate(ambiguous_residues):
                    indlist = table[(table.resSeq == sel_res)].index.tolist()
                    rmsd = sum((Smean.xyz[t][indlist] - B.xyz[0][indlist])**2)
                
                    # Compare new distance value to old one
                    # and than replace the residue in Smean
                    # and in the optimized B
                    if rmsd_list[i] > rmsd:
                        rmsd_list[i] = rmsd
                        B_opt.xyz[0][indlist] = B[0].xyz[0][indlist]

                        
            # calculate the new optimal rotation axis
            for rot in range(Cn):
                B.xyz[rot] = B_opt.xyz[0]
                B.xyz[rot] = Rotfunc(rot,B.xyz[rot],chains)


            rot_axis = calc_rot_axis(B,Cn)

            #if rot_axis is None:
            #    continue
                
            # calculate the S value from the optimal value
            # We rotate around the axis and average over the rotations
            # to get the closest symmetric structure.
            for rot in range(1,Cn):
                # Calculate rotation matrix around rot_axis by an angle of rot*360/Cn
                RMat = RMatrix_Axis_Angle([rot_axis],[ rot * 360./Cn ],deg=True)
                # Do the rotation using this matrix
                B.xyz[rot] = einsum("tnc,tcp->tnp",B[rot].xyz,RMat,casting='same_kind')
            # And average over it
            Smean.xyz[t] = B.xyz.mean(axis=0)
                    
        
        # END OPT
        
        # d_sq: Square of root mean square size of the object.
        #       Used to calculate the correct normalization of the CSM.
        d_sq = sum(B.xyz[0] ** 2)
        
        # CSM: Is the distance of the original object to the closest
        #      symmetric object (Smean) with a normalization.
        #      
        S = [t, 100./d_sq * sum((B.xyz[0] - Smean.xyz[t]) ** 2)]
        
        # The CSM can also be calculated based on the individual subunits.
        # By this the individual contributions to the asymmetry can be detected.
        for ch in chains:
            S.append(100./d_sq * sum((B.xyz[0,ch[0]:ch[1]] - Smean.xyz[t,ch[0]:ch[1]]) ** 2))
        Sn.append(S)
            


        # If ProgressBar works: advance it.
        if pbb:
            pb.advance()
    if pbb:
        pb.finish()
    return array(Sn),Smean
    
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f","--file", dest="filename", help="xtc file to process")
    parser.add_option("-s","--structure", dest="structure", help="pdb file")
    parser.add_option("-n","--subunits", dest="subunits", help="number of subunits")
    parser.add_option("-o","--outtxt", dest="outtxt", default="asym.txt", help="output file of asymmetry measure")
    parser.add_option("-e","--outxtc", dest="outxtc", default="asym.xtc", help="output file for symmetric trajectory")
    parser.add_option("-c","--chains", dest="chains", default=None, help="list of chains if sorting is non standard from atom to atom usage: -c '[[0,10],[20,30],[10,20]]'")
    parser.add_option("-p","--permutations", dest="permutations", action="store_true", default=False, help="Use permeation optimization of ambiguous residues.")
    (options,args) = parser.parse_args()
 
    try:
        sub = int(options.subunits)
    except:
        raise IOError("-s option has to be an integer")
 
    print "loading trajectory ..."
    traj = md.load(options.filename, top=options.structure)
    
    print "calculating symmetry measure ..."    
    if options.chains == None:
        csm = CSM(traj, sub,perm_opt=options.permutations)
    else:
        try:
            chains = eval(options.chains)
        except:
            raise IOError("Could not convert chains to list")
        csm = CSM(traj, sub, eval(options.chains), perm_opt=options.permutations)

    print "writing CSM measure to " + options.outtxt        
    savetxt(options.outtxt, csm[0],header="# CSM overall, CSM for subunits")
    print "writing symmetric coorindates to " + options.outxtc
    csm[1].save_xtc(options.outxtc)
    
    
    
    
    
    
    
    
    
