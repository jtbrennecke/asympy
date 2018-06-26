from pylab import *
from sklearn.cross_decomposition import PLSRegression
from scipy import stats
import mdtraj as md
from optparse import OptionParser

class pls:
    def __init__(self,Sobject,n_components,data):
        """The pls method implements the partial least squares based functional 
        mode analysis. 

        Parameters
        ----------
        Sobject:
            mdtraj trajectory. This is used to calculate the pls vectors.
        n_components:
            Number of dimensions used for the FMA.
            A test should be run using different dimensions to see which is 
            the most suitable for the dimension. In the suitable case the 
            cross-correlation value should be rather high.
        data:
            One dimensional input data. Has to match the number of data points
            (time) in trj_array.
        """
        try:
            trj_array = Sobject.xyz
        except:
            raise TypeError("Sobject has to be a mdtraj trajectory.")
        if not len(trj_array.shape) == 3:
            raise IOError("Shape of trj_array has to be 3d. See help.")
        if not type(data) == ndarray or not len(data.shape) == 1:
            raise TypeError("data has to be a 1d numpy array.")
        if not trj_array.shape[0] == data.shape[0]:
            raise IOError("First dimension of trj_array and data are time and have to match.")
        
        self._obj = PLSRegression(n_components=n_components,scale=False)
        
        # Create the object for full fitting.
        self._objfit = self._obj.fit(trj_array.reshape(trj_array.shape[0],-1), data)
        self.n_components=n_components
        
        # Create the object for ensemble weighted fitting.
        self._ew = PLSRegression(n_components=1,scale=False)
        self._ewfit = self._ew.fit(trj_array.reshape(trj_array.shape[0],-1), data)
        
    def proj(self, Sobject):
        """Projects the array onthe the calculated pls data
        
        Parameters
        ----------
        Sobject:
            mdtraj trajectory which is projected onto the pls vec.
        
        Returns
        -------
            numpy.array of shape (time, n_components) containing
            projection onto the pls space.
        """
        
        try:
            trj_array = Sobject.xyz
        except:
            raise TypeError("Sobject has to be a mdtraj trajectory.")
        if not len(trj_array.shape) == 3:
            raise IOError("Shape of trj_array has to be 3d. See help.")
        
        return self._obj.predict(trj_array.reshape(trj_array.shape[0],-1))
    
    @staticmethod
    def evaluate(Sobject, data, nframes=None, n_components_min=1, n_components_max=10, plot_b=False):    
        """Tries the pls  with components from n_components_min to n_components_max.
        The correlation coefficient for the data and the predicted projection is 
        calculated and returned as an array.
        
        To subsequently use the pls model n_components should be choosen to maximise
        the correlation between data and predicted projection in the validation set.
        
        Parameters
        ----------
        Sobject:
            traj of the mdtraj class. Frames in this model are used to 
            calculate the extremes and averages. These are used to interpolate.
        data:
            One dimensional input data. Has to match the number of data points
            in the Sobject.
        nframes:
            Number of frames in the Sobject used for building the pls model.
            By default half the frames of the trajectory are used and half 
            are used for validation.
        n_components_min:
            Minimum number of components used for model building.
        n_components_max:
            Maximum number of components used for model building.
        plot_b:
            If set to true a plot is generated showing the training data as 
            well as the validation data cross-correlation in comparision to 
            the number of components used.
        """        
        
        try:
            trj_array = Sobject.xyz
        except:
            raise TypeError("Sobject has to be a mdtraj trajectory.")
        if not len(trj_array.shape) == 3:
            raise IOError("Shape of trj_array has to be 3d. See help.")
            
        if not type(data) == ndarray or not len(data.shape) == 1:
            raise TypeError("data has to be a 1d numpy array.")
        if not trj_array.shape[0] == data.shape[0]:
            raise IOError("First dimension of trj_array and data are time and have to match.")

        if nframes==None:
            nframes=trj_array.shape[0]/2
        if n_components_min>n_components_max:
            raise ValueError("n_components_min has to be smaler than n_components_max.")

        # Create the array which will be filled with values and 
        # returned at the end.
        return_array = zeros((1+n_components_max-n_components_min,3))
        
        for i in range(n_components_min, n_components_max+1):
            print "Evaluating PLS component " + str(i) + " / " + str(n_components_max)
            testobj = pls(Sobject[:nframes], i, data[:nframes])    
            proj_data = testobj.proj(Sobject)
            #store values in array
            return_array[i-n_components_min,0] = i
            return_array[i-n_components_min,1] = corrcoef(proj_data[:nframes,0],data[:nframes])[0,1]
            return_array[i-n_components_min,2] = corrcoef(proj_data[nframes:,0],data[nframes:])[0,1]
        if plot_b:
            rc("figure", figsize=(8,6))
            title("Evaluate n_components", size=17, weight="bold")
            xlabel("n_components", size=14)
            ylabel("R", size=14)
            plot(return_array[:,0], return_array[:,1], label="training")
            plot(return_array[:,0], return_array[:,2], label="model")
            rcParams['legend.loc']='best'
            rcParams['legend.shadow'] = True
            rcParams['legend.fancybox'] = True
            rcParams['legend.frameon'] = True
            rcParams['legend.framealpha'] = 0.5
            rcParams['xtick.labelsize'] = 13
            rcParams['ytick.labelsize'] = 13
            legend(prop={'size':14})
        return return_array
    
    
    def extr(self, Sobject, nframes=10, filename=None):
        """Calculates an interpolation between the minimum and maximum along the
        pls vectors. This can be used to visualize the motion described by the pls.
        
        Parameters
        ----------
        Sobject:
            traj of the mdtraj class. Frames in this model are used to 
            calculate the extremes and averages. These are used to interpolate.
        nframes:
            Specifies the number of frames interpolating between the extremes.
        filename:
            If a filename is specified the interpolated trajectory is writen to the
            file.
            
        Returns
        -------
            mdtraj containing nframes interpolating between extremes in pls space.
        """
        
        try:
            Sobject.xyz
        except:
            raise TypeError("Sobject has to be a mdtraj trajectory.")
        
        _result = PLS._obj.coefs / sqrt( sum( PLS._obj.coefs**2 ))
        
        # Create the average structure of the trajectory
        M_avv = Sobject[0]
        M_avv.xyz = Sobject.xyz.mean(axis=0)
        
        # Remove the average structure 
        M_centered = Sobject[:Sobject.xyz.shape[0]]
        M_centered.xyz = Sobject.xyz - M_avv.xyz
        
        minval = min(dot(M_centered.xyz.reshape(Sobject.xyz.shape[0],-1),_result))
        maxval = max(dot(M_centered.xyz.reshape(Sobject.xyz.shape[0],-1),_result))
        diff = (maxval - minval)/(nframes - 1)
        
        fit_model = Sobject[:nframes]
        _result_reshape = _result.reshape(-1,3)
        for i in range(nframes):
            fit_model.xyz[i] = M_avv.xyz + (minval+diff*i) * _result_reshape
        
        if not filename==None:
            fit_model.save_pdb(filename)
        return fit_model 
    
    def ew_extr(self, Sobject, nframes=10, filename=None):
        """Calculates an interpolation between the minimum and maximum along the
        pls ensemble weighted vectors. This can be used to visualize the motion
        described by the pls.
        
        Parameters
        ----------
        Sobject:
            traj of the mdtraj class. Frames in this model are used to 
            calculate the extremes and averages. These are used to interpolate.
        nframes:
            Specifies the number of frames interpolating between the extremes.
        filename:
            If a filename is specified the interpolated trajectory is writen to the
            file.
            
        Returns
        -------
            mdtraj containing nframes interpolating between extremes in pls space.
        """
        
        try:
            Sobject.xyz
        except:
            raise TypeError("Sobject has to be a mdtraj trajectory.")
        
        _result = PLS._ew.coefs / sqrt( sum( PLS._ew.coefs**2 ))
        
        # Create the average structure of the trajectory
        M_avv = Sobject[0]
        M_avv.xyz = Sobject.xyz.mean(axis=0)
        
        # Remove the average structure 
        M_centered = Sobject[:Sobject.xyz.shape[0]]
        M_centered.xyz = Sobject.xyz - M_avv.xyz
        
        minval = min(dot(M_centered.xyz.reshape(Sobject.xyz.shape[0],-1),_result))
        maxval = max(dot(M_centered.xyz.reshape(Sobject.xyz.shape[0],-1),_result))
        diff = (maxval - minval)/(nframes - 1)
        
        fit_model = Sobject[:nframes]
        _result_reshape = _result.reshape(-1,3)
        for i in range(nframes):
            fit_model.xyz[i] = M_avv.xyz + (minval+diff*i) * _result_reshape
        
        if not filename==None:
            fit_model.save_pdb(filename)
        return fit_model 
    
    
    
    def _partial_predict(self, trj_array, Cn):
        """Predict only the parts of the total motion. 
        For this the projection of the PLS is split into parts.
        
        Parameters
        ----------
        trj_array:
            Numpy array consisting of a time dimension and two dimension of the 
            positions (xyz).
            
        Cn:
            (int) Number of chains which are used for symmetry calculations.
            
        Returns
        -------
        Prediction, [Partial Predictions]:
            (tuple) first element, the overall prediction is returned.
            As the second element, the partial predictions are returned for
            the subparts (Cn long list).
        """
        
        if not type(trj_array) == ndarray:
            raise TypeError("trj_array has to be a numpy array encoding the trajectory.")
        if not len(trj_array.shape) == 3:
            raise IOError("Shape of trj_array has to be 3d. See help.")        

        oo = self._obj
        X = copy(trj_array.reshape(trj_array.shape[0],-1))
        X -= oo.x_mean_
        X /= oo.x_std_

        Ypred = dot(X[:,:],oo.coefs[:])
        Ypred += oo.y_mean_

        Ypred_list = []
        for i in range(Cn):
            Ypred_list.append( dot(X[:,X.shape[1]*i/Cn:X.shape[1]*(i+1)/Cn],oo.coefs[oo.coefs.shape[0]*i/Cn:oo.coefs.shape[0]*(i+1)/Cn]) )
            Ypred_list[-1] += oo.y_mean_/Cn

        return Ypred, Ypred_list
    
    def _pred_contribution(self,tot,su1,vmin=0.,vmax=1., epsilon=1e-5):
        """Auxiliary function to predict the contributions"""
        if abs(tot * vmax) < abs(su1):
            return nan
        elif abs(tot * vmin) > abs(su1):
            return nan
        else:
            counter = 0
            while (vmax - vmin) > epsilon and counter<1e7:
                counter += 1
                tmp = (vmax + vmin) / 2.
                if abs(tot * tmp) < abs(su1):
                    vmin = tmp
                else:
                    vmax = tmp
            return (vmax + vmin) / 2.
        
    def _contributions(self,tot,su1,vmin=0.,vmax=1., epsilon=1e-5):
        """Auxiliary function to predict the contributions"""
        if type(tot) != np.ndarray or type(su1) != np.ndarray:
            print type(tot)
            print type(su1)
            raise TypeError("Inputs must be numpy arrays.")
        elif tot.shape != su1.shape:
            raise IndexError("Shapes of inputs do not match.")
        rotations = []
        for to,su in zip(tot.flatten(),su1.flatten()):
            rotations.append(self._pred_contribution(to,su,vmin,vmax,epsilon))
        return array(rotations)
    
    def contributions(self,Sobject, Cn,epsilon=1e-5):
        """Predict the contribution of each of the subunits to the overall symmetry.
        

        Parameters
        ----------
        Sobject:
            traj of the mdtraj class. Used to predict the contributions.
            
        Cn:
            (int) Number of chains which are used for symmetry calculations.
            
            
        Returns
        -------
        [contributions]:
            list of the contributions of each subunit at any time.
        """        
            
        try:
            trj_array = Sobject.xyz
        except:
            raise TypeError("Sobject has to be a mdtraj trajectory.")
        
        # Do the partial predictions first.
        # The result is needed to predict the contributions of each of
        # subunits to the overall motion.
        Ypred, Ypred_list = self._partial_predict(trj_array, Cn)
        
        contributions = []
        for Ypred_el in Ypred_list:
            contributions.append(self._contributions(Ypred,Ypred_el,epsilon=epsilon))
        
        return contributions,Ypred, Ypred_list    

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f","--file", dest="filename", help="xtc file to process")
    parser.add_option("-s","--structure", dest="structure", help="pdb file")
    parser.add_option("-d","--data", dest="data", help="data file 1 or 2d input")
    parser.add_option("-n","--subunits", dest="subunits", help="number of subunits")
    parser.add_option("-c","--components", dest="components", default=-10, help="Number of PLS components, if negative will evaluate number of components to be used.")
    parser.add_option("-o","--out", dest="outfile",default="comp.txt", help="file of the asymmetric contributions")
    parser.add_option("-w","--prediction", dest="plsfma",default="model.dat", help="predictions for the functional input")
    parser.add_option("-e","--ew", dest="ensemble",default="extremes_ew.pdb", help="pdb output of the ensemble weighted pls")
    parser.add_option("-p","--pls", dest="plsout",default="extremes.pdb", help="pdb output of the pls motion")
    (options,args) = parser.parse_args()
    

    
    try:
        comp = int(options.components)
    except:
        raise IOError("-c option has to be an integer")
       
    try:
        dat = loadtxt(options.data)
    except:
        raise IOError("-d file could not be loaded")

        
    traj = md.load(options.filename, top=options.structure)

       
    if comp < 1:
        
        if len(dat.shape) == 1:
            print "evaluating components for pls"
            res = pls.evaluate(traj, dat, n_components_max=abs(comp) )
        elif len(dat.shape) == 2:
            print "evaluating components for pls"
            res = pls.evaluate(traj, dat[:,1], n_components_max=abs(comp) )
        else:
            print "data file has to be 1 or 2d input data."
            exit(1)
            
        print "writing evaluation data to " + options.outfile       
        savetxt(options.outfile , res,header="# Nr. Components, training, cross-val")
        
    else:       
        try:
            sub = int(options.subunits)
        except:
            raise IOError("-n option has to be an integer")
    
        if len(dat.shape) == 1:
            print "Calculating pls model"
            PLS = pls(traj,comp, dat)
        elif len(dat.shape) == 2:
            print "Calculating pls model"
            PLS = pls(traj,comp, dat[:,1])
        else:
            print "data file has to be 1 or 2d input data."
            exit(1)
            
        contr,PLSprediction,PLSpredictionPartial = PLS.contributions(traj,sub)
        print "writing contributions to " + options.outfile
        savetxt(options.outfile, array(contr).T, header="# Contributions of subunits")
        print "writing prediction to " + options.plsfma
        savetxt(options.plsfma, concatenate((array([PLSprediction]),array(PLSpredictionPartial)),axis=0)[:,:,0].T, header="# prediction | prediction SU1 | ...")
        print "writing extremes to " + options.plsout
        PLS.extr(traj, filename=options.plsout)
        print "writing extremes_ew to " + options.ensemble
        PLS.ew_extr(traj, filename=options.ensemble)
        print "finished successfully"
    
    
    
    
    
    
    
    
    
    
