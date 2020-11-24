
import pyfftw
from math import *
import numpy as np
from scipy.special import erfinv
from scipy import misc
from scipy.signal import convolve2d
# from scipy.special import kv as Kv
# from itertools import product
import os, sys, csv
from time import time
from multiprocessing import Pool, Process
# from joblib import Parallel, delayed
from tqdm import tqdm
import matplotlib.pyplot as plt
# from pyevtk.hl import imageToVTK 

from RandomFieldModule.GaussianRandomField import GaussianRandomField
from RandomFieldModule.utilities.common import *
from RandomFieldModule.utilities.ErrorMessages import *
from RandomFieldModule.utilities.Exports import exportVTK




#######################################################################################################

def levelcut(field, level=0):
    phase = np.where(field > level, 1, 0)
    return phase.astype(np.intc)

def get_vf(phase):
    return np.mean(phase)


#######################################################################################################
#	Random material generator class
#######################################################################################################

class RandomMaterial:

    def __init__(self, grid_level, ndim=2, verbose=0, levelcut_strategy="abs", **kwargs):

        self.verbose = verbose

        self.ndim = int(ndim) # dimension 2D or 3D
        if not (self.ndim==2 or self.ndim==3): msgDimError(self.ndim)

        self.N  = int(2**grid_level)
        self.Nd = self.N * np.ones(self.ndim, dtype=np.intc)
        self.nvoxels = np.prod(self.Nd)

        ### Level-cut strategy:
        ###   'abs' - level-cut of the abs(field)
        ###   'sym' - standard level-cut
        self.levelcut_strategy = levelcut_strategy

        vf = kwargs['vf']
        self.set_level(vf)

        ### Gaussian Random Field (Intensity)
        if self.verbose: print('\nBuilding GRF...\n')
        self.GRF = GaussianRandomField(grid_level=grid_level, ndim=ndim, verbose=verbose, **kwargs)




    #--------------------------------------------------------------------------
    #   Updates
    #--------------------------------------------------------------------------
    
    def set_level(self, vf=None, tau=None):
        if vf is not None:
            self.vf  = vf
            self.tau = vf2tau(vf, strategy=self.levelcut_strategy)
        elif tau is not None:
            self.tau = tau
            self.vf  = tau2vf(tau, strategy=self.levelcut_strategy)
        else:
            print('Either expected volume fraction or level set value must be given.')
    
    def reseed(self, seed=None):
        self.GRF.reseed(seed)

    #--------------------------------------------------------------------------
    #   Sampling
    #--------------------------------------------------------------------------
    
    ### Generate a realization

    def sample(self, noise=None):
        field = self.GRF.sample(noise=noise)
        if self.levelcut_strategy is 'abs': field = np.abs(field)
        field = levelcut(field, self.tau)
        return field


    ### Generate a family of realizations
    def generate_samples(self, nsamples=1, path=None, output_format="png", append=False):
        output = False if path is None else True

        if output:
            if not append or not hasattr(self, 'sample_count'):
                os.system('rm -f ' + path + 'sample_*')
                self.sample_count = 0

        time_start = time()

        expected_vf = 0
        for isample in tqdm(range(nsamples)):
            phase = self.sample()
            expected_vf += phase.mean()
            if output:
                self.sample_count += 1
                filename = path + 'sample_{0:d}'.format(self.sample_count)
                if   self.ndim==2 and output_format is "png": self.save_png(phase, filename)
                elif self.ndim==3 or  output_format is "vtk": self.save_vtk(phase, filename)
        expected_vf /= nsamples

        print('All samples generation time: {0} s'.format(time()-time_start))
        print('Volume fraction: {0}'.format(expected_vf))


    #--------------------------------------------------------------------------
    #   EXPORTS
    #--------------------------------------------------------------------------

    def save_png(self, phase, filename):
        misc.imsave(filename, (1-phase))

    def save_vtk(self, phase, filename):
        exportVTK(filename, cellData = {'phase' : phase})

    #--------------------------------------------------------------------------
    #   TESTS
    #--------------------------------------------------------------------------

    def test_Covariance(self, nsamples=1000):        
        return self.GRF.test_Covariance(nsamples)


    def test_VolumeFraction(self, nsamples=1000):
        from RandomMaterial.utilities import compute_ProbaDist, fit_ProbaDist, vf2tau, tau2vf

        data = []
        for isample in tqdm(range(nsamples)):
            phase = self.sample()
            vf = phase.mean()
            data.append(vf)

        p, x = compute_ProbaDist(data)
        m, sigma, p_fit = fit_ProbaDist(x, p, type='LogNormal')
        print('m={0}, s={1}'.format(exp(m), sigma))
        print('expected vf = ', np.mean(data))

        plt.plot(x,p)
        plt.plot(x,p_fit)
        plt.axvline(x=self.vf, color='black')
        plt.legend(['vf pdf', 'Fit'])
        plt.show()
        return 0


    def test_TwoPointProbability(self, nsamples=1000):
        from RandomMaterial.utilities import autocorrelation
        from RandomMaterial.reconstruction import RadialAverage
        from utilities.image_statistics import LevelSetGaussian as LSG

        C = np.zeros(self.Nd)
        for isample in tqdm(range(nsamples)):
            phase = self.sample()
            C += autocorrelation(phase)
        C /= nsamples

        S2, Var = RadialAverage(C)

        h = 1/self.N
        Slope = (S2[1]-S2[0])/h

        nu  = self.GRF.Covariance.nu
        rho = self.GRF.Covariance.corr_len[0]
        Slope0 = LSG.S2_slope_at_zero_Matern(self.tau, nu, rho)

        print("nu = ", nu)
        print("Slope = ", Slope)
        print("Slope (analytic) = ", Slope0)

        return 0

#######################################################################################################


#######################################################################################################
#                                        Run as main (for testing)
#######################################################################################################
if __name__ == "__main__":

    import importlib
    config = importlib.import_module(sys.argv[1])

    RM = RandomMaterial(config)


        
