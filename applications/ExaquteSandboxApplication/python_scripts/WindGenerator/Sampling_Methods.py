
from math import *
import numpy as np
import scipy.fftpack as fft
from time import time
# from WindGeneration.utilities.common import FourierOfGaussian, SpacialCovariance, autocorrelation
# 



METHOD_DST  = 'dst'
METHOD_DCT  = 'dct'
METHOD_FFT  = 'fft'
METHOD_FFTW = 'fftw'
METHOD_VF_FFTW = 'vf_fftw'
METHOD_NFFT = 'nfft'
METHOD_H2   = 'H2'
METHOD_H2_hlibpro = 'H2_hlibpro'
METHOD_H2_h2tools = 'H2_h2tools'
METHOD_ODE = 'ODE'
METHOD_RAT = 'Rational'
METHOD_VF_RAT_HALFSPACE_RAPID_DISTORTION = 'vf_rat_halfspace_rapid_distortion'


#######################################################################################################

class Sampling_method_base:

    def __init__(self, RandomField):
        self.verbose = RandomField.verbose
        self.L, self.Nd, self.ndim = RandomField.L, RandomField.ext_grid_shape, RandomField.ndim
        self.DomainSlice = RandomField.DomainSlice


class Sampling_method_freq(Sampling_method_base):

    def __init__(self, RandomField):
        super().__init__(RandomField)
        L, Nd, d = self.L, self.Nd, self.ndim
        self.Frequences = [(2*pi/L[j])*(Nd[j]*fft.fftfreq(Nd[j])) for j in range(d)]
        self.TransformNorm = np.sqrt(L.prod())   #/(2*np.pi)**d)
        self.Spectrum = RandomField.Covariance.precompute_Spectrum(self.Frequences)



#######################################################################################################
#	Fourier Transform (FFTW)
#######################################################################################################
### - Only stationary covariance
### - Uses the Fastest Fourier Transform on the West

class Sampling_FFTW(Sampling_method_freq):

    def __init__(self, RandomField):
        super().__init__(RandomField)

        import pyfftw
        shpR = RandomField.ext_grid_shape
        shpC = shpR.copy()
        shpC[-1] = int(shpC[-1] // 2)+1
        axes = np.arange(self.ndim)
        flags=('FFTW_MEASURE', 'FFTW_DESTROY_INPUT', 'FFTW_UNALIGNED')
        self.fft_x     = pyfftw.empty_aligned(shpR, dtype='float64')
        self.fft_y 	   = pyfftw.empty_aligned(shpC, dtype='complex128')
        self.fft_plan  = pyfftw.FFTW(self.fft_x, self.fft_y, axes=axes, direction='FFTW_FORWARD',  flags=flags)
        self.ifft_plan = pyfftw.FFTW(self.fft_y, self.fft_x, axes=axes, direction='FFTW_BACKWARD', flags=flags)
        self.Spectrum_half = self.Spectrum[...,:shpC[-1]] * np.sqrt(self.Nd.prod())

    def __call__(self, noise):
        self.fft_x[:] = noise
        self.fft_plan()
        self.fft_y[:] *= self.Spectrum_half 
        self.ifft_plan()
        return self.fft_x[self.DomainSlice] / self.TransformNorm


#######################################################################################################
#	Vector Field Fourier Transform (VF_FFTW)
#######################################################################################################
### - Random vector fields
### - Only stationary covariance
### - Uses the Fastest Fourier Transform in the West

class Sampling_VF_FFTW(Sampling_method_freq):

    def __init__(self, RandomField):
        super().__init__(RandomField)

        import pyfftw
        n_cpu = 4
        try:
            n_cpu = int(os.environ['OMP_NUM_THREADS'])
        except:
            pass
        shpR = RandomField.ext_grid_shape
        shpC = shpR.copy()
        shpC[-1] = int(shpC[-1] // 2)+1
        axes = np.arange(self.ndim)
        flags=('FFTW_MEASURE', 'FFTW_DESTROY_INPUT', 'FFTW_UNALIGNED')
        self.fft_x     = pyfftw.empty_aligned(shpR, dtype='float32')
        self.fft_y 	   = pyfftw.empty_aligned(shpC, dtype='complex64')
        self.fft_plan  = pyfftw.FFTW(self.fft_x, self.fft_y, axes=axes, direction='FFTW_FORWARD',  flags=flags, threads=n_cpu)
        self.ifft_plan = pyfftw.FFTW(self.fft_y, self.fft_x, axes=axes, direction='FFTW_BACKWARD', flags=flags, threads=n_cpu)
        self.Spectrum_half = self.Spectrum[...,:shpC[-1]] * np.sqrt(self.Nd.prod())
        self.hat_noise = np.stack([np.zeros(shpC, dtype='complex64') for _ in range(3)], axis=-1)
        self.shpC = shpC

    def __call__(self, noise):
        tmp = np.zeros(noise.shape)
        for i in range(noise.shape[-1]):
            self.fft_x[:] = noise[...,i]
            self.fft_plan()
            self.hat_noise[...,i] = self.fft_y[:]
        self.hat_noise  = np.einsum('kl...,...l->...k' , self.Spectrum_half, self.hat_noise)
        for i in range(noise.shape[-1]):
            self.fft_y[:] = self.hat_noise[...,i]
            self.ifft_plan()
            tmp[...,i] = self.fft_x[:]
        return tmp[self.DomainSlice] / self.TransformNorm


#######################################################################################################
#	Fourier Transform
#######################################################################################################
### - Only stationary covariance
### - Uses sscipy.fftpack (non the fastest solution)

class Sampling_FFT(Sampling_method_freq):

    def __init__(self, RandomField):
        super().__init__(RandomField)

    def __call__(self, noise):
        # noise_hat = FourierOfGaussian(noise)
        noise_hat = fft.ifftn(noise)
        y = self.Spectrum * noise_hat
        y = fft.fftn(y)
        return y.real[self.DomainSlice] / self.TransformNorm



#######################################################################################################
#	Sine Transform
#######################################################################################################
### - Only stationary covariance
### - Uses sscipy.fftpack (non the fastest solution)

class Sampling_DST(Sampling_method_freq):

    def __init__(self, RandomField):
        super().__init__(RandomField)

    def __call__(self, noise):
        y = self.Spectrum * noise
        for j in range(self.ndim):
            y = fft.dst(y, axis=j, type=1)
        return y[self.DomainSlice] / self.TransformNorm



#######################################################################################################
#	Cosine Transform
#######################################################################################################
### - Only stationary covariance
### - Uses sscipy.fftpack (non the fastest solution)

class Sampling_DCT(Sampling_method_freq):

    def __init__(self, RandomField):
        super().__init__(RandomField)

    def __call__(self, noise):
        y = self.Spectrum * noise 
        for j in range(self.ndim):
            y = fft.dct(y, axis=j, type=2)
        return y[self.DomainSlice] / self.TransformNorm



#######################################################################################################
#	Non-Uniform FFT
#######################################################################################################
### - Only stationary covariance
### - Non-regular grid

class Sampling_NFFT(Sampling_method_freq):

    def __init__(self, RandomField):
        super().__init__(RandomField)

        # Prepare NFFT objects
        from pynfft.nfft import NFFT
        x = RandomField.nodes
        M = x.shape[1]
        self.nfft_obj = NFFT(self.Nd, M)#, n=self.Nd, m=1)
        self.nfft_obj.x = (x - 0.5) / np.tile(self.L, [M, 1]).T
        self.nfft_obj.precompute()


    def __call__(self, noise):
        self.nfft_obj.f_hat = self.Spectrum * FourierOfGaussian(noise)
        y = self.nfft_obj.trafo()
        # assert (abs(y.imag) < 1.e-8).all(), np.amax(abs(y.imag))
        y = np.array(y.real, dtype=np.float) / self.TransformNorm
        return y



#######################################################################################################
#	Hierarchical matrix H2 (Fast Multipole)
#######################################################################################################
### - Non-stationnary covariances
### - Non-regular grid

class Sampling_H2(Sampling_method_base):

    def __init__(self, RandomField, lib=METHOD_H2, **kwargs):
        super().__init__(RandomField)
        L, Nd, d = self.L, self.Nd, self.ndim
        if self.verbose: print('\nSetting up H2-matrix...\n')
        
        t0 = time()

        axes = (np.arange(N).astype(np.float)/N,)*d
        position = np.meshgrid(*axes)
        position = np.vstack(list(map(np.ravel, position)))
        nvoxels = position.shape[1]

        if lib in (METHOD_H2_hlibpro,):
            import os
            os.system("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/khristen/ThirdPartyCode/hlibpro-2.7.2/lib")
            os.system("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/khristen/ThirdPartyCode/hlibpro-2.7.2/aux/lib") 
            from source.hlibpro.wrapper_Covariance import H2matrix
            self.M_H2 = H2matrix(position)

        elif lib in (METHOD_H2_h2tools, METHOD_H2):    
            from h2tools import Problem, ClusterTree
            from h2tools.mcbh import mcbh
            from h2tools.collections import particles

            Covariance = RandomField.Covariance
            nu, corrlen = Covariance.nu, Covariance.corrlen
            nu_mod = (nu-d/2)/2
            sigma = sqrt( gamma(nu + d/2) / gamma(nu) * (nu/(2*pi))**(d/2) / np.prod(corrlen) )
            sigma *= gamma(nu_mod)/gamma(nu_mod + d/2)

            if 'angle_field' in kwargs.keys():
                anis_angle = kwargs['angle_field']
            else:
                anis_angle = Covariance.angle_anis

            from source.Covariance.wrapper_Covariance import py_Matern_block_func
            def block_func(row_data, rows, col_data, cols):
                submatrix = sigma * py_Matern_block_func(row_data, rows, col_data, cols, nu_mod, corrlen * sqrt(nu_mod/nu), anis_angle, d) / N**(d/2)
                return submatrix
            
            data = particles.Particles(d, nvoxels, position)
            tree = ClusterTree(data, block_size=10)
            problem = Problem(block_func, tree, tree, symmetric=1, verbose=self.verbose)
            self.M_H2 = mcbh(problem, tau=1e-1, iters=0, verbose=self.verbose, random_init=0)

        print('Build H2:', time()-t0)

        ### Test
        if self.verbose:
            import matplotlib.pyplot as plt
            N = self.Nd[0]
            k = N
            r = np.arange(k)/N
            x = np.zeros(nvoxels)
            # x[0] = 1
            x[:k] = Covariance.eval_sqrt(r)/N**(d/2)
            t0 = time()
            Mx = self.M_H2.dot(x)
            print('dot', time()-t0)
            plt.figure()
            plt.plot(r, Mx[:k], 'o-')
            plt.plot(r, Covariance.eval(r), 'o--')


    def __call__(self, noise):
        z = noise.flatten()
        y = self.M_H2.dot(z)
        y = y.reshape(self.Nd)
        return y[self.DomainSlice]


#######################################################################################################
#	ODE-basd method
#######################################################################################################
### - Non-stationnary covariances
### - Non-regular grid

class Sampling_ODE(Sampling_method_base):

    def __init__(self, RandomField, lib="h2tools", **kwargs):
        super().__init__(RandomField)
        L, Nd, d = self.L, self.Nd, self.ndim

        from source.ODE_based.TransientPower import Problem, TransientPower
        from source.CovarianceKernels import set_ShapeOperator

        Covariance = RandomField.Covariance
        nu, corrlen = Covariance.nu, Covariance.corrlen

        if 'angle_field' in kwargs.keys():
            anis_angle = kwargs['angle_field']
        else:
            anis_angle = Covariance.angle_anis

        coef, detTheta = set_ShapeOperator(corrlen, anis_angle, ndim=self.ndim)
        coef /= 2*nu

        t0 = time()
        self.pb = Problem(N, d, coef)
        if self.verbose: print('Build problem time', time()-t0)
        self.tpow = TransientPower(self.pb)
        alpha = (nu + d/2)/2
        eta = sqrt( gamma(nu + d/2)/gamma(nu) * ((4*pi)/(2*nu))**(d/2) * sqrt(detTheta) )
        slope = 1.e6 * detTheta**(1/d)/(2*nu)
        self.apply = lambda x: self.tpow(alpha, x, nts=10, theta=4, slope=40) * eta * np.sqrt(self.Nd.prod())


    def __call__(self, noise):
        y = self.apply(noise.flatten())
        y = y.reshape(self.Nd)
        return y[self.DomainSlice]




#######################################################################################################
#	Best Rational Approximation
#######################################################################################################
### - Non-stationnary covariances
### - Non-regular grid
### - Opimal method

class Sampling_Rational(Sampling_method_base):

    def __init__(self, RandomField, **kwargs):
        super().__init__(RandomField)
        L, Nd, d = self.L, self.Nd, self.ndim

        from source.RationalApproximation import RationalApproximation
        from source.ODE_based.TransientPower import Problem
        from source.CovarianceKernels import set_ShapeOperator

        Covariance = RandomField.Covariance
        nu, corrlen = Covariance.nu, Covariance.corrlen

        if 'angle_field' in kwargs.keys():
            anis_angle = kwargs['angle_field']
        else:
            anis_angle = Covariance.angle_anis

        coef, detTheta = set_ShapeOperator(corrlen, anis_angle, ndim=self.ndim)
        coef /= 2*nu

        t0 = time()
        self.pb = Problem(N, d, coef)
        if self.verbose: print('Assemble problem time', time()-t0)
        alpha = (nu + d/2)/2
        eta = sqrt( gamma(nu + d/2)/gamma(nu) * ((4*pi)/(2*nu))**(d/2) * sqrt(detTheta) )
        self.RA = RationalApproximation(self.pb, alpha, niter=4)
        self.apply = lambda x: self.RA(x) * eta * np.sqrt(self.Nd.prod())


    def __call__(self, noise):
        y = self.apply(noise.flatten())
        y = y.reshape(self.Nd)
        return y[self.DomainSlice]


#######################################################################################################
#	Best Rational Approximation for Von Karman wind field with blocking
#######################################################################################################
### - Non-stationary covariances
### - Optimal method
        
class Sampling_Rational_VK_Wind_Blocking(Sampling_method_base):
    '''
    This function assumes \mcL u = -\nabla\cdot ( L(x)^2 \nabla u )
    '''
    def __init__(self, RandomField, **kwargs):
        super().__init__(RandomField)

        from fde_solve import fde_solve

        corrlen = kwargs['correlation_length']
        # eps = kwargs['viscous_dissipation_rate']
        # C   = kwargs['kolmogorov_constant']
        self.corrlen = corrlen

        self.z_grid = kwargs.get('z_grid', None)
        if self.z_grid is not None:
            self.h = np.diff(self.z_grid)
        else:
            self.h = self.L[2] / self.Nd[2] * np.ones(self.Nd[2]-1)

        L, Nd = self.L, self.Nd
        self.Frequencies = [(2*pi/L[j])*(Nd[j]*fft.fftfreq(Nd[j])) for j in range(3)]
        self.TransformNorm = np.sqrt(L[0]*L[1]) / (2*pi)
        
        # NOTE: the following values are fixed for the VK model
        self.corr_len_fun = lambda z: corrlen**2 #* np.tanh(z)
        self.alpha = 17/12

        # self.factor = np.sqrt( C * (eps**(2/3)) * (corrlen**(17/3)) / (4*np.pi) )
        self.factor = np.sqrt( kwargs['E0'] / (4*np.pi) )
        # self.factor = np.sqrt( C * (eps**(2/3)) / (4*np.pi) )

        # instantiate rational approximation object
        self.fde_solve = fde_solve(Nd[2], self.alpha, self.corr_len_fun, domain_height=L[2], z_grid=self.z_grid)

    ### z-derivative
    def Dz(self, f, adjoint=False):
        h = self.h
        dzf = np.zeros_like(f)
        if adjoint:
            dzf[...,0] = 0
            dzf[...,1] = (f[...,0]/h[0] - f[...,2]/(h[1] + h[2]))
            dzf[...,-2] = f[...,-3]/(h[-2] + h[-3])
            dzf[...,-1] = f[...,-2]/(h[-1] + h[-2])
            dzf[...,2:-2] = -(f[...,3:-1]/(h[2:-1] + h[3:]) - f[...,1:-3]/(h[:-3] + h[1:-2]))
        else:
            dzf[:,:,2:-1] = (f[:,:,3:] - f[:,:,1:-2]) / (h[2:] + h[1:-1])
            dzf[:,:,0] = f[:,:,1] / h[0]
            dzf[:,:,1] = f[:,:,2] / (h[1] + h[0])
            dzf[:,:,-1] = 0
        return dzf


    ### Apply curl
    def curl(self, f_hat, adjoint=False):
        k1, k2, _ = np.meshgrid(*self.Frequencies, indexing='ij')

        tmp1 = f_hat[...,0].copy()
        tmp2 = f_hat[...,1].copy()
        tmp3 = f_hat[...,2].copy()
        if not adjoint:
            dzf1 = self.Dz(f_hat[...,0], adjoint=False)
            dzf2 = self.Dz(f_hat[...,1], adjoint=False)
            f_hat[...,0] = 1j*k2*tmp3 - dzf2
            f_hat[...,1] = dzf1       - 1j*k1*tmp3
            f_hat[...,2] = 1j*k1*tmp2 - 1j*k2*tmp1
        else:
            dzf1 = self.Dz(f_hat[...,0], adjoint=True)
            dzf2 = self.Dz(f_hat[...,1], adjoint=True)
            f_hat[...,0] = 1j*k2*tmp3 + dzf2
            f_hat[...,1] = -dzf1      - 1j*k1*tmp3
            f_hat[...,2] = 1j*k1*tmp2 - 1j*k2*tmp1

        return f_hat


    def __call__(self, noise, **kwargs):
        Robin_const = [np.infty, np.infty, kwargs.get('Robin_const')] ### 3rd component is None if is not in kwargs

        f_hat = noise / np.sqrt(self.L[2] / self.Nd[2])

        ## Fourier transform of RHS in x- and y-coordinates
        f_hat = fft.fft(fft.fft(f_hat, axis = 1), axis = 0)

        ### STEP 1

        # define frequencies in 2D domain (x & y)
        k1, k2 = np.meshgrid(*self.Frequencies[:2], indexing='ij')

        # solve rational approx problem (in z) for each k1, k2
        t0=time()

        rhs = np.zeros(f_hat.shape[2], dtype=f_hat.dtype)
        for l in range(3):
            for i in range(self.Nd[0]):
                for j in range(self.Nd[1]):
                    rhs[:] = f_hat[i,j,:,l]
                    f_hat[i,j,:,l] = self.fde_solve(rhs, k1[i,j], k2[i,j], Robin_const=Robin_const[l], component=l, t=0, adjoint=adjoint, jac=jac, grad_coef=grad_coef)

        # Apply curl
        f_hat = self.curl(f_hat, adjoint=False)

        # Inverse Fourier transform of RHS in x- and y-coordinates
        f = fft.ifft(fft.ifft(f_hat, axis = 1), axis = 0)
        f = self.factor * f.real[self.DomainSlice]# / self.TransformNorm

        return f 


#######################################################################################################
#	Best Rational Approximation for rapid distortion wind field with blocking
#######################################################################################################
### - Non-stationary covariances
### - Optimal method

class Sampling_Rational_Rapid_Distortion_Wind_Blocking(Sampling_Rational_VK_Wind_Blocking):
    '''
    Solves (1 -\nabla\cdot ( \Theta_t \nabla ))^\alpha u = f_t
    '''
    def __init__(self, RandomField, **kwargs):
        super().__init__(RandomField, **kwargs)

        from fde_solve import fde_solve

        self.alpha = 17/12
        self.t = 1.0

        # instantiate rational approximation object
        L, Nd = self.L, self.Nd
        self.fde_solve = fde_solve(Nd[2], self.alpha, self.corr_len_fun, domain_height=L[2], t=self.t)

    def __call__(self, noise, **kwargs):
        Robin_const = [np.infty, np.infty, kwargs.get('Robin_const')] ### 3rd component is None if is not in kwargs
        t = kwargs.get('t', 1.0)

        ## Fourier transform of RHS in x- and y-coordinates
        f_hat = noise
        f_hat = fft.fft(fft.fft(f_hat, axis = 1), axis = 0)

        ## Correlate noise
        f_hat = self.distort_noise(f_hat, t)

        # define frequencies in 2D domain (x & y)
        k1, k2 = np.meshgrid(*self.Frequencies[:2], indexing='ij')

        # solve rational approx problem (in z) for each k1, k2
        t0=time()

        rhs = np.zeros(f_hat.shape[2], dtype=f_hat.dtype)
        for l in range(3):
            for i in range(self.Nd[0]):
                for j in range(self.Nd[1]):
                    rhs[:] = f_hat[i,j,:,l]
                    f_hat[i,j,:,l] = self.fde_solve(rhs, k1[i,j], k2[i,j], Robin_const=Robin_const[l], component=l, t=t, adjoint=False, jac=False, grad_coef=False)

        # Apply curl
        f_hat = self.curl(f_hat, adjoint=False)

        # Inverse Fourier transform of RHS in x- and y-coordinates
        f = fft.ifft(fft.ifft(f_hat, axis = 1), axis = 0)
        f = self.factor * f.real[self.DomainSlice]# / self.TransformNorm

        return f

    def distort_noise(self, f_hat, beta):

        f_hat = fft.fft(f_hat, axis = 2)

        k = np.array(list(np.meshgrid(*self.Frequencies, indexing='ij')))

        with np.errstate(divide='ignore', invalid='ignore'):

            k1  = k[0,...]
            k2  = k[1,...]
            k3  = k[2,...]
            k30  = k3 + beta*k1

            kk = k1**2 + k2**2 + k3**2
            kk0 = k1**2 + k2**2 + k30**2

            #### RDT

            s = k1**2 + k2**2
            C1  =  beta * k1**2 * (kk0 - 2 * k30**2 + beta * k1 * k30) / (kk * s)
            tmp =  beta * k1 * np.sqrt(s) / (kk0 - k30 * k1 * beta)
            C2  =  k2 * kk0 / s**(3/2) * np.arctan (tmp)

            zeta1_by_zeta3 =  (C1 - k2/k1 * C2)*kk/kk0
            zeta2_by_zeta3 =  (k2/k1 *C1 + C2)*kk/kk0
            one_by_zeta3 =  kk/kk0

            # deal with divisions by zero (2/2)
            zeta1_by_zeta3 = np.nan_to_num(zeta1_by_zeta3)
            zeta2_by_zeta3 = np.nan_to_num(zeta2_by_zeta3)
            one_by_zeta3 = np.nan_to_num(one_by_zeta3)

            f_hat[...,2] = -f_hat[...,0]*zeta1_by_zeta3 - f_hat[...,1]*zeta2_by_zeta3 + f_hat[...,2]*one_by_zeta3

            return fft.ifft(f_hat, axis = 2)



#######################################################################################################


       

        
            