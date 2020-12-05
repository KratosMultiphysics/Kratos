from math import *
import numpy as np
import scipy.fftpack as fft

class GenerateCN:

    def __init__(self, corrlen, grid_dimensions=None, grid_shape=None, grid_level=4, ndim=2, const=1e2, **kwargs):

        self.corrlen = corrlen
        self.ndim = ndim
        self.const=const

        if grid_dimensions is None:
            self.grid_dimensions = [1] * ndim
        else:
            self.grid_dimensions = grid_dimensions

        if grid_shape is None: ### Unit square (cube)
            N = 2**grid_level + 1
            self.grid_shape = [N] * ndim
        elif np.isscalar(grid_shape):
            self.grid_shape = [grid_shape] * ndim
        else:
            self.grid_shape = grid_shape[:ndim]

        self.grid_dimensions = np.array(self.grid_dimensions)
        self.grid_shape = np.array(self.grid_shape)

        self.Frequencies = [(2*pi/self.grid_dimensions[j])*(self.grid_shape[j]*fft.fftfreq(self.grid_shape[j])) for j in range(ndim)]
        self.noise_std = np.sqrt(np.prod(self.grid_dimensions/self.grid_shape))
        self.prng = np.random.RandomState()

    def __call__(self, seed=None):

        if seed is not None:
            self.prng = np.random.RandomState(seed=seed)

        if self.ndim == 3:
            # GENERATE NOISE
            noise = np.stack([self.prng.normal(0, 1, self.grid_shape) for _ in range(3)], axis=-1)
            noise *= self.noise_std
            noise = fft.fftn(noise,axes=(0,1,2)) # DFT FORWARD

            # GENERATE FREQUENCIES
            k = np.array(list(np.meshgrid(*self.Frequencies, indexing='ij')))
            k1 = k[0,...]
            k2 = k[1,...]
            k3 = k[2,...]
            kk = np.sum(k**2,axis=0)

            # APPLY SPECTRUM
            spectrum = self.const / (1 + (self.corrlen**2) * kk)**(17/12)
            for i in range(3):
                noise[...,i] *= spectrum

            # APPLY CURL
            tmp1 = noise[...,0].copy()
            tmp2 = noise[...,1].copy()
            tmp3 = noise[...,2].copy()
            noise[...,0] = k2*tmp3 - k3*tmp2
            noise[...,1] = k3*tmp1 - k1*tmp3
            noise[...,2] = k1*tmp2 - k2*tmp1
            vf = fft.ifftn(1j*noise,axes=(0,1,2)).real # DFT BACKWARD

        else:
            # GENERATE NOISE
            noise = self.prng.normal(0, 1, self.grid_shape)
            noise *= self.noise_std
            noise = fft.fftn(noise) # DFT FORWARD

            # GENERATE FREQUENCIES
            k = np.array(list(np.meshgrid(*self.Frequencies, indexing='ij')))
            k1 = k[0,...]
            k2 = k[1,...]
            kk = np.sum(k**2,axis=0)

            # APPLY SPECTRUM
            spectrum = self.const / (1 + (self.corrlen**2) * kk)**(17/12)
            noise *= spectrum

            # APPLY CURL
            vf = np.zeros((noise.shape + (2,)),dtype=np.complex128)
            vf[...,0] = k2*noise
            vf[...,1] =-k1*noise
            vf = fft.ifftn(1j*vf,axes=(0,1)).real # DFT BACKWARD

        return vf


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    vf_gen = GenerateCN(corrlen=0.1, grid_level=8)
    vf = vf_gen()

    fig, (ax1, ax2) = plt.subplots(ncols=2)

    ax1.set_title("x-component")
    im = ax1.imshow(vf[:,:,0].T, origin='lower')

    ax2.set_title("y-component")
    im = ax2.imshow(vf[:,:,1].T, origin='lower')

    plt.show()




