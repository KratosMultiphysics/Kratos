import numpy as np

# TODO: move in an external file and create a generator class, e.g. MonteCarlo.generator.Generate()

'''
function generatng the random variable
'''
def GenerateSample():

    # Poisson equation
    # alpha = 2.0
    # beta = 6.0
    # number_samples = 1
    # sample = np.random.beta(alpha,beta,number_samples)
    # return sample

    # Compressible potential flow problem
    # sample = []
    # mean_Mach = 0.3
    # std_deviation_Mach = 0.01
    # number_samples = 1
    # sample.append(np.random.normal(mean_Mach,std_deviation_Mach,number_samples))
    # mean_angle_attack = 0.0 # [rad] = 0 [degrees] airfoil already has 5 degrees
    # std_deviation_angle_attack = np.deg2rad(0.1)
    # sample.append(np.random.normal(mean_angle_attack,std_deviation_angle_attack,number_samples))
    # if sample[0] >= 1.0 or sample[0] <= 0.0 :
    #     raise Exception ("stochastic Mach number computed > 1 or < 0")

    # Problem zero problem
    alpha0 = 0.12 # corresponds to open terrain from annexure to Euro Code DIN EN 1991-1-4 NA
    salpha = 0.012 # a standard deviation of 10 % is considered
    sample = [np.random.normal(alpha0,salpha)]

    return sample
