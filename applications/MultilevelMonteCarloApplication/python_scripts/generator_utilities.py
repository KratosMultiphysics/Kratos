import numpy as np

# TODO: move in an external file and create a generator class, e.g. MonteCarlo.generator.Generate()

'''
function generatng the random variable
'''
def GenerateSample():
    alpha = 2.0
    beta = 6.0
    number_samples = 1
    sample = list(np.random.beta(alpha,beta,number_samples))
    return sample