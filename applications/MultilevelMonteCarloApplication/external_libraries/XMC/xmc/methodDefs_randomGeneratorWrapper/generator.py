import numpy as np

def returnStandardNormal(*args):
    """
    Return a normally distributed variable with given parameters
    """
    return [np.random.normal(args[0], args[1], 1)]

def returnStandardNormalandUniform(*args):
    """
    Return one standard normal random variable and one uniformly distributed
    random variable in the interval [0,1]
    """
    return [np.random.normal(0, 1, 1), np.random.uniform(0,1,1)]

def returnBeta(*args):
    """
    Return a beta distributed random variables with given parameters
    """
    return [np.random.beta(args[0],args[1],1)]

def returnUniform(*args):
    """
    Return one uniformly distributed random variable in the interval [0,1]
    """
    return [np.random.uniform(args[0], args[1], 1)]

def returnRayleigh(*args):
    """
    Return one rayleigh random variable with given parameters
    use numpy.random.rayleigh(scale, size)
    already returns as a numpy array
    """
    return np.random.rayleigh(args[0],size=1)

def returnRayleighAndUniform(*args):
    """
    Return one rayleigh and uniformly distributed random variables
    with given parameters
    """
    rayleigh_rv = np.random.rayleigh(args[0],size=1)[0]
    uniform_rv = np.random.uniform(args[1], args[2], 1)
    return [rayleigh_rv, uniform_rv]