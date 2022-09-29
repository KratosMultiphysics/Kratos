import numpy as np

# TODO Remove 'return' in function names

def normal(mean : float, deviation : float, shape : tuple = None):
    """
    Draws samples from the normal distribution defined by arguments.

    Input parameters:
    - mean: float, mean of the normal distribution.
    - deviation: float, standard deviation of the normal distribution.
    - shape (optional): tuple, shape of the desired output. If omitted, the output is a float.

    Output: float if shape is omitted. Otherwise, nested list of required shape, whose every element is a float.
    These floats are drawn independantly from the normal distribution specified by the input arguments.

    Note: differences with returnStandardNormal.
    - random number generator: this method uses numpy.random.Generator.normal, whereas returnStandardNormal uses numpy.random.normal.
    - output type and shape: this method returns either a float or a nested list of floats of arbitrary shape, whereas returnStandardNormal returns a list containing a single float.
    """

    # Create random number generator
    rng = np.random.default_rng()
    # Draw samples
    sample = rng.normal(mean, deviation, shape)
    # Ensure output type
    if isinstance(sample,np.ndarray):
        sample = sample.astype(float).tolist()
    return sample


def returnStandardNormal(*args):
    """
    Return a normally-distributed variable with given parameters
    """
    return [np.random.normal(args[0], args[1])]

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

def returnUniformAndNormal(*args):
    """
    Return one integer uniformly distributed random variable and
    one normal random variable
    """
    return [int(np.random.uniform(args[0],args[1])),np.random.normal(args[2], args[3], 1)]

def returnUniformAndTwoNormal(*args):
    """
    Return one integer uniformly distributed random variable and
    two normal random variables
    """
    return [int(np.random.uniform(args[0],args[1])),np.random.normal(args[2], args[3], 1),np.random.normal(args[4], args[5], 1)]

def returnIntegerUniform(*args):
    """
    Return one integer uniformly distributed random variable
    """
    return [int(np.random.uniform(args[0],args[1]))]

def returnIntegerUniformAndNormalAndUniform(*args):
    """
    Return one integer uniformly distributed random variable,
    one normal random variable and
    one uniformly distributed random variable.
    """
    return [int(np.random.uniform(args[0],args[1])),np.random.normal(args[2], args[3], 1),np.random.uniform(args[4], args[5], 1)]