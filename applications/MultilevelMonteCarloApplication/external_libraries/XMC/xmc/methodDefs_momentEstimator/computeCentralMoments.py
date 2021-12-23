from exaqute import *

def centralMomentWrapper(dimension,order,*args):
    if dimension==0:
        if order==1:
            return computeCentralMomentsOrderOneDimensionZero(*args)
        elif order==2:
            return computeCentralMomentsOrderTwoDimensionZero(*args)
        elif order==3:
            return computeCentralMomentsOrderThreeDimensionZero(*args)
        elif order==4:
            return computeCentralMomentsOrderFourDimensionZero(*args)
        else:
            raise ValueError('Moments of order > 4 are not supported yet for dimension 0.')
    elif dimension==1:
        if order==1:
            return computeCentralMomentsOrderOneDimensionOne(*args)
        elif order==2:
            return computeCentralMomentsOrderTwoDimensionOne(*args)
        else:
            raise ValueError('Moments of order > 2 are not supported yet for dimension 1.')
    else:
        raise ValueError('Index sets of dimension > 1 are not supported yet.')

@task(keep=True, returns=1)
def centralMomentWrapper_Task(*args):
    return centralMomentWrapper(*args)

def centralCombinedMomentWrapper(dimension,order,*args):
    if dimension==0:
        if order==1:
            return computeCentralMomentsOrderOneDimensionZero(*args)
        elif order==2:
            return computeCentralMomentsOrderTwoDimensionZeroBiased(*args)
        else:
            raise ValueError('Moments of order > 2 are not supported yet for dimension 0.')
    else:
        raise ValueError('Index sets of dimension > 0 are not supported yet.')

@task(keep=True, returns=1)
def centralCombinedMomentWrapper_Task(*args):
    return centralCombinedMomentWrapper(*args)


def computeCentralMomentsOrderOneDimensionZero(*args):
    power_sum_1 = args[0] ; number_samples = args[-1]
    central_moment_1 = power_sum_1/number_samples
    return central_moment_1

def computeCentralMomentsOrderTwoDimensionZero(*args):
    power_sum_1 = args[0] ; power_sum_2 = args[1] ; number_samples = args[-1]
    central_moment_2 = (number_samples*power_sum_2-power_sum_1**2) / ((number_samples-1)*number_samples)
    return central_moment_2


def computeCentralMomentsOrderThreeDimensionZero(*args):
    power_sum_1 = args[0] ; power_sum_2 = args[1] ; power_sum_3 = args[2] ; number_samples = args[-1]
    central_moment_3 = (number_samples**2*power_sum_3-3*number_samples*power_sum_2*power_sum_1+2*power_sum_1**3) / \
    ((number_samples-2)*(number_samples-1)*number_samples)
    return central_moment_3


def computeCentralMomentsOrderFourDimensionZero(*args):
    power_sum_1 = args[0] ; power_sum_2 = args[1] ; power_sum_3 = args[2]
    power_sum_4 = args[3] ; number_samples = args[-1]
    central_moment_4 = ((-4*number_samples**2+8*number_samples-12)*power_sum_3*power_sum_1+ \
    (number_samples**3-2*number_samples**2+3*number_samples)*power_sum_4+ \
    6*number_samples*power_sum_2*power_sum_1**2+(9-6*number_samples)*power_sum_2**2-3*power_sum_1**4) / \
    ((number_samples-3)*(number_samples-2)*(number_samples-1)*number_samples)
    return central_moment_4


def computeCentralMomentsOrderTwoDimensionZeroBiased(*args):
    power_sum_1 = args[0] ; power_sum_2 = args[1] ; number_samples = args[-1]
    central_moment_2 = (number_samples*power_sum_2-power_sum_1**2) / ((number_samples)**2)
    return central_moment_2


def computeCentralMomentsOrderOneDimensionOne(*args):
    power_sum_10 = args[0] ; power_sum_01 = args[1] ; number_samples = args[-1]
    partial_central_moment_1 = power_sum_01/number_samples
    return partial_central_moment_1


def computeCentralMomentsOrderTwoDimensionOne(*args):
    power_sum_10 = args[0] ; power_sum_01 = args[1] ; power_sum_20 = args[2]
    power_sum_11 = args[3] ; power_sum_02 = args[4] ; number_samples = args[-1]
    partial_central_moment_2 = (number_samples*power_sum_02-power_sum_01**2) / ((number_samples-1)*number_samples)
    return partial_central_moment_2
