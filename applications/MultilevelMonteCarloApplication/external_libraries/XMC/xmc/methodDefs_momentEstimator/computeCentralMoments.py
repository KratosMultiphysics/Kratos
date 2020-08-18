# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

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
        else:
            raise ValueError('Moments of order > 0 are not supported yet for dimension 1.')
    else:
        raise ValueError('Index sets of dimension > 0 are not supported yet.')

@ExaquteTask(returns=1)
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

@ExaquteTask(returns=1)
def centralCombinedMomentWrapper_Task(*args):
    return centralCombinedMomentWrapper(*args)


def computeCentralMomentsOrderOneDimensionZero(power_sum_1,number_samples):
    central_moment_1 = power_sum_1/number_samples
    return central_moment_1

def computeCentralMomentsOrderTwoDimensionZero(power_sum_1,power_sum_2,number_samples):
    central_moment_2 = (number_samples*power_sum_2-power_sum_1**2) / ((number_samples-1)*number_samples)
    return central_moment_2


def computeCentralMomentsOrderThreeDimensionZero(power_sum_1,power_sum_2,power_sum_3,number_samples):
    central_moment_3 = (number_samples**2*power_sum_3-3*number_samples*power_sum_2*power_sum_1+2*power_sum_1**3) / \
    ((number_samples-2)*(number_samples-1)*number_samples)
    return central_moment_3


def computeCentralMomentsOrderFourDimensionZero(power_sum_1,power_sum_2,power_sum_3,power_sum_4,number_samples):
    central_moment_4 = ((-4*number_samples**2+8*number_samples-12)*power_sum_3*power_sum_1+ \
    (number_samples**3-2*number_samples**2+3*number_samples)*power_sum_4+ \
    6*number_samples*power_sum_2*power_sum_1**2+(9-6*number_samples)*power_sum_2**2-3*power_sum_1**4) / \
    ((number_samples-3)*(number_samples-2)*(number_samples-1)*number_samples)
    return central_moment_4


def computeCentralMomentsOrderTwoDimensionZeroBiased(power_sum_1,power_sum_2,number_samples):
    central_moment_2 = (number_samples*power_sum_2-power_sum_1**2) / ((number_samples)**2)
    return central_moment_2


def computeCentralMomentsOrderOneDimensionOne(power_sum_10,power_sum_01,number_samples):
    partial_central_moment_1 = power_sum_01/number_samples
    return partial_central_moment_1


def computeCentralMomentsOrderTwoDimensionOne(power_sum_10,power_sum_01,power_sum_20,power_sum_11,power_sum_02,number_samples):
    partial_central_moment_2 = (number_samples*power_sum_02-power_sum_01**2) / ((number_samples-1)*number_samples)
    return partial_central_moment_2
