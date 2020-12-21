from xmc.distributedEnvironmentFramework import *

def centralMomentErrorWrapper(dimension,order,*args):
    if dimension==0:
        if order==1:
            return computeCentralMomentsErrorEstimationOrderOneDimensionZero(*args)
        elif order==2:
            return computeCentralMomentsErrorEstimationOrderTwoDimensionZero(*args)
        else:
            raise ValueError('Moments of order > 2 are not supported yet for dimension 0.')
    elif dimension==1:
        if order==1:
            return computeCentralMomentsErrorEstimationOrderOneDimensionOne(*args)
        else:
            raise ValueError('Moments of order > 1 are not supported yet for dimension 1.')
    else:
        raise ValueError('Index set of dimension > 2 are not supported yet.')

@ExaquteTask(returns=1)
def centralMomentErrorWrapper_Task(*args):
    return centralMomentErrorWrapper(*args)


def centralCombinedMomentErrorWrapper(dimension,order,*args):
    if dimension==0:
        if order==1:
            return computeCentralMomentsErrorEstimationOrderOneDimensionZero(*args)
        else:
            raise ValueError('Moments of order > 1 are not supported yet for dimension 0.')
    else:
        raise ValueError('Index set of dimension > 0 are not supported yet.')

@ExaquteTask(returns=1)
def centralCombinedMomentErrorWrapper_Task(*args):
    return centralCombinedMomentErrorWrapper(*args)


def computeCentralMomentsErrorEstimationOrderOneDimensionZero(power_sum_1,power_sum_2,number_samples):
    central_moment_2 = (number_samples*power_sum_2-power_sum_1**2) / ((number_samples-1)*number_samples)
    estimation_error_1 = central_moment_2/number_samples
    return estimation_error_1


def computeCentralMomentsErrorEstimationOrderTwoDimensionZero(power_sum_1,power_sum_2,power_sum_3,power_sum_4,number_samples):
    central_moment_2 = (number_samples*power_sum_2-power_sum_1**2) / ((number_samples-1)*number_samples)
    central_moment_4 = ((-4*number_samples**2+8*number_samples-12)*power_sum_3*power_sum_1+ \
    (number_samples**3-2*number_samples**2+3*number_samples)*power_sum_4+ \
    6*number_samples*power_sum_2*power_sum_1**2+(9-6*number_samples)*power_sum_2**2-3*power_sum_1**4) / \
    ((number_samples-3)*(number_samples-2)*(number_samples-1)*number_samples)
    estimation_error_2 = central_moment_4/number_samples - \
        (central_moment_2**2)*(number_samples-3)/(number_samples*(number_samples-1))
    return estimation_error_2


def computeCentralMomentsErrorEstimationOrderOneDimensionOne(power_sum_10,power_sum_01,power_sum_20,power_sum_11,power_sum_02,number_samples):
    partial_central_moment_2 = (number_samples*power_sum_02-power_sum_01**2) / ((number_samples-1)*number_samples)
    partial_estimation_error_1 = partial_central_moment_2/number_samples
    return partial_estimation_error_1
