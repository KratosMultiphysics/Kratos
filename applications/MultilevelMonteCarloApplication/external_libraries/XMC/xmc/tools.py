import importlib
from operator import mul
from functools import reduce
from math import sqrt, log, factorial
from collections.abc import Iterable

from xmc.distributedEnvironmentFramework import *


def dynamicImport(fullName):
    """
    This function returns the requested attribute of a module; it manages the necessary imports.
    The input argument is the full name as a string, e.g. 'module.package.method'.
    Alternatively, the input argument may be a list of such full names. The output will then be a list of same size containing the requested attributes in the same order.
    This only works for modern Python (at least Python > 3).
    """
    if fullName is None:  # safety
        return None
    elif fullName is list or fullName is tuple:
        requested_attributes = []
        for name in fullName:
            requested_attributes.append(dynamicImport(name))
        return requested_attributes
    else:
        split_name = fullName.split(".")
        module_name = ".".join(split_name[0:-1])
        attribute_name = split_name[-1]
        module = importlib.import_module(module_name)
        requested_attribute = getattr(module, attribute_name)
        return requested_attribute


def instantiateObject(fullName, **kwargs):
    object_type = dynamicImport(fullName)
    return object_type(**kwargs)


def doNothing(*args, **kwArgs):
    """
    This is a function that does nothing. Useful to hold a place without taking any action.
    """
    pass


def returnInput(*args):
    """
    This is a function that returns exactly what is passed to it. Useful to hold a place in the
    dataflow without modifying the data.
    """
    if len(args) == 1:
        # A single argument has been passed and is now packaged inside a tuple because of *.
        # Unpack it:
        args = args[0]
    return args


@ExaquteTask(returns=1)
def returnInput_Task(*args):
    return returnInput(*args)


# For backward compatibility
convertObjectToFuture = returnInput_Task

# TODO Unused. remove?
def getUnionAndMap(listOfLists):
    """
    Given a list of lists, computes union of all lists, as well as map from list of lists
    to union lists.
    """
    union_list = []
    union_map = []
    for i in range(len(listOfLists)):
        union_map.append([None] * len(listOfLists[i]))
        for j in range(len(listOfLists[i])):
            found = False
            for k in range(len(union_list)):
                if listOfLists[i][j] == union_list[k]:
                    union_map[i][j] = k
                    found = True
                    break
            if found is False:
                union_list.append(listOfLists[i][j])
                union_map[i][j] = len(union_list) - 1
    return union_list, union_map


def splitOneListIntoTwo(inputList):
    """
    Function that takes a list, each of whose entries is a list
    containing two entries, and assembles two separate lists, one
    for each of these entries
    """
    list1 = []
    list2 = []
    for entry in inputList:
        list1.append(entry[0])
        list2.append(entry[1])
    return list1, list2


def mergeTwoListsIntoOne(list1, list2):
    """
    Function that takes two lists and assembles them into one list
    such that each entry of the final list is a list containing the
    two corresponding entries of the two input lists
    """
    if len(list1) != len(list2):
        raise ValueError("Input lists not of same length")

    output_list = [[list1[i], list2[i]] for i in range(len(list1))]
    return output_list


def strictlyPositiveBoundaryBooleans(indexSet):
    """
    Accepts a convex set of non-negative multi-indices and returns
    a list of booleans of the same size that states whether the
    corresponding index lies on the surface of the set where none
    of the indices have zero-values.
    """
    max_indices = [
        max([indexSet[i][j] for i in range(len(indexSet))]) for j in range(len(indexSet[0]))
    ]
    is_index_in_bias = [None for _ in indexSet]
    for i in range(len(indexSet)):
        index = indexSet[i]
        is_index_max = [index[i] == max_indices[i] for i in range(len(index))]
        is_index_nonzero = [index[i] != 0 for i in range(len(index))]
        is_index_in_bias[i] = any(is_index_max) and all(is_index_nonzero)
    return is_index_in_bias


# TODO is this function useful?
def summation(*args):
    return sum(args)


@ExaquteTask(returns=1)
def summation_Task(*args):
    return sum(args)


def splitList(listToSplit, num_sublists=1):
    """
    Function that splits a list into a given number of sublists
    Reference: https://stackoverflow.com/questions/752308/split-list-into-smaller-lists-split-in-half
    """
    length = len(listToSplit)
    return [
        listToSplit[i * length // num_sublists : (i + 1) * length // num_sublists]
        for i in range(num_sublists)
    ]


def normalInverseCDF(y0):
    """
    Computes inverse of normal distribution function (percent point function)
    Sources:
    https://github.com/scipy/scipy/blob/59347ae8b86bcc92c339efe213128f64ab6df98c/scipy/special/cephes/ndtri.c,
    https://stackoverflow.com/questions/41338539/how-to-calculate-a-normal-distribution-percent-point-function-in-python
    References:
    Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, in press.
    """

    def polevl(x, coef):
        accum = 0
        for c in coef:
            accum = x * accum + c
        return accum

    s2pi = 2.50662827463100050242e0

    P0 = [
        -5.99633501014107895267e1,
        9.80010754185999661536e1,
        -5.66762857469070293439e1,
        1.39312609387279679503e1,
        -1.23916583867381258016e0,
    ]

    Q0 = [
        1,
        1.95448858338141759834e0,
        4.67627912898881538453e0,
        8.63602421390890590575e1,
        -2.25462687854119370527e2,
        2.00260212380060660359e2,
        -8.20372256168333339912e1,
        1.59056225126211695515e1,
        -1.18331621121330003142e0,
    ]

    P1 = [
        4.05544892305962419923e0,
        3.15251094599893866154e1,
        5.71628192246421288162e1,
        4.40805073893200834700e1,
        1.46849561928858024014e1,
        2.18663306850790267539e0,
        -1.40256079171354495875e-1,
        -3.50424626827848203418e-2,
        -8.57456785154685413611e-4,
    ]

    Q1 = [
        1,
        1.57799883256466749731e1,
        4.53907635128879210584e1,
        4.13172038254672030440e1,
        1.50425385692907503408e1,
        2.50464946208309415979e0,
        -1.42182922854787788574e-1,
        -3.80806407691578277194e-2,
        -9.33259480895457427372e-4,
    ]

    P2 = [
        3.23774891776946035970e0,
        6.91522889068984211695e0,
        3.93881025292474443415e0,
        1.33303460815807542389e0,
        2.01485389549179081538e-1,
        1.23716634817820021358e-2,
        3.01581553508235416007e-4,
        2.65806974686737550832e-6,
        6.23974539184983293730e-9,
    ]

    Q2 = [
        1,
        6.02427039364742014255e0,
        3.67983563856160859403e0,
        1.37702099489081330271e0,
        2.16236993594496635890e-1,
        1.34204006088543189037e-2,
        3.28014464682127739104e-4,
        2.89247864745380683936e-6,
        6.79019408009981274425e-9,
    ]

    if y0 <= 0 or y0 >= 1:
        raise ValueError("normalInverseCDF(x) needs 0 < x < 1")
    negate = True
    y = y0
    if y > 1.0 - 0.13533528323661269189:
        y = 1.0 - y
        negate = False

    if y > 0.13533528323661269189:
        y = y - 0.5
        y2 = y * y
        x = y + y * (y2 * polevl(y2, P0) / polevl(y2, Q0))
        x = x * s2pi
        return x

    x = sqrt(-2.0 * log(y))
    x0 = x - log(x) / x

    z = 1.0 / x
    if x < 8.0:
        x1 = z * polevl(z, P1) / polevl(z, Q1)
    else:
        x1 = z * polevl(z, P2) / polevl(z, Q2)
    x = x0 - x1
    if negate:
        x = -x
    return x


def doubleFactorial(n):
    """
    Returns double factorial of an integer.
    """
    return reduce(mul, range(n, 0, -2))


def binomial(n, k):
    """
    Returns binomial coefficient (n,k) or 'n choose k' with n ≥ k.
    """
    if n < k:
        raise ValueError("First argument must be greater than second.")
    return factorial(n) // (factorial(k) * factorial(n - k))


def flatten(sequence):
    """
    Returns a generator of all non-iterable elements in the input
    sequence. This can be used to flatten nested heterogeneous arrays.
    It may crash on deeply-nested array.

    Source: https://stackoverflow.com/a/2158532
    """
    for element in sequence:
        if isinstance(element, Iterable) and not isinstance(element, (str, bytes)):
            yield from flatten(element)
        else:
            yield element
