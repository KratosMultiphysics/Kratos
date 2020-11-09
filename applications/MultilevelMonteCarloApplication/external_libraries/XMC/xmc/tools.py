import importlib
import math
import warnings

from xmc.distributedEnvironmentFramework import *

def dynamicImport(fullName):
    """
    This function returns the requested attribute of a module; it manages the necessary imports.
    The input argument is the full name as a string, e.g. 'module.package.method'.
    Alternatively, the input argument may be a list of such full names. The output will then be a list of same size containing the requested attributes in the same order.
    This only works for modern Python (at least Python > 3).
    """
    if fullName is None: #safety
        return None
    elif fullName is list or fullName is tuple:
        requested_attributes = []
        for name in fullName:
            requested_attributes.append(dynamicImport(name))
        return requested_attributes
    else:
        split_name = fullName.split('.')
        module_name = '.'.join(split_name[0:-1])
        attribute_name = split_name[-1]
        module = importlib.import_module(module_name)
        requested_attribute = getattr(module,attribute_name)
        return requested_attribute

def instantiateObject(fullName, **kwargs):
    object_type = dynamicImport(fullName)
    return object_type(**kwargs)

def doNothing(*args,**kwArgs):
    """
    This is a function that does nothing. Useful to hold a place without taking any action.
    """
    pass

def returnInput(*args):
    """
    This is a function that returns exactly what is passed to it. Useful to hold a place in the dataflow without modifying the data.
    """
    return args

@ExaquteTask(returns=1)
def returnInput_Task(*args):
    return returnInput(*args)

# For backward compatibility
convertObjectToFuture = returnInput_Task

#TODO Unused. remove?
def getUnionAndMap(listOfLists):
    """
    Given a list of lists, computes union of all lists, as well as map from list of lists
    to union lists.
    """
    union_list = []
    union_map = []
    for i in range(len(listOfLists)):
        union_map.append([None]*len(listOfLists[i]))
        for j in range(len(listOfLists[i])):
            found = False
            for k in range(len(union_list)):
                if listOfLists[i][j] == union_list[k]:
                    union_map[i][j] = k
                    found = True
                    break
            if found is False:
                union_list.append(listOfLists[i][j])
                union_map[i][j] = len(union_list)-1
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

def mergeTwoListsIntoOne(list1,list2):
    """
    Function that takes two lists and assembles them into one list
    such that each entry of the final list is a list containing the
    two corresponding entries of the two input lists
    """
    if (len(list1)!=len(list2)):
        raise ValueError('Input lists not of same length')

    output_list = [[list1[i],list2[i]] for i in range(len(list1))]
    return output_list

def strictlyPositiveBoundaryBooleans(indexSet):
    """
    Accepts a convex set of non-negative multi-indices and returns
    a list of booleans of the same size that states whether the
    corresponding index lies on the surface of the set where none
    of the indices have zero-values.
    """
    max_indices = [max([indexSet[i][j] for i in range(len(indexSet))]) for j in range(len(indexSet[0])) ]
    is_index_in_bias = [None for _ in indexSet]
    for i in range(len(indexSet)):
        index = indexSet[i]
        is_index_max = [index[i]==max_indices[i] for i in range(len(index))]
        is_index_nonzero = [index[i]!=0 for i in range(len(index))]
        is_index_in_bias[i] = any(is_index_max) and all(is_index_nonzero)
    return is_index_in_bias

#TODO is this function useful?
def summation(*args):
    return sum(args)

@ExaquteTask(returns=1)
def summation_Task(*args):
    return sum(args)

def packedList(obj):
    """
    Adapted from PyCOMPS module exaqute.ExaquteTask
    """

    warnings.warn(('packedList is deprecated. '
                   'Use COLLECTION type in task parameters with COMPSs version ≥ 2.6. '
                   'Retro-compatibility is ensured only until 2020-08.'),
                   DeprecationWarning)
    # Ensure that obj is not a tuple
    if isinstance(obj,tuple):
        obj = list(obj)
    index = [i for i, x in enumerate(obj) if x == "##"]
    if len(index) == 0:
        # The original nest list must have been of length 1,
        # e.g. [[1]], then obj=[1] and so we return [[1]]=[obj]
        return [obj]
    new_vector = []
    new_vector.append(list(obj[0:index[0]]))
    for i in range(len(index) - 1):
        new_vector.append(list(obj[(index[i] + 1):index[i + 1]]))
    new_vector.append(list(obj[(index[-1] + 1):len(obj)]))
    return new_vector

def unpackedList(obj):
    """
    Adapted from PyCOMPS module exaqute.ExaquteTask.
    """

    warnings.warn(('unpackedList is deprecated. '
                   'Use COLLECTION type in task parameters with COMPSs version ≥ 2.6. '
                   'Retro-compatibility is ensured only until 2020-08.'),
                   DeprecationWarning)
    if not any(isinstance(l,list) for l in obj):
        return obj
    new_vector = []
    new_vector.extend(obj[0])
    for i in range(1, len(obj)):
        new_vector.append("##")
        new_vector.extend(obj[i])
    return new_vector

def splitList(listToSplit, num_sublists=1):
    """
    Function that splits a list into a given number of sublists
    Reference: https://stackoverflow.com/questions/752308/split-list-into-smaller-lists-split-in-half
    """
    length = len(listToSplit)
    return [ listToSplit[i*length // num_sublists: (i+1)*length // num_sublists]
             for i in range(num_sublists) ]

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

    s2pi = 2.50662827463100050242E0

    P0 = [
        -5.99633501014107895267E1,
        9.80010754185999661536E1,
        -5.66762857469070293439E1,
        1.39312609387279679503E1,
        -1.23916583867381258016E0,
    ]

    Q0 = [
        1,
        1.95448858338141759834E0,
        4.67627912898881538453E0,
        8.63602421390890590575E1,
        -2.25462687854119370527E2,
        2.00260212380060660359E2,
        -8.20372256168333339912E1,
        1.59056225126211695515E1,
        -1.18331621121330003142E0,
    ]

    P1 = [
        4.05544892305962419923E0,
        3.15251094599893866154E1,
        5.71628192246421288162E1,
        4.40805073893200834700E1,
        1.46849561928858024014E1,
        2.18663306850790267539E0,
        -1.40256079171354495875E-1,
        -3.50424626827848203418E-2,
        -8.57456785154685413611E-4,
    ]

    Q1 = [
        1,
        1.57799883256466749731E1,
        4.53907635128879210584E1,
        4.13172038254672030440E1,
        1.50425385692907503408E1,
        2.50464946208309415979E0,
        -1.42182922854787788574E-1,
        -3.80806407691578277194E-2,
        -9.33259480895457427372E-4,
    ]

    P2 = [
        3.23774891776946035970E0,
        6.91522889068984211695E0,
        3.93881025292474443415E0,
        1.33303460815807542389E0,
        2.01485389549179081538E-1,
        1.23716634817820021358E-2,
        3.01581553508235416007E-4,
        2.65806974686737550832E-6,
        6.23974539184983293730E-9,
    ]

    Q2 = [
        1,
        6.02427039364742014255E0,
        3.67983563856160859403E0,
        1.37702099489081330271E0,
        2.16236993594496635890E-1,
        1.34204006088543189037E-2,
        3.28014464682127739104E-4,
        2.89247864745380683936E-6,
        6.79019408009981274425E-9,
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

    x = math.sqrt(-2.0 * math.log(y))
    x0 = x - math.log(x) / x

    z = 1.0 / x
    if x < 8.0:
        x1 = z * polevl(z, P1) / polevl(z, Q1)
    else:
        x1 = z * polevl(z, P2) / polevl(z, Q2)
    x = x0 - x1
    if negate:
        x = -x
    return x
