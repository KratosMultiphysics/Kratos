def interpretationStructure():
    """
    This function returns the typical output of a MultiCriterion's interpreter,
    with default values. Such output is meant to be created by first calling
    this function, then changing specific values according the interpreter.
    The point is to always have the same data structure.
    """
    structure = {
        "stop": False,  # stop algorithm
        "updateHierarchy": True,
        "updateTolerance": True,  # for continuation-type algorithms
        "updateIndexSpace": False,
        "updateSampleNumberSpace": False,
    }
    return structure


def interpretAsStoppingFlag(flag):
    """
    This is the trivial interpreter for an input which is a boolean answering
    the question `Must I stop the algorithm now'?
    """
    interpretation = interpretationStructure()
    interpretation["stop"] = flag[0]
    return interpretation


def interpretAsConvergenceAndIterationBounds(flags):
    """
    This is the common interpreter which expects three booleans as input:
    1. Is convergence achieved?
    2. Is the number of iterations greater than the lower bound?
    3. Is the number of iterations greater than the uppder bound?
    """
    flag = flags[2] or (flags[0] and flags[1])
    interpretation = interpretAsStoppingFlag([flag])
    return interpretation


def interpretAsMultipleRequiredConvergencesAndIterationBounds(flags):
    """
    This is the common interpreter which expects N booleans as input:
    + 1 to N-2: convergence flags, all necessary (e.g. split tolerance)
    + N-1: is the number of iteration high enough?
    + N: is the number of iteration low enough?
    """
    flag = [all(flags[0:-2]), flags[-2], flags[-1]]
    interpretation = interpretAsConvergenceAndIterationBounds([flag])
    return interpretation


def interpretAsMultipleAlternativeConvergencesAndIterationBounds(flags):
    """
    This is the common interpreter which expects N booleans as input:
    + 1 to N-2: alternative convergence flags (at least one necessary)
    + N-1: is the number of iteration high enough?
    + N: is the number of iteration low enough?
    """
    flag = [any(flags[0:-2]), flags[-2], flags[-1]]
    interpretation = interpretAsConvergenceAndIterationBounds([flag])
    return interpretation
