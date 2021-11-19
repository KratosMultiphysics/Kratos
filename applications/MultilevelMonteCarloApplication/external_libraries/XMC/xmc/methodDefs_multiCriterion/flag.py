from exaqute import *


def plainFlag(values, criteria, valuesForCriterion, interpreter):
    """
    Returns the result of the criterion for the input information provided.
    This result is a dictionary whose entries match possible choices based on that criterion. See the definition of the interpreter for details.
    """
    # TODO do we keep this input argument as an array?
    elementaryFlags = []
    for i, criterion in enumerate(criteria):
        val = [values[j] for j in valuesForCriterion[i]]
        elementaryFlags.append(criterion.flag(*val))
    return interpreter(elementaryFlags)


# TODO Providing depth of collection improves performance
@task(keep=True, returns=1, values=COLLECTION_IN)
def plainFlag_Task(values, *args):
    """
    This method is identical to plainFlag but for the task decorator.
    """
    flag = plainFlag(values, *args)
    return flag
