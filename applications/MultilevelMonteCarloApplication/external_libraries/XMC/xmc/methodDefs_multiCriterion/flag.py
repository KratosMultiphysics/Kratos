# Import PyCOMPSs
# from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
from exaqute.ExaquteTaskLocal import *      # to execute with python3

def plainFlag(self,values):
    """
    Returns the result of the criterion for the input information provided.
    This result is a dictionary whose entries match possible choices based on that criterion. See the definition of the interpreter for details.
    """
    #TODO do we keep this input argument as an array?
    elementary_flags = []
    for i in range(len(self.criteria)):
        current_values = []
        for j in range(len(self.inputsForCriterion[i])):
            current_values.append(values[self.inputsForCriterion[i][j]])
        new_elementary_flag = self.criteria[i].flag(*current_values)
        elementary_flags.append(new_elementary_flag)
    return self._interpreter(elementary_flags)

@ExaquteTask(returns=1)
def plainFlag_Task(*args):
    """
    This method is identical to plainFlag but for the task decorator.
    """
    return plainFlag(*args)
