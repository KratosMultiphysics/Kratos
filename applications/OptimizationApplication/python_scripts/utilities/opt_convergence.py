import KratosMultiphysics as Kratos

def CreateConvergenceCriteria(parameters: Kratos.Parameters):
    type = parameters["type"].GetString()
    if type == "max_iter":
        return MaxIterConvCriterium(parameters)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class MaxIterConvCriterium(object):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "max_iter",
            "max_iter"          : 0,
        }""")

    def __init__(self, parameters: Kratos.Parameters):
        self.__max_iter = parameters["max_iter"].GetInt()

    def CheckConvergence(self, iter):
        conv = True if iter >= self.__max_iter else False
        msg = f"""\t Convergence info: 
            type          : {"max_iter"} 
            value         : {iter} of {self.__max_iter}
            status        : {conv}"""
        print(msg)
        return conv