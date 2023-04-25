import KratosMultiphysics as Kratos
import standardized_objective as objective
import standardized_constraint as constraint

class CommonResponseFunction(object):

    def __init__(self) -> None:
        self.name
        self.response_type
        self.function = objective # or constraint, depend on type
        self.model

    def ComputeResponseValue(des_var, MC):
        MC.ComputePrimal(des_var, self.model)
        value = self.function.ComputeResponseValue() # we don't need to parse the des_var here anymore

    def ComputeResponseGradients(des_var, gradients, MC):
        MC.ComputePrimal(des_var)
        MC.ComputeAdjoint(des_var)


    