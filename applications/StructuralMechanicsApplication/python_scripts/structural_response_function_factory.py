from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics.StructuralMechanicsApplication import structural_response
from KratosMultiphysics.StructuralMechanicsApplication.structural_response import StrainEnergyResponseFunction
from KratosMultiphysics.StructuralMechanicsApplication.structural_response import MassResponseFunction
from KratosMultiphysics.StructuralMechanicsApplication.structural_response import Custom_StructuralMechanicsAnalysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_response import UpdateCoordinates_MassResponse
import KratosMultiphysics.ShapeOptimizationApplication as KSO

class Custom_Strain_Response(StrainEnergyResponseFunction):

    def __init__(self, response_id, response_settings, model):
        super(Custom_Strain_Response, self).__init__(response_id, response_settings, model)

    def SetCoordinatesUpdate(self, x, y, z):
        self.primal_analysis.SetCoordinatesUpdate(x, y, z, self.primal_model_part)

#============================================================================================================
class Custom_Mass_Response(MassResponseFunction):

    def __init__(self, response_id, response_settings, model):
        super(Custom_Mass_Response, self).__init__(response_id, response_settings, model)

    def SetCoordinatesUpdate(self, x, y, z):
        self.update.SetCoordinatesUpdate(x, y, z, self.model_part)


#============================================================================================================
class Custom_AdjointStrain_Response(structural_response.AdjointResponseFunction):

    def __init__(self, response_id, response_settings, model):
        super(Custom_AdjointStrain_Response, self).__init__(response_id, response_settings, model)
        self.model = model

    def SetCoordinatesUpdate(self, x, y, z):
        self.primal_analysis.SetCoordinatesUpdate(x, y, z, self.primal_model_part)


#============================================================================================================
def CreateResponseFunction(response_id, response_settings, model):
    response_type = response_settings["response_type"].GetString()

    if response_type == "strain_energy":
        return Custom_Strain_Response(response_id, response_settings, model)
        # return structural_response.StrainEnergyResponseFunction(response_id, response_settings, model)

    elif response_type == "mass":
        return Custom_Mass_Response(response_id, response_settings, model)
        # return structural_response.MassResponseFunction(response_id, response_settings, model)

    elif response_type == "eigenfrequency":
        return structural_response.EigenFrequencyResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_nodal_displacement":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_linear_strain_energy":
        return Custom_AdjointStrain_Response(response_id, response_settings, model)

    elif response_type == "adjoint_local_stress":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_max_stress":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    elif response_type == "adjoint_nodal_reaction":
        return structural_response.AdjointResponseFunction(response_id, response_settings, model)

    else:
        raise NameError("The type of the following response function is not specified: "+ response_id +
                        ".\nAvailable types are: 'mass', 'strain_energy', 'eigenfrequency', 'adjoint_nodal_displacement', 'adjoint_linear_strain_energy', 'adjoint_local_stress', 'adjoint_max_stress', 'adjoint_nodal_reaction'." )