# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
import KratosMultiphysics as KM

# Additional imports
from .analyzer_base import AnalyzerBaseClass
from .response_functions import response_function_factory as sho_response_factory
try:
    from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory as csm_response_factory
except ImportError:
    csm_response_factory = None
import time as timer

# ==============================================================================
class KratosInternalAnalyzer( AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, specified_responses, model_part_controller ):
        self.model_part_controller = model_part_controller
        self.response_functions = self.__CreateResponseFunctions(specified_responses, model_part_controller.GetModel())

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop( self ):
        for response in self.response_functions.values():
            response.Initialize()

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        optimization_model_part = self.model_part_controller.GetOptimizationModelPart()

        time_before_analysis = optimization_model_part.ProcessInfo.GetValue(KM.TIME)
        step_before_analysis = optimization_model_part.ProcessInfo.GetValue(KM.STEP)
        delta_time_before_analysis = optimization_model_part.ProcessInfo.GetValue(KM.DELTA_TIME)

        for identifier, response in self.response_functions.items():

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            optimization_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0)

            response.InitializeSolutionStep()

            # response values
            if communicator.isRequestingValueOf(identifier):
                response.CalculateValue()
                communicator.reportValue(identifier, response.GetValue())

            # response gradients
            if communicator.isRequestingGradientOf(identifier):
                response.CalculateGradient()
                communicator.reportGradient(identifier, response.GetShapeGradient())

            response.FinalizeSolutionStep()

            # Clear results or modifications on model part
            optimization_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_analysis)

            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop( self ):
        for response in self.response_functions.values():
            response.Finalize()

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateResponseFunctions( specified_responses, model ):
        response_functions = {}

        sho_response_functions = ["plane_based_packaging", "mesh_based_packaging"]
        csm_response_functions = ["strain_energy", "mass", "eigenfrequency", "adjoint_local_stress", "adjoint_max_stress"]

        for (response_id, response_settings) in specified_responses:
            if response_id in response_functions.keys():
                raise NameError("There are multiple response functions with the following identifier: " + response_id)

            response_type = response_settings["response_type"].GetString()

            if response_type in csm_response_functions:
                if csm_response_factory is None:
                    raise RuntimeError("ShapeOpt: {} response function requires StructuralMechanicsApplication.".format(response_type))
                response_functions[response_id] = csm_response_factory.CreateResponseFunction(response_id, response_settings, model)
            elif response_type in sho_response_functions:
                response_functions[response_id] = sho_response_factory.CreateResponseFunction(response_id, response_settings, model)
            else:
                raise NameError("The response function '{}' of type '{}' is not available.".format(response_id, response_type))

        return response_functions