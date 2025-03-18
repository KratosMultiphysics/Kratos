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

import os, pathlib, sys

# Kratos Core and Apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication.analyzers.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication.response_functions import response_function_factory as sho_response_factory
try:
    from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory as csm_response_factory
    import KratosMultiphysics.StructuralMechanicsApplication as KSM
except ImportError:
    csm_response_factory = None
    KSM = None
try:
    from KratosMultiphysics.ConvectionDiffusionApplication.response_functions import convection_diffusion_response_function_factory as convdiff_response_factory
except ImportError:
    convdiff_response_factory = None
try:
    from KratosMultiphysics.CompressiblePotentialFlowApplication import potential_flow_response_function_factory as potential_flow_response_factory
except ImportError:
    potential_flow_response_factory = None

import time as timer

class IterationScope:
    def __init__(self, response_id, iteration_number, is_evaluated_in_folder):
        self.is_evaluated_in_folder = is_evaluated_in_folder
        if (self.is_evaluated_in_folder):
            self.currentPath = pathlib.Path.cwd()
            output_path = pathlib.Path("Design_Iterations")
            response_text = "{:}/{:d}".format(response_id, iteration_number)
            self.scope = output_path / response_text

    def __enter__(self):
        if (self.is_evaluated_in_folder):
            self.scope.mkdir(parents=True, exist_ok=True)
            sys.path.insert(0, str(self.scope.absolute()))
            os.chdir(str(self.scope))

    def __exit__(self, exc_type, exc_value, traceback):
        if (self.is_evaluated_in_folder):
            os.chdir(self.currentPath)
            sys.path.remove(str(self.scope.absolute()))


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

            # now we scope in to the directory where response operations are done
            with IterationScope(identifier, optimizationIteration, response.IsEvaluatedInFolder()):
                response.UpdateDesign(optimization_model_part, KM.SHAPE_SENSITIVITY)

                response.InitializeSolutionStep()

                # response values
                if communicator.isRequestingValueOf(identifier):
                    response.CalculateValue()
                    communicator.reportValue(identifier, response.GetValue())

                # response gradients
                if communicator.isRequestingGradientOf(identifier) or communicator.isRequestingThicknessGradientOf(identifier):
                    response.CalculateGradient()

                # response shape gradients
                if communicator.isRequestingGradientOf(identifier):
                    communicator.reportGradient(identifier, response.GetNodalGradient(KM.SHAPE_SENSITIVITY))

                # response thickness gradients
                if communicator.isRequestingThicknessGradientOf(identifier):
                    if KSM is None:
                        raise RuntimeError(f"ThicknessOpt: {identifier} response function requires StructuralMechanicsApplication.")
                    communicator.reportThicknessGradient(identifier, response.GetElementalGradient(KSM.THICKNESS_SENSITIVITY))

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

        sho_response_functions = [
            "plane_based_packaging",
            "mesh_based_packaging",
            "surface_normal_shape_change",
            "face_angle",
            "water_drain",
            "directional_derivative",
            "airfoil_angle_of_attack",
            "airfoil_chord_length",
            "airfoil_perimeter"
        ]
        csm_response_functions = ["strain_energy", "mass", "eigenfrequency", "adjoint_local_stress", "adjoint_max_stress", "adjoint_nodal_displacement", "adjoint_linear_strain_energy", "adjoint_nodal_reaction"]
        cps_response_functions = ["adjoint_lift_potential_jump", "stochastic_adjoint_lift_potential_jump"]
        convdiff_response_functions = ["point_temperature"]

        for (response_id, response_settings) in specified_responses:
            if response_id in response_functions.keys():
                raise NameError("There are multiple response functions with the following identifier: " + response_id)

            response_type = response_settings["response_type"].GetString()

            if response_type in csm_response_functions:
                if csm_response_factory is None:
                    raise RuntimeError("ShapeOpt: {} response function requires StructuralMechanicsApplication.".format(response_type))
                response_functions[response_id] = csm_response_factory.CreateResponseFunction(response_id, response_settings, model)
            elif response_type in convdiff_response_functions:
                if convdiff_response_factory is None:
                    raise RuntimeError("ShapeOpt: {} response function requires ConvectionDiffusionApplication.".format(response_type))
                response_functions[response_id] = convdiff_response_factory.CreateResponseFunction(response_id, response_settings, model)
            elif response_type in cps_response_functions:
                if potential_flow_response_factory is None:
                    raise RuntimeError("ShapeOpt: {} response function requires CompressiblePotentialFlowApplication.".format(response_type))
                response_functions[response_id] = potential_flow_response_factory.CreateResponseFunction(response_id, response_settings, model)
            elif response_type in sho_response_functions:
                response_functions[response_id] = sho_response_factory.CreateResponseFunction(response_id, response_settings, model)
            else:
                raise NameError("The response function '{}' of type '{}' is not available.".format(response_id, response_type))

        return response_functions
