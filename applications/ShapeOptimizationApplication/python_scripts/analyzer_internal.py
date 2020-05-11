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
import KratosMultiphysics.ShapeOptimizationApplication as KSO

from KratosMultiphysics.ShapeOptimizationApplication import model_part_controller_factory
import KratosMultiphysics.kratos_utilities as kratos_utilities

# Additional imports
from .analyzer_base import AnalyzerBaseClass
from .response_functions import response_function_factory as sho_response_factory
try:
    from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory as csm_response_factory
except ImportError:
    csm_response_factory = None
try:
    from KratosMultiphysics.ConvectionDiffusionApplication.response_functions import convection_diffusion_response_function_factory as convdiff_response_factory
except ImportError:
    convdiff_response_factory = None


import time as timer
import shutil
import glob, os

# ==============================================================================
class KratosInternalAnalyzer( AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def __init__( self, specified_responses, model_part_controller ):
        self.model_part_controller = model_part_controller
        self.specified_responses = specified_responses
        self.model = model_part_controller.GetModel()
        self.original_directory = os.getcwd()
  
    #---------------------------------------------------------------------------    
    def InitializeBeforeOptimizationLoop( self ):
        print("::InitOpti Ignored")
    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        
        optimization_model_part = self.model_part_controller.GetOptimizationModelPart()
        model_part_nodes = optimization_model_part.Nodes

        # Copying Opti ModelPart Nodes to Primal ModelPart
        x = []
        y = []
        z = []
        for node in model_part_nodes:
            x.append(node.X)
            y.append(node.Y)
            z.append(node.Z)
       
        time_before_analysis = optimization_model_part.ProcessInfo.GetValue(KM.TIME)
        step_before_analysis = optimization_model_part.ProcessInfo.GetValue(KM.STEP)
        delta_time_before_analysis = optimization_model_part.ProcessInfo.GetValue(KM.DELTA_TIME)

        if optimizationIteration == 1:

            self.response_functions = {} 

            for (response_id, response_settings) in self.specified_responses:
                # Clean OLD Response Folder
                kratos_utilities.DeleteDirectoryIfExisting('Response_'+response_id)

                self.response_functions[response_id] = csm_response_factory.CreateResponseFunction(response_id, response_settings, self.model) 

        else:    
            for identifier, response in self.response_functions.items():
                if identifier == "mass":
                    response.model.DeleteModelPart(response.model_part.Name)        # Other than Opti ITR 1, delete ModelPart
                    print("::ModelPart Deleted::", response.model_part.Name)
                else:
                    response.model.DeleteModelPart(response.primal_model_part.Name) # Other than Opti ITR 1, delete ModelPart
                    print("::ModelPart Deleted::", response.primal_model_part.Name)

            for (response_id, response_settings) in self.specified_responses:
                self.response_functions[response_id] = csm_response_factory.CreateResponseFunction(response_id, response_settings, self.model)


        
        for identifier, response in self.response_functions.items():    

            # Reset step/time iterators such that they match the optimization iteration after calling CalculateValue (which internally calls CloneTimeStep)
            optimization_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis-1)
            optimization_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0)

            response.SetCoordinatesUpdate(x, y, z)  #Transfer the Opti ITR Node Coordinates

            response.Initialize()

            ## For structural anlaysis of modified shape, here you can create .txt file of modified nodes at ith iteration
            
            # if optimizationIteration == 20:
            #     # Create Optimization Iteration Folder
            #     os.makedirs('Response_'+identifier +'/Opti_ITR_'+str(optimizationIteration))
            #     # Change Directory to ~/Response[i]/Opti_ITR_[i]/
            #     os.chdir('Response_'+identifier + '/Opti_ITR_'+str(optimizationIteration))
            #     file = open("Modified_XYZ.txt","w+")
            #     for node in model_part_nodes:
            #         file.write('\n\t' + str(node.Id) +'  '+ str(node.X) +' '+ str(node.Y) +' '+ str(node.Z))
            #     file.close()
            
            response.InitializeSolutionStep()

            # response values
            if communicator.isRequestingValueOf(identifier):
                response.CalculateValue()
                communicator.reportValue(identifier, response.GetValue())
            
            # response gradients
            if communicator.isRequestingGradientOf(identifier):
                response.CalculateGradient()
                communicator.reportGradient(identifier, response.GetNodalGradient(KM.SHAPE_SENSITIVITY))

            response.FinalizeSolutionStep()

            # Clear results or modifications on model part
            optimization_model_part.ProcessInfo.SetValue(KM.STEP, step_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(KM.TIME, time_before_analysis)
            optimization_model_part.ProcessInfo.SetValue(KM.DELTA_TIME, delta_time_before_analysis)  

            # Change to Original folder
            os.chdir(self.original_directory)               

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
            elif response_type in sho_response_functions:
                response_functions[response_id] = sho_response_factory.CreateResponseFunction(response_id, response_settings, model)
            else:
                raise NameError("The response function '{}' of type '{}' is not available.".format(response_id, response_type))

        return response_functions
