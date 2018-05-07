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
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from algorithm_base import OptimizationAlgorithm
from custom_math import NormInf3D, DotProduct
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToList
from custom_timer import Timer
import mapper_factory

# ==============================================================================
class AlgorithmTrustRegion( OptimizationAlgorithm ) :
    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings, model_part_controller, analyzer, communicator ):
        self.model_part_controller = model_part_controller
        self.analyzer = analyzer
        self.communicator = communicator

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()
        self.algorithm_settings = optimization_settings["optimization_algorithm"]

        self.objective_ids, self.equality_constraint_ids, self.inequality_constraint_ids = self.__DetermineResponseIds( optimization_settings )

        self.Mapper = mapper_factory.CreateMapper( model_part_controller, optimization_settings )

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop( self ):
        self.analyzer.InitializeBeforeOptimizationLoop()
        self.model_part_controller.InitializeMeshController()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop( self ):
        timer = Timer()
        timer.StartTimer()

        for self.optimizationIteration in range(1,self.algorithm_settings["max_iterations"].GetInt()):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.optimizationIteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            # Initialize new shape
            self.model_part_controller.UpdateMeshAccordingInputVariable( SHAPE_UPDATE )
            self.model_part_controller.SetReferenceMeshToMesh()

            # Analyze shape
            self.communicator.initializeCommunication()

            for id, idx in self.objective_ids:
                self.communicator.requestValueOf(id)
                self.communicator.requestGradientOf(id)

            for id, idx in self.equality_constraint_ids:
                self.communicator.requestValueOf(id)
                self.communicator.requestGradientOf(id)

            for id, idx in self.inequality_constraint_ids:
                self.communicator.requestValueOf(id)
                self.communicator.requestGradientOf(id)

            self.analyzer.AnalyzeDesignAndReportToCommunicator( self.design_surface, self.optimizationIteration, self.communicator )

            # Store values and gradients from analysis
            obj_values = []
            for id, idx in self.objective_ids:
                obj_values.append(self.communicator.getStandardizedValue(id))
                obj_gradients_dict = self.communicator.getStandardizedGradient(id)

                nodal_variable = KratosGlobals.GetVariable("DF"+str(idx+1)+"DX")
                WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

            eq_values = []
            for id, idx in self.equality_constraint_ids:
                eq_values.append(self.communicator.getStandardizedValue(id))
                eq_gradients_dict = self.communicator.getStandardizedGradient(id)

                nodal_variable = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX")
                WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

            ineq_values = []
            for id, idx in self.inequality_constraint_ids:
                ineq_values.append(self.communicator.getStandardizedValue(id))
                ineq_gradients_dict = self.communicator.getStandardizedGradient(id)

                nodal_variable = KratosGlobals.GetVariable("DC"+str(idx)+"DX")
                WriteDictionaryDataOnNodalVariable(obj_gradients_dict, self.optimization_model_part, nodal_variable)

            # Reset possible shape modifications during analysis
            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

            # Convert anylsis results to length direction format considering mapping and damping
            dir_obj = []
            for id, idx in self.objective_ids:
                nodal_variable = KratosGlobals.GetVariable("DF"+str(idx+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DF"+str(idx+1)+"DX_MAPPED")

                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

                obj_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                obj_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_obj_i = self.__ConvertToLengthDirectionFormat(obj_gradient, obj_gradient_mapped)
                dir_obj.append(dir_obj_i)

            len_eq = []
            dir_eq = []
            for id, idx in self.equality_constraint_ids:
                nodal_variable = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX_MAPPED")

                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

                eq_value = self.communicator.getStandardizedValue(id)
                eq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                eq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_eq_i, len_eq_i = self.__ConvertToLengthDirectionFormat(eq_gradient, eq_gradient_mapped, eq_value)

                dir_eq.append(dir_eq_i)
                len_eq.append(len_eq_i)

            len_ineq = []
            dir_ineq = []
            for id, idx in self.inequality_constraint_ids:
                nodal_variable = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX_MAPPED")

                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

                ineq_value = self.communicator.getStandardizedValue(id)
                ineq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                ineq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_ineq_i, len_ineq_i = self.__ConvertToLengthDirectionFormat(ineq_gradient, ineq_gradient_mapped, ineq_value)

                dir_ineq.append(dir_ineq_i)
                len_ineq.append(len_ineq_i)






            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __DetermineResponseIds(self, optimization_settings):
        objective_ids = []
        for idx in range(optimization_settings["objectives"].size()):
            objective_ids.append([optimization_settings["objectives"][idx]["identifier"].GetString(), idx])

        equality_constraint_ids = []
        inequality_constraint_ids = []
        for idx in range(optimization_settings["constraints"].size()):
            if(optimization_settings["constraints"][idx]["type"].GetString()=="="):
                equality_constraint_ids.append([optimization_settings["constraints"][idx]["identifier"].GetString(),idx])
            else:
                inequality_constraint_ids.append([optimization_settings["constraints"][idx]["identifier"].GetString(), idx])

        return objective_ids, equality_constraint_ids, inequality_constraint_ids

    # --------------------------------------------------------------------------
    def __PerformMapping(self,nodal_variable, nodal_variable_mapped):
        self.Mapper.MapToDesignSpace(nodal_variable, nodal_variable_mapped)
        self.Mapper.MapToGeometrySpace(nodal_variable_mapped, nodal_variable_mapped)

    # --------------------------------------------------------------------------
    def __ConvertToLengthDirectionFormat(self, gradient, modified_gradient, value=None):
        norm_inf = NormInf3D(modified_gradient)
        direction = [-modified_gradient[itr]/norm_inf for itr in range(len(modified_gradient))]

        if value == None:
            return direction
        else:
            grad_dot_dir = DotProduct(gradient, direction)
            length = -value/grad_dot_dir
            return direction, length

# ==============================================================================
