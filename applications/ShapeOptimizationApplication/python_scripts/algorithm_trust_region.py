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
class AlgorithmTrustRegion(OptimizationAlgorithm):
    # --------------------------------------------------------------------------
    def __init__(self, optimization_settings, analyzer, communicator, model_part_controller):
        default_algorithm_settings = Parameters("""
        {
            "name"            : "trust_region",
            "max_step_length" : 1,
            "max_iterations"  : 10
        }""")
        self.algorithm_settings =  optimization_settings["optimization_algorithm"]
        self.algorithm_settings.RecursivelyValidateAndAssignDefaults(default_algorithm_settings)

        self.analyzer = analyzer
        self.communicator = communicator
        self.model_part_controller = model_part_controller

        self.optimization_model_part = model_part_controller.GetOptimizationModelPart()
        self.design_surface = model_part_controller.GetDesignSurface()

        self.objective_ids, self.equality_constraint_ids, self.inequality_constraint_ids = self._DetermineResponseIds(optimization_settings)

        self.mapper = mapper_factory.CreateMapper(self.design_surface, optimization_settings["design_variables"]["filter"])

    # --------------------------------------------------------------------------
    def CheckApplicability(self):
        if len(self.objective_ids) > 1:
            raise RuntimeError("Trust-region algorithm only supports one objective function!")

    # --------------------------------------------------------------------------
    def InitializeOptimizationLoop(self):
        self.model_part_controller.ImportOptimizationModelPart()
        self.model_part_controller.InitializeMeshController()
        self.mapper.InitializeMapping()
        self.analyzer.InitializeBeforeOptimizationLoop()

    # --------------------------------------------------------------------------
    def RunOptimizationLoop(self):
        timer = Timer()
        timer.StartTimer()

        for self.optimizationIteration in range(1,self.algorithm_settings["max_iterations"].GetInt()):
            print("\n>===================================================================")
            print("> ",timer.GetTimeStamp(),": Starting optimization iteration ",self.optimizationIteration)
            print(">===================================================================\n")

            timer.StartNewLap()

            # Initialize new shape
            self.model_part_controller.UpdateMeshAccordingInputVariable(SHAPE_UPDATE)
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

            self.analyzer.AnalyzeDesignAndReportToCommunicator(self.design_surface, self.optimizationIteration, self.communicator)

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
            dir_objs = []
            for id, idx in self.objective_ids:
                nodal_variable = KratosGlobals.GetVariable("DF"+str(idx+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DF"+str(idx+1)+"DX_MAPPED")

                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

                obj_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                obj_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_obj = self.__ConvertToLengthDirectionFormat(obj_gradient, obj_gradient_mapped)
                dir_objs.append(dir_obj)

            len_eqs = []
            dir_eqs = []
            for id, idx in self.equality_constraint_ids:
                nodal_variable = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX_MAPPED")

                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

                eq_value = self.communicator.getStandardizedValue(id)
                eq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                eq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_eq, len_eq = self.__ConvertToLengthDirectionFormat(eq_gradient, eq_gradient_mapped, eq_value)

                dir_eqs.append(dir_eq)
                len_eqs.append(len_eq)

            len_ineqs = []
            dir_ineqs = []
            for id, idx in self.inequality_constraint_ids:
                nodal_variable = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX")
                nodal_variable_mapped = KratosGlobals.GetVariable("DC"+str(idx+1)+"DX_MAPPED")

                self.__PerformMapping(nodal_variable, nodal_variable_mapped)

                ineq_value = self.communicator.getStandardizedValue(id)
                ineq_gradient = ReadNodalVariableToList(self.design_surface, nodal_variable)
                ineq_gradient_mapped = ReadNodalVariableToList(self.design_surface, nodal_variable_mapped)

                dir_ineq, len_ineq = self.__ConvertToLengthDirectionFormat(ineq_gradient, ineq_gradient_mapped, ineq_value)

                dir_ineqs.append(dir_ineq)
                len_ineqs.append(len_ineq)

            # Determine step length
            max_step_length = self.algorithm_settings["max_step_length"].GetDouble()
            step_length = self.__ApplyStepLengthRule(max_step_length, objective_values[0], len_eqs, len_ineqs)

            len_bar_eqs = [l/step_length for l in len_eqs]
            len_bar_ineqs = [l/step_length for l in len_ineqs]





            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop(self):
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __PerformMapping(self,nodal_variable, nodal_variable_mapped):
        self.mapper.MapToDesignSpace(nodal_variable, nodal_variable_mapped)
        self.mapper.MapToGeometrySpace(nodal_variable_mapped, nodal_variable_mapped)

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
