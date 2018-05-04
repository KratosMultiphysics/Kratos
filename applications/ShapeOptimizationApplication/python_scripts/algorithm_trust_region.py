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
import mapper_factory
import data_logger_factory
from custom_timer import Timer

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

            for id in self.objective_ids:
                self.communicator.requestValueOf(id)
                self.communicator.requestGradientOf(id)

            for id in self.equality_constraint_ids:
                self.communicator.requestValueOf(id)
                self.communicator.requestGradientOf(id)

            for id in self.inequality_constraint_ids:
                self.communicator.requestValueOf(id)
                self.communicator.requestGradientOf(id)

            self.analyzer.AnalyzeDesignAndReportToCommunicator( self.design_surface, self.optimizationIteration, self.communicator )

            # Get response from communicator and store values as list and gradients on nodes and as list of lists (matrix)
            obj_values = []
            obj_gradients = [[]]
            for id in self.objective_ids:
                obj_value = self.communicator.getStandardizedValue(id)
                obj_gradients_dict = self.communicator.getStandardizedGradient(id)

                nodal_variable = KratosGlobals.GetVariable("OBJECTIVE_SENSITIVITY")
                self.__StoreGradientDataOnNodalVariable(obj_gradients_dict,nodal_variable)

                mapped_nodal_variable = KratosGlobals.GetVariable("MAPPED_OBJECTIVE_SENSITIVITY")
                self.Mapper.MapToDesignSpace(nodal_variable, mapped_nodal_variable)
                self.Mapper.MapToGeometrySpace(mapped_nodal_variable, mapped_nodal_variable)

                obj_gradients_list = self.__ReadNodalVariableToList(mapped_nodal_variable)

                obj_values.append(obj_value)
                obj_gradients.append(obj_gradients_list)


















            eq_values = []
            eq_gradients = [[]]
            for id in self.equality_constraint_ids:
                eq_values.append(self.communicator.getStandardizedValue(id))
                eq_gradients_dict = self.communicator.getStandardizedGradient(id)
                eq_gradients_list = self.__ReadNodalVariableToList(eq_gradients_dict)
                eq_gradients.append(eq_gradients_list)

            ineq_values = []
            ineq_gradients = [[]]
            for id in self.inequality_constraint_ids:
                ineq_values.append(self.communicator.getStandardizedValue(id))
                ineq_gradients_dict = self.communicator.getStandardizedGradient(id)
                ineq_gradients_list = self.__ReadNodalVariableToList(ineq_gradients_dict)
                ineq_gradients.append(ineq_gradients_list)

            # Reset possible shape modifications during analysis
            self.model_part_controller.SetMeshToReferenceMesh()
            self.model_part_controller.SetDeformationVariablesToZero()

            # Convert to length direction format
            dir_obj,l_eq,dir_eq,l_ineq,dir_ineq = self.__ConvertToLengthDirectionFormat(obj_gradients,eq_values,eq_gradients,ineq_values,ineq_gradients)






            print("\n> Time needed for current optimization step = ", timer.GetLapTime(), "s")
            print("> Time needed for total optimization so far = ", timer.GetTotalTime(), "s")

    # --------------------------------------------------------------------------
    def FinalizeOptimizationLoop( self ):
        self.analyzer.FinalizeAfterOptimizationLoop()

    # --------------------------------------------------------------------------
    def __DetermineResponseIds(self, optimization_settings):
        objective_ids = []
        for itr in range(optimization_settings["objectives"].size()):
            objective_ids.append(optimization_settings["objectives"][itr]["identifier"].GetString())

        equality_constraint_ids = []
        inequality_constraint_ids = []
        for itr in range(optimization_settings["constraints"].size()):
            if(optimization_settings["constraints"][itr]["type"].GetString()=="="):
                equality_constraint_ids.append(optimization_settings["constraints"][itr]["identifier"].GetString())
            else:
                inequality_constraint_ids.append(optimization_settings["constraints"][itr]["identifier"].GetString())

        return objective_ids, equality_constraint_ids, inequality_constraint_ids

    # --------------------------------------------------------------------------
    def __StoreGradientDataOnNodalVariable(self, gradient_dict, nodal_variable):
        for nodeId, tmp_gradient in gradient_dict.items():
            self.optimization_model_part.Nodes[nodeId].SetSolutionStepValue(nodal_variable,0,tmp_gradient)

    # --------------------------------------------------------------------------
    def __ReadNodalVariableToList(self, nodal_variable):
        variable_values_list = []
        for node in self.design_surface.Nodes:
            tmp_values = node.GetSolutionStepValue(nodal_variable)
            variable_values_list.append(tmp_values[0])
            variable_values_list.append(tmp_values[1])
            variable_values_list.append(tmp_values[2])
        return variable_values_list

    # --------------------------------------------------------------------------
    def __ConvertToLengthDirectionFormat(self, obj_gradients, eq_values, eq_gradients, ineq_values, ineq_gradients):

        lInequality = [0 for _ in inequality]
        eInequality = [0 for _ in inequality]
        lEquality = [0 for _ in equality]
        eEquality = [0 for _ in equality]

        dir_obj = self.__ValueFormatToLengthFormat(obj_gradients)

        for i in range(self.ni):
            lInequality[i],eInequality[i] = self.__ValueFormatToLengthFormat(inequality[i],gradInequality[i])
        for i in range(self.ne):
            lEquality[i],eEquality[i] = self.__ValueFormatToLengthFormat(equality[i],gradEquality[i])

        return dir_obj,lInequality,eInequality,lEquality,eEquality

    # --------------------------------------------------------------------------
    def __ValueFormatToLengthFormat(self, gradient, value=None):
        gradientOriginal = deepcopy(gradient)
        if self.dampingIsSpecified:
            if self.dampOnlyAfterMapping:
                print("CAREFUL: DAMP ONLY ONCE")
            else:
                gradient = self.__DampVector(gradient)
                print("damp also before mapping")
        gradient = self.__MapVector(gradient)
        if self.dampingIsSpecified:
            gradient = self.__DampVector(gradient)
        ninf = norminf3d(gradient)
        eDir = [-gradient[i]/ninf for i in range(len(gradient))]
        gradL = sum(a*b for a,b in zip(gradientOriginal,eDir)) # dot product
        l = -value/gradL

        if value == None:
            return eDir
        else:
            return eDir,l

# ==============================================================================
