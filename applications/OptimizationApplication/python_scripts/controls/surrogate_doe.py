from typing import Optional
from importlib import import_module
import KratosMultiphysics as Kratos
from KratosMultiphysics.analysis_stage import AnalysisStage
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import numpy as np
def Factory(models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"MaterialPropertiesControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaterialPropertiesControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return MaterialPropertiesControl(parameters["name"].GetString(), models, parameters["settings"])

class MaterialPropertiesControl(Control):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    TODO: Extend with filtering techniques when they are implemented.

    """
    def __init__(self, name: str, models: 'dict[Kratos.Model]', parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
           "model_name" : "",
           "root_model_part_name" : "",
           "controlVars" : {"analysis_name"     : "",
                            "variable_name"     : "",
                            "doeDomain"         : {}
                          },
            "outputDir"    : {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        parameters["controlVars"].ValidateAndAssignDefaults(default_settings["controlVars"])
        self.models = models
        self.model_name = parameters["model_name"].GetString()
        self.root_model_part_name = parameters["root_model_part_name"].GetString()
        self.analysis_name = parameters["controlVars"]["analysis_name"].GetString()
        self.output_dir = parameters["outputDir"]
        self.control_vars = parameters["controlVars"]
        self.variable_utils = Kratos.VariableUtils()
        self.scalar_vars = ["PRESSURE", "REACTION_WATER_PRESSURE"]
        self.vec_vars = ["VELOCITY", "MESH_DISPLACEMENT", "MESH_VELOCITY", "MESH_ACCELERATION", "MESH_REACTION", "ACCELERATION","DISPLACEMENT", "ROTATION", "REACTION", "FORCE", "POINT_LOAD"]

    def Initialize(self) -> None:
        pass
        #self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass
    
    
    def GetPhysicalKratosVariables(self):
        pass

    
    def GetEmptyField(self):
        pass

    
    def GetControlField(self):
        pass

    
    def MapGradient(self, physical_gradient_variable_container_expression_map):
        pass

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        pass

    def NextTrainingSet(self, step) -> bool:

        self.ResetModels()
        
        #fluid_model_part.ReduceTimeStep(fluid_model_part,0.0)                       ################# probably reduces just one step to particular time
        
        #for i in range(40):
        #    solid_model_part.OverwriteSolutionStepData(i+1,0)                         ################ leads to NaN values
        #    fluid_model_part.OverwriteSolutionStepData(i+1,0)

        for model_part in self.models[self.model_name].GetModelPartNames():
            mdp = self.models[self.model_name].GetModelPart(model_part)
            mdp.Properties[0].SetValue( Kratos.KratosGlobals.GetVariable( "YOUNG_MODULUS" ), self.control_var_values_train[step] )
        #for k,v in self.models.items():
        #    for model_part in v.GetModelPartNames():
        #       mdp = v.GetModelPart(model_part)
        #       self.variable_utils.UpdateCurrentToInitialConfiguration(mdp.GetNodes())
        #       for node in mdp.GetNodes():
        #           for vec_field in self.vec_vars:
        #                if mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(vec_field)):
        #                    node.SetSolutionStepValue(Kratos.KratosGlobals.GetVariable(vec_field), [0.0, 0.0, 0.0])
                   #for scala_field in self.scalar_vars:
                   #     if mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(scala_field)):
                   #         node.SetSolutionStepValue(Kratos.KratosGlobals.GetVariable(scala_field), 0)
               #mdp.ProcessInfo.ClearHistory(mdp.GetBufferSize())
               #for node in mdp.GetNodes():
               #    node.Clear()
        #        if mdp.Name != "AutomaticInlet2D_Inlet":
        #            self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("VELOCITY"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("MESH_DISPLACEMENT"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("MESH_VELOCITY"),mdp.GetNodes())
        #        #if mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable("PRESSURE")):
        #        #    for node in mdp.GetNodes():
        #        #        print(node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("PRESSURE"),40),"   ","for steo 40",'\n')
        #        #        print(node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("PRESSURE"),41),'\n')
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("PRESSURE"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("DISPLACEMENT"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("REACTION"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("FORCE"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("POINT_LOAD"),mdp.GetNodes())
        #        self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable("ROTATION"),mdp.GetNodes())

        #exec_policy.analysis: AnalysisStage = getattr(import_module(exec_policy.analysis_full_module), exec_policy.analysis_type)(exec_policy.analysis_settings.Clone(), exec_policy.models)
        #exec_policy.analysis.Initialize()
    
    def NextValidationSet(self, step) -> bool:
        
        self.ResetModels()
        
        #fluid_model_part.ReduceTimeStep(fluid_model_part,0.0)                       ################# probably reduces just one step to particular time
        
        #for i in range(40):
        #    solid_model_part.OverwriteSolutionStepData(i+1,0)                         ################ leads to NaN values
        #    fluid_model_part.OverwriteSolutionStepData(i+1,0)

        for model_part in self.models[self.model_name].GetModelPartNames():
            mdp = self.models[self.model_name].GetModelPart(model_part)
            mdp.Properties[0].SetValue( Kratos.KratosGlobals.GetVariable( "YOUNG_MODULUS" ), self.control_var_values_valid[step] )
        
    def GetControlVarName(self):
        return self.control_vars["variable_name"].GetString()
    
    def SetTrainDOEValues(self,values):
        self.control_var_values_train = values
    
    def SetValidDOEValues(self,values):
        self.control_var_values_valid = values
    
    def ResetModels(self):
        for name,model in self.models.items():
            mdpNames = model.GetModelPartNames()
            mdp = model.GetModelPart(mdpNames[0])
            root_mdp = mdp.GetRootModelPart()
            #model_part_names = model.GetModelPartNames()
            #for model_part in model_part_names:
            #    mdp = model.GetModelPart(model_part)
            #    mdp.ProcessInfo.SetValue(Kratos.IS_RESTARTED, False)
            #    mdp.ProcessInfo.SetValue(Kratos.STEP, 0)
            #    mdp.ProcessInfo.SetValue(Kratos.TIME, 0.0)
            #    mdp.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.0)
            #    if name == 'fluid':
            #        mdp.ProcessInfo.SetValue(Kratos.DELTA_TIME, 0.1)
            #step = root_mdp.ProcessInfo.GetValue(Kratos.STEP)
            #if step == 0:
            #    continue
            #solid_model_part = self.models["structure"].GetModelPart("Structure")
            #self.variable_utils.UpdateCurrentToInitialConfiguration(root_mdp.GetNodes())
            #for node in root_mdp.GetNodes():
            #    for vec_field in self.vec_vars:
            #         if root_mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(vec_field)):
            #             node.SetSolutionStepValue(Kratos.KratosGlobals.GetVariable(vec_field), [0.0, 0.0, 0.0])
            #solid_model_part.ReduceTimeStep(solid_model_part,0.0)
            #fluid_model_part = self.models["fluid"].GetModelPart("FluidModelPart")
            self.variable_utils.UpdateCurrentToInitialConfiguration(root_mdp.GetNodes())
            #for node in root_mdp.GetNodes():
                #if name == 'fluid':
                #    vel = [node.GetSolutionStepValue(Kratos.VELOCITY)[0], node.GetSolutionStepValue(Kratos.VELOCITY)[1], node.GetSolutionStepValue(Kratos.VELOCITY)[2]]
                #    acc = [node.GetSolutionStepValue(Kratos.ACCELERATION)[0], node.GetSolutionStepValue(Kratos.ACCELERATION)[1], node.GetSolutionStepValue(Kratos.ACCELERATION)[2]]
                #    meshVel = [node.GetSolutionStepValue(Kratos.MESH_VELOCITY)[0], node.GetSolutionStepValue(Kratos.MESH_VELOCITY)[1], node.GetSolutionStepValue(Kratos.MESH_VELOCITY)[2]]
                #    press = node.GetSolutionStepValue(Kratos.PRESSURE)
                #    isStruct = node.GetSolutionStepValue(Kratos.IS_STRUCTURE)
                #    disp = [node.GetSolutionStepValue(Kratos.DISPLACEMENT)[0], node.GetSolutionStepValue(Kratos.DISPLACEMENT)[1], node.GetSolutionStepValue(Kratos.DISPLACEMENT)[2]]
                #    bodyForce = [node.GetSolutionStepValue(Kratos.BODY_FORCE)[0], node.GetSolutionStepValue(Kratos.BODY_FORCE)[1], node.GetSolutionStepValue(Kratos.BODY_FORCE)[2]]
                #    nodal_area = node.GetSolutionStepValue(Kratos.NODAL_AREA)
                #    nodal_h = node.GetSolutionStepValue(Kratos.NODAL_H)
                #    advproj = [node.GetSolutionStepValue(Kratos.ADVPROJ)[0], node.GetSolutionStepValue(Kratos.ADVPROJ)[1], node.GetSolutionStepValue(Kratos.ADVPROJ)[2]]
                #    divproj = node.GetSolutionStepValue(Kratos.DIVPROJ)
                #    reaction = [node.GetSolutionStepValue(Kratos.REACTION)[0], node.GetSolutionStepValue(Kratos.REACTION)[1], node.GetSolutionStepValue(Kratos.REACTION)[2]]
                #    reacWater = node.GetSolutionStepValue(Kratos.REACTION_WATER_PRESSURE)
                #    extPress = node.GetSolutionStepValue(Kratos.EXTERNAL_PRESSURE)
                #    normal = [node.GetSolutionStepValue(Kratos.NORMAL)[0], node.GetSolutionStepValue(Kratos.NORMAL)[1], node.GetSolutionStepValue(Kratos.NORMAL)[2]]
                #    ywall = node.GetSolutionStepValue(Kratos.Y_WALL)
                #    qvalue = node.GetSolutionStepValue(KratosCFD.Q_VALUE)
            for vec_field in self.vec_vars:
                 if root_mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(vec_field)):
                    #ValueX = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("{}_X".format(vec_field)),0)
                    #ValueY = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("{}_Y".format(vec_field)),0)
                    #ValueZ = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("{}_Z".format(vec_field)),0)
                    #for i in range(0,20):
                    #    node.SetSolutionStepValue(Kratos.KratosGlobals.GetVariable(vec_field), i, [0.0, 0.0, 0.0])
                    #ValueX = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("{}_X".format(vec_field)),0)
                    #ValueY = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("{}_Y".format(vec_field)),0)
                    #ValueZ = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable("{}_Z".format(vec_field)),0)
                    self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(vec_field),root_mdp.GetNodes())
            for scala_field in self.scalar_vars:
                if root_mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(scala_field)):
                    #Value = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable(scala_field),2)
                    #for i in range(0,20):
                    #    node.SetSolutionStepValue(Kratos.KratosGlobals.GetVariable(scala_field), i, 0.0)
                    #Value = node.GetSolutionStepValue(Kratos.KratosGlobals.GetVariable(scala_field),2)
                    self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(scala_field),root_mdp.GetNodes())


    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"
