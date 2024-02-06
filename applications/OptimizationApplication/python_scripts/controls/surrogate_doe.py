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
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA
import numpy as np
import math

def Factory(models: 'dict[Kratos.Model]', parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"MaterialPropertiesControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaterialPropertiesControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    if "aoa_free_stream_mach_control" == parameters["name"].GetString():
        return AngleOfAttackFreeStreamMachControl(parameters["name"].GetString(), models, parameters["settings"])
    elif "materials_control" == parameters["name"].GetString():
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
           "control_model_name" : "",
           "names_of_models_to_reset"    : [],
           "analysis_name"     : "",
           "variable_names_list"     : [],
           "property_ID"       : 0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.models = models
        self.control_model_name = parameters["control_model_name"].GetString()
        self.models_to_reset = parameters["names_of_models_to_reset"].GetStringArray()
        self.analysis_name = parameters["analysis_name"].GetString()
        self.control_vars = parameters["variable_names_list"].GetStringArray()
        self.property_ID = parameters["property_ID"].GetInt()
        self.variable_utils = Kratos.VariableUtils()
        self.scalar_vars = ["PRESSURE", "REACTION_WATER_PRESSURE"]
        self.vec_vars = ["VELOCITY", "MESH_DISPLACEMENT", "MESH_VELOCITY", "MESH_ACCELERATION", "MESH_REACTION", "ACCELERATION","DISPLACEMENT", "ROTATION", "REACTION", "FORCE", "POINT_LOAD"]

    def Initialize(self) -> None:
        pass

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

        for model_part in self.models[self.control_model_name].GetModelPartNames():
            mdp = self.models[self.control_model_name].GetModelPart(model_part)
            for var in self.control_vars:
                mdp.Properties[self.property_ID].SetValue( Kratos.KratosGlobals.GetVariable( var ), self.control_var_values_train[var][step] )
    
    def NextValidationSet(self, step) -> bool:
        
        self.ResetModels()

        for model_part in self.models[self.control_model_name].GetModelPartNames():
            mdp = self.models[self.control_model_name].GetModelPart(model_part)
            for var in self.control_vars:
                mdp.Properties[self.property_ID].SetValue( Kratos.KratosGlobals.GetVariable( var ), self.control_var_values_valid[var][step] )
        
    def GetControlVarNames(self):
        return self.control_vars
    
    def SetTrainDOEValues(self,values):
        self.control_var_values_train = values
    
    def SetValidDOEValues(self,values):
        self.control_var_values_valid = values
    
    def ResetModels(self):
        for mdl in self.models_to_reset:
            model = self.models[mdl]
            mdpNames = model.GetModelPartNames()
            mdp = model.GetModelPart(mdpNames[0])
            root_mdp = mdp.GetRootModelPart()

            self.variable_utils.UpdateCurrentToInitialConfiguration(root_mdp.GetNodes())

            for vec_field in self.vec_vars:
                 if root_mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(vec_field)):
                    self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(vec_field),root_mdp.GetNodes())

            for scala_field in self.scalar_vars:
                if root_mdp.HasNodalSolutionStepVariable(Kratos.KratosGlobals.GetVariable(scala_field)):
                    self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(scala_field),root_mdp.GetNodes())


    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"


class AngleOfAttackFreeStreamMachControl(Control):
    """Material properties control

    This is a generic material properties control which creates a control
    for the specified control variable. This does not do any filtering.

    TODO: Extend with filtering techniques when they are implemented.

    """
    def __init__(self, name: str, models: 'dict[Kratos.Model]', parameters: Kratos.Parameters):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
           "control_model_name" : "",
           "names_of_models_to_reset"    : [],
           "analysis_name"     : "",
           "variable_names_list"     : [],
           "angle_of_attack" : 0.0
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.models = models
        self.control_model_name = parameters["control_model_name"].GetString()
        self.models_to_reset = parameters["names_of_models_to_reset"].GetStringArray()
        self.analysis_name = parameters["analysis_name"].GetString()
        self.control_vars = parameters["variable_names_list"].GetStringArray()
        self.angle_of_attack = parameters["angle_of_attack"].GetDouble()
        self.variable_utils = Kratos.VariableUtils()
        self.scalar_vars = ["PRESSURE", "REACTION_WATER_PRESSURE", "VELOCITY_POTENTIAL", "AUXILIARY_VELOCITY_POTENTIAL",
                            "NEGATIVE_FACE_PRESSURE", "POSITIVE_FACE_PRESSURE"]
        self.vec_vars = ["VELOCITY", "MESH_DISPLACEMENT", "MESH_VELOCITY", "MESH_ACCELERATION", "MESH_REACTION", "ACCELERATION",
                         "VOLUME_ACCELERATION","DISPLACEMENT", "ROTATION", "REACTION", "FORCE", "POINT_LOAD", "LINE_LOAD", "SURFACE_LOAD",
                         "REACTION_MOMENT", "POINT_MOMENT"]

    def Initialize(self) -> None:
        model_parts = self.models[self.control_model_name].GetModelPartNames()
        mdp = self.models[self.control_model_name].GetModelPart(model_parts[0])
        root_model_part = mdp.GetRootModelPart()
        self.domain_size = root_model_part.ProcessInfo.GetValue(Kratos.DOMAIN_SIZE)

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

        model_parts = self.models[self.control_model_name].GetModelPartNames()
        mdp = self.models[self.control_model_name].GetModelPart(model_parts[0])
        root_model_part = mdp.GetRootModelPart()
        if  "FREE_STREAM_MACH" in self.control_vars:
            root_model_part.ProcessInfo.SetValue( Kratos.KratosGlobals.GetVariable( "FREE_STREAM_MACH" ), self.control_var_values_train["FREE_STREAM_MACH"][step] )
        if "ANGLE_OF_ATTACK" in self.control_vars:
            self.angle_of_attack = self.control_var_values_train["ANGLE_OF_ATTACK"][step]
        self.u_inf = root_model_part.ProcessInfo.GetValue(Kratos.KratosGlobals.GetVariable("FREE_STREAM_MACH")) * root_model_part.ProcessInfo.GetValue(Kratos.KratosGlobals.GetVariable("SOUND_VELOCITY"))
        self.free_stream_velocity = Kratos.Vector(3)
        if self.domain_size == 3:
            self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
            self.free_stream_velocity[1] = 0.0
            self.free_stream_velocity[2] = round(self.u_inf*math.sin(self.angle_of_attack),8)
        if self.domain_size == 2:
            self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
            self.free_stream_velocity[1] = round(self.u_inf*math.sin(self.angle_of_attack),8)
            self.free_stream_velocity[2] = 0.0
        self.free_stream_velocity_direction = self.free_stream_velocity / self.u_inf
        root_model_part.ProcessInfo.SetValue(Kratos.KratosGlobals.GetVariable("FREE_STREAM_VELOCITY"),self.free_stream_velocity)
        root_model_part.ProcessInfo.SetValue(Kratos.KratosGlobals.GetVariable("FREE_STREAM_VELOCITY_DIRECTION"),self.free_stream_velocity_direction)
    
    def NextValidationSet(self, step) -> bool:
        
        self.ResetModels()

        model_parts = self.models[self.control_model_name].GetModelPartNames()
        mdp = self.models[self.control_model_name].GetModelPart(model_parts[0])
        root_model_part = mdp.GetRootModelPart()
        if  "FREE_STREAM_MACH" in self.control_vars:
            root_model_part.ProcessInfo.SetValue( Kratos.KratosGlobals.GetVariable( "FREE_STREAM_MACH" ), self.control_var_values_valid["FREE_STREAM_MACH"][step] )
        if "ANGLE_OF_ATTACK" in self.control_vars:
            self.angle_of_attack = self.control_var_values_valid["ANGLE_OF_ATTACK"][step]
        self.u_inf = root_model_part.ProcessInfo.GetValue(Kratos.KratosGlobals.GetVariable("FREE_STREAM_MACH")) * root_model_part.ProcessInfo.GetValue(Kratos.KratosGlobals.GetVariable("SOUND_VELOCITY"))
        self.free_stream_velocity = Kratos.Vector(3)
        if self.domain_size == 3:
            self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
            self.free_stream_velocity[1] = 0.0
            self.free_stream_velocity[2] = round(self.u_inf*math.sin(self.angle_of_attack),8)
        if self.domain_size == 2:
            self.free_stream_velocity[0] = round(self.u_inf*math.cos(self.angle_of_attack),8)
            self.free_stream_velocity[1] = round(self.u_inf*math.sin(self.angle_of_attack),8)
            self.free_stream_velocity[2] = 0.0
        self.free_stream_velocity_direction = self.free_stream_velocity / self.u_inf
        root_model_part.ProcessInfo.SetValue(Kratos.KratosGlobals.GetVariable("FREE_STREAM_VELOCITY"),self.free_stream_velocity)
        root_model_part.ProcessInfo.SetValue(Kratos.KratosGlobals.GetVariable("FREE_STREAM_VELOCITY_DIRECTION"),self.free_stream_velocity_direction)
        
    def GetControlVarNames(self):
        return self.control_vars
    
    def SetTrainDOEValues(self,values):
        self.control_var_values_train = values
    
    def SetValidDOEValues(self,values):
        self.control_var_values_valid = values
    
    def ResetModels(self):
        for mdl in self.models_to_reset:
            model = self.models[mdl]
            mdpNames = model.GetModelPartNames()
            mdp = model.GetModelPart(mdpNames[0])
            root_mdp = mdp.GetRootModelPart()
            self.variable_utils.UpdateCurrentToInitialConfiguration(root_mdp.GetNodes())
            histVars = root_mdp.GetHistoricalVariablesNames()
            for vec_field in self.vec_vars:
                if vec_field in histVars:
                    self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(vec_field),root_mdp.GetNodes())
                else:
                    self.variable_utils.SetNonHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(vec_field),root_mdp.GetNodes())

            for scala_field in self.scalar_vars:
                if scala_field in histVars:
                    self.variable_utils.SetHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(scala_field),root_mdp.GetNodes())
                else:
                    self.variable_utils.SetNonHistoricalVariableToZero(Kratos.KratosGlobals.GetVariable(scala_field),root_mdp.GetNodes())
                


    def __str__(self) -> str:
        return f"Control [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}, control variable = {self.controlled_physical_variable.Name()}"
