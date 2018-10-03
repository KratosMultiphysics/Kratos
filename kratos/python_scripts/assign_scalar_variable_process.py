import KratosMultiphysics
import sys
from math import *


def Factory(settings, current_model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarVariableProcess(current_model, settings["Parameters"])

class AssignScalarVariableProcess(KratosMultiphysics.Process):
    def __init__(self, current_model, settings ):
        KratosMultiphysics.Process.__init__(self)

        #The value can be a double or a string (function)
        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"            : "This process sets a given scalar value for a certain variable in all the nodes of a submodelpart",
                "mesh_id"         : 0,
                "model_part_name" : "please_specify_model_part_name",
                "variable_name"   : "SPECIFY_VARIABLE_NAME",
                "interval"        : [0.0, 1e30],
                "constrained"     : true,
                "value"           : 0.0,
                "local_axes"      : {}
            }
            """
            )

        #assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        #here i do a trick, since i want to allow "value" to be a string or a double value
        if(settings.Has("value")):
            if(settings["value"].IsString()):
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if(type(self.variable) != KratosMultiphysics.Array1DComponentVariable and type(self.variable) != KratosMultiphysics.DoubleVariable and type(self.variable) != KratosMultiphysics.VectorVariable):
            msg = "Error in AssignScalarToNodesProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a scalar or a component"
            raise Exception(msg)

        self.model_part = current_model[settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.is_fixed = settings["constrained"].GetBool()

        self.value_is_numeric = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:
            self.function_string = settings["value"].GetString()
            self.aux_function = KratosMultiphysics.PythonGenericFunctionUtility(self.function_string, settings["local_axes"])

            if(self.aux_function.DependsOnSpace()):
                self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(self.mesh.Nodes, self.aux_function )

        #construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.step_is_active = False

    def ExecuteBeforeSolutionLoop(self):
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.interval.IsInInterval(current_time)):

            self.step_is_active = True

            if(self.is_fixed):
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)

            if self.value_is_numeric:
                self.variable_utils.SetScalarVar(self.variable, self.value, self.mesh.Nodes)
            else:
                if self.aux_function.DependsOnSpace() == False: #depends on time only
                    self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time)
                    self.variable_utils.SetScalarVar(self.variable, self.value, self.mesh.Nodes)
                else: #most general case - space varying function (possibly also time varying)
                    self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)

    def ExecuteFinalizeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if(self.step_is_active):
            #here we free all of the nodes in the mesh
            if(self.is_fixed):
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

        self.step_is_active = False
