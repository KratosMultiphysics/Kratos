from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

import sys
from math import *

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignScalarVariableProcess(Model, settings["Parameters"])

class AssignScalarVariableProcess(KratosMultiphysics.Process):
    """This process sets a given scalar value for a certain variable in all the nodes of a submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings ):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

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

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Here i do a trick, since i want to allow "value" to be a string or a double value
        if settings.Has("value"):
            if settings["value"].IsString():
                default_settings["value"].SetString("0.0")

        settings.ValidateAndAssignDefaults(default_settings)

        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_name"].GetString())
        if not isinstance(self.variable, KratosMultiphysics.Array1DComponentVariable) and not isinstance(self.variable, KratosMultiphysics.DoubleVariable) and not isinstance(self.variable, KratosMultiphysics.VectorVariable):
            msg = "Error in AssignScalarToNodesProcess. Variable type of variable : " + settings["variable_name"].GetString() + " is incorrect . Must be a scalar or a component"
            raise Exception(msg)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.mesh = self.model_part.GetMesh(settings["mesh_id"].GetInt())
        self.is_fixed = settings["constrained"].GetBool()

        self.value_is_numeric = False
        if settings["value"].IsNumber():
            self.value_is_numeric = True
            self.value = settings["value"].GetDouble()
        else:
            self.function_string = settings["value"].GetString()
            self.aux_function = KratosMultiphysics.PythonGenericFunctionUtility(self.function_string, settings["local_axes"])

            if self.aux_function.DependsOnSpace():
                self.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(self.mesh.Nodes, self.aux_function )

        # Construct a variable_utils object to speedup fixing
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.step_is_active = False

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed in before initialize the solution step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):

            self.step_is_active = True

            if self.is_fixed:
                self.variable_utils.ApplyFixity(self.variable, self.is_fixed, self.mesh.Nodes)

            if self.value_is_numeric:
                self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
            else:
                if self.aux_function.DependsOnSpace() == False: #depends on time only
                    self.value = self.aux_function.CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                    self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
                else: #most general case - space varying function (possibly also time varying)
                    self.cpp_apply_function_utility.ApplyFunction(self.variable, current_time)

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        if self.step_is_active:
            # Here we free all of the nodes in the mesh
            if self.is_fixed:
                fixity_status  = False
                self.variable_utils.ApplyFixity(self.variable, fixity_status, self.mesh.Nodes)

        self.step_is_active = False
