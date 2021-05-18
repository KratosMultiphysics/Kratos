import KratosMultiphysics
from KratosMultiphysics import assign_vector_variable_process
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignTimeDerivativeProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AssignTimeDerivativeProcess(assign_vector_variable_process.AssignVectorVariableProcess):
    '''this process fixes the components of the value named "variable_to_be_solved_for"
    in case the components of "variable_name" is fixed.
    this is needed when the variable to  be solved for is different than the one being fixed'''

    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        if(settings.Has("variable_to_be_solved_for") == False):
            raise Exception("settings must contain an indication of the variable to be solved for")

        self.dof_components = []
        if(not settings["value"][0].IsNull()):
            self.dof_components.append( KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_to_be_solved_for"].GetString()+"_X") )
        if(not settings["value"][1].IsNull()):
            self.dof_components.append( KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_to_be_solved_for"].GetString()+"_Y") )
        if(not settings["value"][2].IsNull()):
            self.dof_components.append( KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_to_be_solved_for"].GetString()+"_Z") )

        settings.RemoveValue("variable_to_be_solved_for") #remove this value from the settings in order to be able to use the settings in constructing the base object

        #construct the base class
        super(AssignTimeDerivativeProcess, self).__init__(Model,settings)


    def ExecuteInitializeSolutionStep(self):
        super(AssignTimeDerivativeProcess, self).ExecuteInitializeSolutionStep()

        fixed = True
        for component in range(len(self.dof_components)):
            step_is_active =self.aux_processes[component].step_is_active
            is_fixed = self.aux_processes[component].step_is_active
            if(step_is_active == True and is_fixed==True):
                self.aux_processes[component].variable_utils.ApplyFixity(self.dof_components[component], fixed, self.model_part.Nodes)


    def ExecuteFinalizeSolutionStep(self):
        fixed = False
        for component in range(len(self.dof_components)):
            step_is_active =self.aux_processes[component].step_is_active
            is_fixed = self.aux_processes[component].step_is_active
            if(step_is_active == True and is_fixed==True):
                self.aux_processes[component].variable_utils.ApplyFixity(self.dof_components[component], fixed, self.model_part.Nodes)

        super(AssignTimeDerivativeProcess, self).ExecuteFinalizeSolutionStep()

