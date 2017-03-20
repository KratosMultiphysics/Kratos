import KratosMultiphysics
import assign_vector_components_to_nodes_process
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignTimeDerivativeProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class AssignTimeDerivativeProcess(KratosMultiphysics.Process):
    '''this process fixes the components of the value named "variable_to_be_solved_for" 
    in case the components of "variable_name" is fixed.
    this is needed when the variable to  be solved for is different than the one being fixed'''
    
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        if(settings.Has("variable_to_be_solved_for") == False):
            raise Exception("settings must contain an indication of the variable to be solved for")
        
        #create a process to fix the variable with name "variable_name" 
        aux_settings = KratosMultiphysics.Parameters("""{}"""
            )
        for key, value in settings.items():
            if(key != "variable_to_be_solved_for"):
                aux_settings.AddValue(key, value)
        self.velocity_application_process = assign_vector_components_to_nodes_process.AssignVectorComponentsToNodesProcess(Model,aux_settings)
        
        #here prepare the data to fix the variable "to be solved for" - note that "aux_settings" is used so that the checking on the defaults is done on the base class.
        self.model_part = Model[aux_settings["model_part_name"].GetString()]
        self.variable_utils = KratosMultiphysics.VariableUtils()
        self.interval = KratosMultiphysics.IntervalUtility(aux_settings)
        
        self.x_component = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_to_be_solved_for"].GetString()+"_X")
        self.y_component = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_to_be_solved_for"].GetString()+"_Y")
        self.z_component = KratosMultiphysics.KratosGlobals.GetVariable(settings["variable_to_be_solved_for"].GetString()+"_Z")
        
        self.fix_x = aux_settings["constrained"][0].GetBool()
        self.fix_y = aux_settings["constrained"][1].GetBool()
        self.fix_z = aux_settings["constrained"][2].GetBool()
        if(aux_settings["value"][0].IsNull()):
            self.fix_x = False
        if(aux_settings["value"][1].IsNull()):
            self.fix_y = False
        if(aux_settings["value"][2].IsNull()):
            self.fix_z = False
            
        self.interval_is_active = False

    def ExecuteInitializeSolutionStep(self):
        self.velocity_application_process.ExecuteInitializeSolutionStep()
        
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        
        self.interval_is_active = self.interval.IsInInterval(current_time)
        
        if(self.interval_is_active):
            fixed = True
            if(self.fix_x):
                self.variable_utils.ApplyFixity(self.x_component, fixed, self.model_part.Nodes)
            if(self.fix_y):
                self.variable_utils.ApplyFixity(self.y_component, fixed, self.model_part.Nodes)
            if(self.fix_z):
                self.variable_utils.ApplyFixity(self.z_component, fixed, self.model_part.Nodes)

    def ExecuteFinalizeSolutionStep(self):
        self.velocity_application_process.ExecuteFinalizeSolutionStep()

        if(self.interval_is_active):
            fixed = False #here we must free the velocity
            if(self.fix_x):
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, fixed, self.model_part.Nodes)
            if(self.fix_y):
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, fixed, self.model_part.Nodes)
            if(self.fix_z):
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, fixed, self.model_part.Nodes)
