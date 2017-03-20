import KratosMultiphysics
import assign_vector_components_to_nodes_process
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignAccelerationForDisplacementDrivenSolutionsProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class AssignAccelerationForDisplacementDrivenSolutionsProcess(KratosMultiphysics.Process):
    '''this process fixes the DISPLACEMENT components in case the VELOCITY components are fixed.
    this is needed when the variable to  be solved for is DISPLACEMENT and the user wants to prescribe VELOCITY'''
    
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        
        settings.AddEmptyValue("variable_name").SetString("ACCELERATION")
        self.variable_utils = KratosMultiphysics.VariableUtils()
        
        self.velocity_application_process = assign_vector_components_to_nodes_process.AssignVectorComponentsToNodesProcess(Model,settings)
        
        self.model_part = Model[settings["model_part_name"].GetString()]
        
        self.interval = KratosMultiphysics.IntervalUtility(settings)
        
        self.fix_x = settings["constrained"][0].GetBool()
        self.fix_y = settings["constrained"][1].GetBool()
        self.fix_z = settings["constrained"][2].GetBool()
        if(settings["value"][0].IsNull()):
            self.fix_x = False
        if(settings["value"][1].IsNull()):
            self.fix_y = False
        if(settings["value"][2].IsNull()):
            self.fix_z = False
            
        self.interval_is_active = False

    def ExecuteInitializeSolutionStep(self):
        self.velocity_application_process.ExecuteInitializeSolutionStep()
        
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        
        self.interval_is_active = self.interval.IsInInterval(current_time)
        
        if(self.interval_is_active):
            fixed = True
            if(self.fix_x):
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_X, fixed, self.model_part.Nodes)
            if(self.fix_y):
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_Y, fixed, self.model_part.Nodes)
            if(self.fix_z):
                self.variable_utils.ApplyFixity(KratosMultiphysics.DISPLACEMENT_Z, fixed, self.model_part.Nodes)

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
