import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager
from KratosMultiphysics import assign_vector_variable_process
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMVelocityBoundaryProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMVelocityBoundaryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "variable_name"             : "VELOCITY",
                "value"                     : [0.0, "0*t", 0.0],
                "interval"                  : [0.0, 1e30],
                "option"                    : "",   
                "local_axes"                : {}
            }  """ )

         # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

       
        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.model_part_name = settings["model_part_name"].GetString()
        self.model = Model
        self.option = settings["option"].GetString()

        # get variable imposed and check
        variable_name = settings["variable_name"].GetString()
        variable_name_list = ["DISPLACEMENT","VELOCITY"]
        if(variable_name in variable_name_list):
            self.variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        else:
            err_msg =  "The given variable \"" + variable_name + "\" is not available to be imposed with this process.\n"
            err_msg += "Available options are: " + ", ".join(variable_name_list)
            raise Exception(err_msg)

        self.value_is_numeric = [False, False, False]
        self.value = KratosMultiphysics.Vector(3)
        self.aux_function = ["0.0","0.0","0.0"]
        self.name = ["0.0","0.0","0.0"]
        # Loop over components X, Y and Z
        for i, variable in enumerate(["_X", "_Y", "_Z"]):
            self.name[i] = settings["variable_name"].GetString() + variable
            if settings["value"][i].IsNumber():
                self.value_is_numeric[i] = True
                self.value[i] = settings["value"][i].GetDouble()
            else:
                self.function_string = settings["value"][i].GetString()
                self.aux_function[i] = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])

        self.model_part = self.model[self.model_part_name]
        self.variable_utils = KratosMultiphysics.VariableUtils()


    #def ApplyBoundaryConditions(self):
     #   super().ApplyBoundaryConditions()

        ### Apply manufactured solution BCs
        # Set the manufactured solution source terms
     #   self.ExecuteInitializeSolutionStep()
     
    def ExecuteBeforeSolutionLoop(self):
     #   """ This method is executed in before initialize the solution step

     #   Keyword arguments:
     #   self -- It signifies an instance of a class.
     #   """
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

        for node in self.model_part.Nodes:
            #mpc_coord = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]
            acceleration = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION,1)
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,1)
            displacement= node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,1)
            aux_displacement = self.value * delta_time - displacement + (0.5 * acceleration * delta_time * delta_time); # 
            #aux_displacement = (self.value * delta_time) # + (0.5 * acceleration * delta_time * delta_time); # 
            #print(delta_time)

            if self.interval.IsInInterval(current_time):
                self.step_is_active = True

                # Loop over components X, Y and Z
                #for i in range(3):
                    #self.variable = self.name[i]
                    #if  not self.value_is_numeric[i]:
                    #    if self.aux_function[i].DependsOnSpace() == False: #depends on time only
                    #        self.value[i] = self.aux_function[i].CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                    #    else: #most general case - space varying function (possibly also time varying)
                    #        self.value[i] = self.aux_function[i].CallFunction(mpc_coord[0],mpc_coord[1],mpc_coord[2],current_time,0.0,0.0,0.0)
            
            #self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,aux_displacement)

            print(aux_displacement)

            #self.variable_utils.SetVariable(KratosMultiphysics.DISPLACEMENT, aux_displacement, self.model_part.Nodes)
            #print(node)

        #for process in self.aux_processes:
        #    process.ExecuteInitializeSolutionStep()
       

    def ExecuteFinalizeSolutionStep(self):
    #    """ This method is executed in order to finalize the current step

    #    Keyword arguments:
    #    self -- It signifies an instance of a class.
    #    """
        #self.ExecuteInitializeSolutionStep()
    #    for process in self.aux_processes:
    #        process.ExecuteFinalizeSolutionStep()
        


        for node in self.model_part.Nodes:
            #mpc_coord = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]
            #self.variable_utils.SetVariable(self.variable, self.value, self.mesh.Nodes)
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            #node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,0,aux_displacement
