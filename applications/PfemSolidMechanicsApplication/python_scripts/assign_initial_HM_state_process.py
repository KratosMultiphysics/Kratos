import KratosMultiphysics
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPFEMSolid

## This process sets the initial value of stress and water pressure to the domain

def Factory( custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetMechanicalInitialStateProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class SetMechanicalInitialStateProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)



        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
         "model_part_name": "Main_Domain",
         "gravity_active": false, 
         "constant_vertical_stress": -50.0,
         "constant_horizontal_stress": -50.0,
         "constant_water_pressure" : -0.0,
         "top_surface_load_bool": false,
         "top_surface_load": 0.0,
         "top_water_pressure": 0.0
        }
        """)
       
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.restarted = False;
        self.the_process_has_been_executed = False;


        ## 

    ##def ExecuteBeforeSolutionLoop(self):
    def ExecuteThisProcess(self):

        self.model_part = self.model_part[self.model_part_name]
        self.restarted = self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]
        if ( self.restarted == True):
            return
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", self.settings["model_part_name"])
        params.AddValue("gravity_active",self.settings["gravity_active"])
        params.AddValue("constant_vertical_stress",self.settings["constant_vertical_stress"])
        params.AddValue("constant_horizontal_stress",self.settings["constant_horizontal_stress"])
        params.AddValue("constant_water_pressure",self.settings["constant_water_pressure"])
        params.AddValue("top_surface_load_bool",self.settings["top_surface_load_bool"])
        params.AddValue("top_surface_load",self.settings["top_surface_load"])
        params.AddValue("top_water_pressure",self.settings["top_water_pressure"])
        initial_state_process = KratosPFEMSolid.SetMechanicalInitialStateProcess(self.model_part, self.settings)
        initial_state_process.Execute()
        self.the_process_has_been_executed = True

    def ExecuteInitializeSolutionStep(self):

        if ( self.the_process_has_been_executed == False):
            self.ExecuteThisProcess() 


