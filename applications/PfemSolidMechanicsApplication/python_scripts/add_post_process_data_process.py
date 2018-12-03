import KratosMultiphysics
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPFEMSolid

## This process sets the initial value of stress and water pressure to the domain

def Factory( custom_settings, Model):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AddAdditionalPostProcessDataProcess(Model, custom_settings["Parameters"])

## All the processes python processes should be derived from "python_process"
class AddAdditionalPostProcessDataProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):
        KratosMultiphysics.Process.__init__(self)



        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
         "model_part_name": "Main_Domain",
         "nodal_data": ["EXCESS_WATER_PRESSURE"],
         "gravity_active": false,
         "constant_water_pressure": -50.0,
         "top_surface_load": false,
         "top_water_pressure": -50.0,
         "top_ref_levelY": 0.0
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.gravity_active = self.settings["gravity_active"].GetBool()
        if (self.gravity_active):
            self.ref_levelY = self.settings["top_ref_levelY"].GetDouble()
            self.top_load_active = self.settings["top_surface_load"].GetBool()
            if (self.top_load_active):
                self.top_wP = self.settings["top_water_pressure"].GetDouble()
            else:
                self.top_wP = 0.0
        else:
            self.constant_wP = self.settings["constant_water_pressure"].GetDouble()


        ## 

    def ExecuteBeforeOutputStep(self):
        self.model_part = self.model_part[self.model_part_name]
        for node in self.model_part.Nodes:
            delta_wp = node.GetSolutionStepValue( KratosMultiphysics.WATER_PRESSURE )
            initial_wp = self._CalucalteInitialWaterPressure(node)
            delta_wp -= initial_wp

            if (node.SolutionStepsDataHas(KratosPFEMSolid.EXCESS_WATER_PRESSURE)):
                node.SetSolutionStepValue(KratosPFEMSolid.EXCESS_WATER_PRESSURE,delta_wp)

    def _CalucalteInitialWaterPressure(self, node):
        if (self.gravity_active):
            wP0 = 10.0*1.0*(node.Y-self.ref_levelY) + self.top_wP
            return wP0
        else:
            return self.constant_water_pressure


