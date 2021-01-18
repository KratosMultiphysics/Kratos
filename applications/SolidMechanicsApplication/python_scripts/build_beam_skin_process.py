import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

def Factory(custom_settings, Model):
    if( not isinstance(custom_settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return BeamBuildSkinProcess(Model, custom_settings["Parameters"])

class BeamBuildSkinProcess(KratosMultiphysics.Process):
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "model_part_name" : "BeamDomain",
	    "number_of_sides" : 4,
	    "diameter": 0.1,
            "echo_level" : 0
        }
        """)

        #overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model
        self.echo_level = self.settings["echo_level"].GetInt()

        #get section sides
        self.sides = self.settings["number_of_sides"].GetInt()

        #get section radius
        self.radius = 0.5 * self.settings["diameter"].GetDouble()

    #
    def ExecuteInitialize(self):
        self.model_part = self.model[self.settings["model_part_name"].GetString()]

        #create skin for the beam string
        self.beam_skin_process = KratosSolid.BuildStringSkinProcess(self.model_part, self.sides, self.radius)

        self.beam_skin_process.ExecuteInitialize()
        print("::[----Skin_Created---]:: (sides:"+str(self.sides),"radius:"+str(self.radius)+")")


    #
    def ExecuteBeforeSolutionLoop(self):
        self.beam_skin_process.ExecuteBeforeSolutionLoop()
        if( self.echo_level > 0 ):
            print("::[----Skin_Created---]::Transfer to output")

    #
    def ExecuteInitializeSolutionStep(self):
        pass

    #
    def ExecuteFinalizeSolutionStep(self):
        self.beam_skin_process.ExecuteFinalizeSolutionStep()
        if( self.echo_level > 0 ):
            print("::[----Skin_Updated---]::", self.settings["section_type"].GetString())

    #
    def Execute(self):
        pass

    #
    def ExecuteFinalize(self):
        pass

    #
    def ExecuteBeforeOutputStep(self):
        self.beam_skin_process.ExecuteBeforeOutputStep()

    #
    def ExecuteAfterOutputStep(self):
        self.beam_skin_process.ExecuteAfterOutputStep()
