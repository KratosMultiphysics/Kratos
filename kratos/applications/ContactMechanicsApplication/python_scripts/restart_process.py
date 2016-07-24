from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
KratosMultiphysics.CheckForPreviousImport()


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RestartProcess(Model, settings["Parameters"])


class RestartProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)
        
        self.main_model_part = Model[custom_settings["main_model_part_name"].GetString()]
    
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "save_restart"        : true,
            "restart_file_name"   : "problem_name",
            "restart_file_label"  : "step",
            "output_control_type" : "step",
            "output_frequency"    : 1.0,
            "json_output"         : false
        }
        """)
 
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.save_restart = self.settings["save_restart"].GetBool()

        # Set up output frequency and format
        self._frequency  = self.settings["output_frequency"].GetDouble()

        self.output_label_is_time = False
        output_file_label = result_file_configuration["restart_file_label"].GetString()
        if(output_file_label == "time"):
            self.output_label_is_time = True
        elif(output_file_label == "step"):
            self.output_label_is_time = False
      
        self.output_control_is_time = False
        output_control_type  = self.settings["output_control_type"].GetString()
        if(output_control_type == "time"):
            self.output_control_is_time = True
        elif(output_control_type == "step"):
            self.output_control_is_time = False

        # restart path
        problem_path = os.getcwd()
        self.restart_path = os.path.join(problem_path, self.settings["restart_file_name"].GetString() )

        self.step_count  = 0
        self.printed_step_count = 0
        self.next_output = 0.0
                       
    #
    def ExecuteInitialize(self):
        pass

    ###

    #
    def ExecuteInitializeSolutionStep(self):
        pass

    #
    def ExecuteFinalizeSolutionStep(self):
        
        if( self.save_restart ):
            self.SaveRestart()

    ###

    #
    def SaveRestart(self):
        
        if( self.echo_level > 0 ):
            print("::[Restart_Process]:: RESTART SAVED...", self.counter)

        # Print the output
        time = self.model_part.ProcessInfo[TIME]
        self.printed_step_count += 1
        if self.output_label_is_time:
            label = time
        else:
            label = self.printed_step_count

        self.restart_path = os.path.join(self.restart_path, label)

        # set serializer flag
        self.serializer_flag = "SERIALIZER_NO_TRACE"      # binary
        # self.serializer_flag = "SERIALIZER_TRACE_ERROR" # ascii
        # self.serializer_flag = "SERIALIZER_TRACE_ALL"   # ascii
        kratos_serializer_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.serializer_flag)
        
        serializer = Serializer(restart_path, kratos_serializer_variable)
        
        serializer.Save(self.main_model_part.GetModelPartName(), self.main_model_part)

        # schedule next output
        if(self.output_frequency >= 0.0):
            if(self.output_control_is_time):
                while(self.next_output <= time):
                    self.next_output += self.output_frequency
            else:
                while(self.next_output <= self.step_count):
                    self.next_output += self.output_frequency



    #
    def GetRestartStep(self):
        return self.counter

    #
    def IsRestartStep(self):

        if(self.output_control_is_time):
            #print( str(self.model_part.ProcessInfo[TIME])+">"+ str(self.next_output) )
            return ( self.model_part.ProcessInfo[TIME] > self.next_output )
        else:
            return ( self.step_count >= self.next_output )
