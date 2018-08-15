from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
KratosMultiphysics.CheckForPreviousImport()


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return RestartProcess(Model, settings["Parameters"])


class RestartProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[custom_settings["model_part_name"].GetString()]

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"     : "DEFINE_MODEL_PART_NAME",
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
        self.output_frequency  = self.settings["output_frequency"].GetDouble()

        self.output_label_is_time = False
        output_file_label = self.settings["restart_file_label"].GetString()
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
        self.problem_path = os.getcwd()
        self.restart_path = os.path.join(self.problem_path, self.settings["restart_file_name"].GetString() )

        self.step_count         = 0
        self.counter            = 1
        self.printed_step_count = 0
        self.next_output        = 0

        self.echo_level  = 1
        if( self.echo_level > 0 ):
            print("::[------Restart------]:: -BUILT-")
    #
    def ExecuteInitialize(self):

        # Set current time parameters
        if( self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.step_count = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
            self.printed_step_count = self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_RESTART_STEP]

            if self.output_control_is_time:
                self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.TIME] + self.output_frequency
            else:
                self.next_output = self.model_part.ProcessInfo[KratosMultiphysics.STEP] + self.output_frequency
        else:
            self.next_output = self.output_frequency

        # Copy to a restart folder the posterior files, delete from problem folder
        if( self.save_restart ):
            if self.output_label_is_time:
                restart_input_label = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
                self.CleanPosteriorFiles(restart_input_label)
            else:
                restart_input_label = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
                self.CleanPosteriorFiles(restart_input_label)
        else:
            self.CleanPreviousFiles()

        # Reconstruct .list and .csv (last build using a restart ID to be rebuild)

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1

    #
    def ExecuteFinalizeSolutionStep(self):
        pass

    #
    def ExecuteAfterOutputStep(self):
        if( self.save_restart ):
            if(self.IsRestartStep()):
                self.SaveRestart()


    ###

    #
    def SaveRestart(self):

        # Print the output
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.printed_step_count += 1
        self.model_part.ProcessInfo[KratosMultiphysics.PRINTED_RESTART_STEP] = self.printed_step_count
        if self.output_label_is_time:
            label = time
        else:
            label = self.printed_step_count

        if( self.echo_level > 0 ):
            print("::[------Restart------]:: Save [ step:", step,"] [ label:", label,"]")

        current_restart_path = self.restart_path + "__" + str(label)

        # set serializer flag
        self.serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_NO_TRACE  # binary
        # self.serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ERROR # ascii
        # self.serializer_flag = KratosMultiphysics.SerializerTraceType.SERIALIZER_TRACE_ALL   # ascii

        serializer = KratosMultiphysics.Serializer(current_restart_path, self.serializer_flag)

        serializer.Save(self.model_part.Name, self.model_part)

        self.counter += 1

        # schedule next output
        if(self.output_frequency > 0.0): # note: if == 0 always active
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
            #print( str(self.model_part.ProcessInfo[KratosMultiphysics.TIME])+">"+ str(self.next_output) )
            return ( self.model_part.ProcessInfo[KratosMultiphysics.TIME] >= self.next_output )
        else:
            return ( self.step_count >= self.next_output )


    #
    def CleanPreviousFiles(self):

        print(" ::[------Restart------]:: Remove Previous Files")

        # remove previous results:
        file_ending_type = ".post.bin"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPreviousFileType(file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous graph files:
        file_ending_type = ".graph.png"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous data files:
        file_ending_type = ".post.csv"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous list files:
        file_ending_type = ".post.lst"
        self.CleanPreviousFileType(file_ending_type)

        # remove previous restart files:
        file_ending_type = ".rest"
        self.CleanPreviousFileType(file_ending_type)

    #
    def CleanPosteriorFiles(self, restart_step):

        print(" ::[------Restart------]:: Clean Post Restart Files) [ STEP:",restart_step," ]")

        # remove posterior results after restart:
        file_ending_type = ".post.bin"
        self.CleanPosteriorFileType(restart_step, file_ending_type)

        file_ending_type = ".post.res"
        self.CleanPosteriorFileType(restart_step, file_ending_type)

        file_ending_type = ".post.msh"
        self.CleanPosteriorFileType(restart_step, file_ending_type)

        # remove posterior graphs after restart:

        file_ending_type = "graph.png"
        self.CleanPosteriorFileType(restart_step+1, file_ending_type)

        # remove posterior restart files:
        file_ending_type = ".rest"
        self.CleanPosteriorFileType(restart_step+1, file_ending_type)


    #
    def CleanPreviousFileType(self, file_ending_type):

        if(os.path.exists(self.problem_path) == False):
            print("::[------Restart------]:: Problem Path does not exist , check the Problem Path selected ")
        else:
            filelist = [f for f in os.listdir(self.problem_path) if f.endswith(file_ending_type)]

            for f in filelist:
                try:
                    os.remove(f)
                except OSError:
                    pass

    #
    def IsInt(self,s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    #
    def CleanPosteriorFileType(self, restart_step, file_ending_type):

        if(os.path.exists(self.problem_path) == False):
            print(" ::[------Restart------]:: Problem Path does not exist , check the Problem Path selected ")
        else:
            filelist = []
            for f in os.listdir(self.problem_path):
                if(f.endswith(file_ending_type)):
                    #if f.name = problem_tested_145.post.bin
                    file_parts = f.split('_')  # you get ["problem","tested","145.post.bin"]
                    num_parts  = len(file_parts)

                    end_parts  = file_parts[num_parts-1].split(".") # you get ["145","post","bin"]
                    print_id   = end_parts[0] # you get "145"

                    if( self.IsInt(print_id) ):
                        if( int(print_id) > int(restart_step) ):
                            filelist.append(f)

            for f in filelist:
                try:
                    os.remove(f)
                except OSError:
                    pass
