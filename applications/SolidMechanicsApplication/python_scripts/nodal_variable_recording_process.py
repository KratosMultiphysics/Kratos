from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
# importing the Kratos Library
import KratosMultiphysics
KratosMultiphysics.CheckForPreviousImport()

from numpy import *

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NodalVariableRecordingProcess(Model, settings["Parameters"])


class NodalVariableRecordingProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"     : "MODEL_PART_NAME",
            "variable"            : "VARIABLE_NAME",
            "output_file_name"    : "FILE_NAME",
            "single_entities"     : true
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # model
        self.model = Model

        # sum all entities
        self.single_entities = self.settings["single_entities"].GetBool()

        # restart path
        self.problem_path = os.getcwd()

        # variable
        self.variable = self.settings["variable"].GetString()
        self.kratos_variable = KratosMultiphysics.KratosGlobals.GetVariable(self.variable)

        self.echo_level  = 1
        if( self.echo_level > 0 ):
            print("::[RecordVariable]:: -BUILT-")
    #
    def ExecuteInitialize(self):

        self.model_part = self.model.GetModelPart(custom_settings["model_part_name"].GetString())

        # Set path and headers
        if( self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):

            current_step = self.GetCurrentStep()

            self.recording_path_list = []

            if( not self.single_entities ):
                self.recording_path = os.path.join(self.problem_path, self.settings["output_file_name"].GetString() + ".post.csv")
                self.CleanPostRestartData(self.recording_path, current_step)

            else:
                for node in self.model_part.GetNodes():
                    self.recording_path = os.path.join(self.problem_path, self.settings["output_file_name"].GetString() + "_node_" + str(node.Id) + ".post.csv")

                    self.recording_path_list.append(self.recording_path)

                    self.CleanPostRestartData(self.recording_path, current_step)

        else:
            self.InitializeRecordingHeaders()


    #
    def InitializeRecordingHeaders(self):

        self.recording_path_list = []

        if( not self.single_entities ):

            first_node = self.first(self.model_part.GetNodes())

            self.recording_path = os.path.join(self.problem_path, self.settings["output_file_name"].GetString() + ".post.csv")

            line_header = "R,STEP,TIME"

            variable_value = first_node.GetSolutionStepValue(self.kratos_variable)

            if( hasattr(variable_value, "__len__") ):
                if( len(variable_value) == 3 ):
                    line_header = line_header + "," + str(self.variable) + "," + str(self.variable) + "_X," + str(self.variable) + "_Y," + str(self.variable) + "_Z"
                else:
                    line_header = line_header + "," + str(self.variable)
            else:
                line_header = line_header + "," + str(self.variable)


            line_header  = line_header + "\n"

            if(os.path.exists(self.recording_path) == True):
                source = self.recording_path
                os.remove(source)

            # print file header
            recording_file = open(self.recording_path, "w")
            recording_file.write(line_header)
            recording_file.close()


        else:

            for node in self.model_part.GetNodes():
                self.recording_path = os.path.join(self.problem_path, self.settings["output_file_name"].GetString() + "_node_" + str(node.Id) + ".post.csv")

                self.recording_path_list.append(self.recording_path)

                line_header = "R,STEP,TIME"

                variable_value = node.GetSolutionStepValue(self.kratos_variable)

                if( hasattr(variable_value, "__len__") ):
                    if( len(variable_value) == 3 ):
                        line_header = line_header + "," + str(self.variable) + "," + str(self.variable) + "_X," + str(self.variable) + "_Y," + str(self.variable) + "_Z"
                    else:
                        line_header = line_header + "," + str(self.variable)
                else:
                    line_header = line_header + "," + str(self.variable)


                line_header  = line_header + "\n"

                if(os.path.exists(self.recording_path) == True):
                    source = self.recording_path
                    os.remove(source)

                # print file header
                recording_file = open(self.recording_path, "w")
                recording_file.write(line_header)
                recording_file.close()


    # Reconstruct .list and .csv (last build using a restart ID to be rebuild)
    #
    def CleanPostRestartData(self, file_path, restart_step):

        if os.path.exists(file_path):

            text = "R," + str(restart_step) + ","

            #if printing all steps:
            #line_number = restart_step+1

            line_number = self.get_line_number(text,file_path)
            print(" Text found: ", text, " in line ", line_number)

            if( line_number != None ):
                self.truncate_file_in_line(line_number,file_path)


    #
    @classmethod
    def truncate_file_in_line(self, line_number, file_path):

        source  = open(file_path, "r+")

        source.seek(0)

        for line in range(line_number):
            source.readline()

        line_number = source.tell()
        source.truncate(line_number)
        source.close()


    #
    @classmethod
    def get_line_number(self, text, file_path):

        with open(file_path) as f:
            for i, line in enumerate(f, 1):
                if text in line:
                    return i

    ###

    #
    def ExecuteFinalizeSolutionStep(self):
        self.RecordStepVariable()



    ###


    #
    def GetCurrentTime(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.TIME]

    #
    def GetCurrentStep(self):

        return self.model_part.ProcessInfo[KratosMultiphysics.STEP]

    #
    def GetVectorVariableValueSum(self, nodal_variable):

        variable_value = []
        counter = 0
        for node in self.model_part.GetNodes():
            nodal_value = node.GetSolutionStepValue(nodal_variable);
            variable_value = variable_value + nodal_value
            counter+=1

        return variable_value

    #
    def GetScalarVariableValueSum(self, nodal_variable):

        variable_value = 0
        for node in self.model_part.GetNodes():
            nodal_value = node.GetSolutionStepValue(nodal_variable);
            variable_value = variable_value + nodal_value

        return variable_value

    #
    @classmethod
    def GetVariableModulus(self, nodal_variable):

        modulus = 0
        for var in nodal_variable:
            modulus = modulus + var * var

        return sqrt(modulus)


    #
    @classmethod
    def first(self, iterable, condition = lambda x: True):

        for item in iterable:
            if condition(item):
                return item

        raise ValueError('No satisfactory value found')

    #
    def RecordStepVariable(self):

        current_time = self.GetCurrentTime()
        current_step = self.GetCurrentStep()

        if( not self.single_entities ):

            first_node = self.first(self.model_part.GetNodes())

            line_value = "R," + str(current_step) + "," + str(current_time)

            variable_value = first_node.GetSolutionStepValue(self.kratos_variable)

            if( hasattr(variable_value, "__len__") ):
                variable_value   = self.GetVectorVariableValueSum(self.kratos_variable)
            else:
                variable_value   = self.GetScalarVariableValueSum(self.kratos_variable)

            if( hasattr(variable_value, "__len__") ):
                if( len(variable_value) == 3 ):
                    variable_modulus = self.GetVariableModulus(variable_value)
                    line_value = line_value + "," + str(variable_modulus) + "," + str(variable_value[0]) + "," + str(variable_value[1]) + "," + str(variable_value[2])
                else:
                    line_value = line_value + "," + str(variable_value)
            else:
                line_value = line_value + "," + str(variable_value)


            line_value  = line_value + "\n"

            recording_file = open(self.recording_path, "a")
            recording_file.write(line_value)
            recording_file.close()

        else:

            counter = 0
            for node in self.model_part.GetNodes():

                self.recording_path = self.recording_path_list[counter]

                line_value = "R," + str(current_step) + "," + str(current_time)

                variable_value = node.GetSolutionStepValue(self.kratos_variable)

                if( hasattr(variable_value, "__len__") ):
                    if( len(variable_value) == 3 ):
                        variable_modulus = self.GetVariableModulus(variable_value)

                        line_value = line_value + "," + str(variable_modulus) + "," + str(variable_value[0]) + "," + str(variable_value[1]) + "," + str(variable_value[2])
                    else:
                        line_value = line_value + "," + str(variable_value)
                else:
                    line_value = line_value + "," + str(variable_value)


                line_value  = line_value + "\n"


                recording_file = open(self.recording_path, "a")
                recording_file.write(line_value)
                recording_file.close()
                counter += 1
