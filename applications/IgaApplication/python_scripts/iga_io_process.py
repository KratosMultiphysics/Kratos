from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IgaApplication as IgaApplication


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IGA_IO_Process(Model, settings)

class IGA_IO_Process(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_results"            : [],
            "integration_point_results": [],
            "output_file_name"         : "",
            "model_part_name"          : "",
            "time_frequency"           : 0.01
            }
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        
        self.model_part = Model[settings["model_part_name"].GetString()]

        self.output_file_name = self.params["output_file_name"].GetString()
        self.nodal_results = self.__generate_variable_list_from_input(self.params["nodal_results"])
        self.integration_point_results = self.__generate_variable_list_from_input(self.params["integration_point_results"])

        self.frequency = self.params["time_frequency"].GetDouble()

    def ExecuteInitialize(self):
        self.time_counter = 0.0
        self.step = 0

        with open(self.output_file_name, 'w') as file:
            file.write("Rhino Post Results File 1.0\n") 

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.step = self.step + 1
        step = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME_STEPS)
        self.time_counter += dt

        if self.time_counter + 0.00001 >= self.frequency:
            self.time_counter = 0.0

            with open(self.output_file_name, 'a') as file:
                for i in range(self.params["nodal_results"].size()):
                    out = self.params["nodal_results"][i]
                    variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )

                    if type(variable) == KratosMultiphysics.DoubleVariable:
                        type_name = "Scalar"
                    else:
                        type_name = "Vector"

                    file.write("Result \"" + out.GetString() + "\" \"Load Case\" " + str(self.step) + " " + type_name + " OnNodes\n")

                    file.write("Values\n")
                    for node in self.sub_model_part.Nodes:
                        value = node.GetSolutionStepValue(variable, 0)
                        if isinstance(value,float):
                            file.write(str(node.Id) + "  " + str(value) + "\n")
                        else: # It is a vector
                            file.write(str(node.Id) + "  " + str(value[0]) + "  " + str(value[1]) + "  " + str(value[2]) + "\n")
                    file.write("End Values\n")

                for i in range(self.params["integration_point_results"].size()):
                    out = self.params["integration_point_results"][i]
                    variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )

                    if type(variable) == KratosMultiphysics.DoubleVariable:
                        type_name = "Scalar"
                    else:
                        type_name = "Vector"

                    file.write("Result \"" + out.GetString() + "\" \"Load Case\" " + str(self.step) + " " + type_name + " OnGaussPoints\n")

                    file.write("Values\n")
                    for element in self.sub_model_part.Elements:
                        value = element.Calculate(variable, self.sub_model_part.ProcessInfo)
                        if isinstance(value,float):
                            file.write(str(element.Id) + "  " + str(value) + "\n")
                        else: # It is a vector
                            file.write(str(element.Id) + "  " + str(value[0]) + "  " +  str(value[1]) + "  " + str(value[2]) + "\n")
                    file.write("End Values\n")

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass

    def __generate_variable_list_from_input(self,param):
        '''Parse a list of variables from input.'''
        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]