from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.IGAStructuralMechanicsApplication import *

KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, model_part):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NurbsBrepProcess(model_part, settings)

class NurbsBrepProcess(KratosMultiphysics.Process):
  
    def __init__(self,model_part,params):

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "nodal_results"            : [],
            "integration_point_results": [],
            "output_file_name"         : "",
            "model_part_name"          : "",
            "sub_model_part_name"      : "",
            "time_frequency"           : 0.10,
            "write_points": {
              "sub_model_part_name": [ ],
              "output_file_name": "convergence.txt"
            }
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)
        
        self.model_part = model_part
        
        self.output_file_name = ""
        self.nodal_results = []
        self.integration_point_results = []
        self.frequency = 0.0
        self.time_counter = 0.0

    def ExecuteInitialize(self):
        self.output_file_name = self.params["output_file_name"].GetString()
        if (len(self.params["sub_model_part_name"].GetString()) > 0):
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()].GetSubModelPart(self.params["sub_model_part_name"].GetString())
        else:
            self.sub_model_part = self.model_part[self.params["model_part_name"].GetString()]
            
        self.nodal_results = self.__generate_variable_list_from_input(self.params["nodal_results"])
        self.integration_point_results = self.__generate_variable_list_from_input(self.params["integration_point_results"])
        self.frequency = self.params["time_frequency"].GetDouble()

        self.step = 0

        with open(self.output_file_name, 'w') as file:
            file.write("Rhino Post Results File 1.0\n") 
            #for i in range(self.params["nodal_results"].size()):
            #    out = self.params["nodal_results"][i]
            #    variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )

            #    file.write("Result \"" + out.GetString() + "\" \"Load Case\" 0  Vector OnNodes\n")
            #    file.write("Values\n")
            #    for node in self.sub_model_part.Nodes:
            #        value = node.GetSolutionStepValue(variable, 0)
            #        if isinstance(value,float):
            #            file.write(str(node.Id) + "  " + str(value))
            #        else: # It is a vector
            #            file.write(str(node.Id) + "  " + str(value[0]) + "  " + str(value[1]) + "  " + str(value[2]) + "\n")

            #    file.write("End Values\n")

            #for i in range(self.params["integration_point_results"].size()):
            #    out = self.params["integration_point_results"][i]
            #    variable = KratosMultiphysics.KratosGlobals.GetVariable( out.GetString() )

            #    value = self.sub_model_part.Elements[0].Calculate(variable, self.sub_model_part.ProcessInfo)
            #    if isinstance(value,float):
            #        type = "Scalar"
            #    else:
            #        type = "Vector"

            #    file.write("Result \"" + out.GetString() + "\" \"Load Case\" 0  " + type + " OnGaussPoints\n")
            #    file.write("Values\n")
            #    for element in self.sub_model_part.Elements:
            #        value = 0.0#element.Calculate(variable, self.sub_model_part.ProcessInfo)

            #        if isinstance(value,float):
            #            file.write(str(element.Id) + "  " + str(value) + "\n")
            #        else: # It is a vector
            #            file.write(str(element.Id) + "  " + str(value[0]) + "  " + str(value[1]) + "  " + str(value[2]) + "\n")
            #    file.write("End Values\n")

        #if(self.params.Has("write_points")):
        #    with open(self.params["write_points"]["output_file_name"].GetString(), 'w') as convergence_file:
        #        convergence_file.write("")
        
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
        
        if True:#self.time_counter + 0.00001 >= self.frequency:
            #self.time_counter = 0.0

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
                        #elementList = self.sub_model_part.ElementsArray(0)
                        #element = elementList[i]
                        value = element.Calculate(variable, self.sub_model_part.ProcessInfo)

                        if isinstance(value,float):
                            file.write(str(element.Id) + "  " + str(value) + "\n")
                        else: # It is a vector
                            file.write(str(element.Id) + "  " + str(value[0]) + "  " +  str(value[1]) + "  " + str(value[2]) + "\n")
                    file.write("End Values\n")

        if(self.params.Has("write_points")):
            with open(self.params["write_points"]["output_file_name"].GetString(), 'a') as convergence_file:
                for i in range(0, self.params["write_points"]["sub_model_part_name"].size()):
                    check_model_part = self.model_part[self.params["model_part_name"].GetString()].GetSubModelPart(self.params["write_points"]["sub_model_part_name"][i].GetString())
                    convergence_file.write(str(self.sub_model_part.NumberOfElements()) + "  " + str(self.sub_model_part.NumberOfNodes()) + "  ")
                    for condition in check_model_part.Conditions:
                        disp = condition.Calculate(KratosMultiphysics.DISPLACEMENT, check_model_part.ProcessInfo)
                        convergence_file.write(str(disp[0]) + "  " + str(disp[1]) + "  " + str(disp[2]) + "  ")
                    convergence_file.write("\n")

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