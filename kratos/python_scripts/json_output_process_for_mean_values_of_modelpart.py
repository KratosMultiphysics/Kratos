from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
#from json_utilities import *
import json
#KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return JsonOutputProcess(Model, settings["Parameters"])

# All the processes python processes should be derived from "Process"

class JsonOutputProcess(KratosMultiphysics.Process):
    """This class is used in order to create a json file containing
    the solution a given model part with a certain frequency

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model_part -- the model part used to construct the process.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model_part, params):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model_part -- the model part used to construct the process.
        settings -- Kratos parameters containing solver settings.
        """
        KratosMultiphysics.Process.__init__(self)

        ## Settings string in json format
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "help"                          : "This process generates a json file containing the solution of a list of variables from a given submodelpart",
            "output_variables"              : [],
            "gauss_points_output_variables" : [],
            "output_file_name"              : "",
            "model_part_name"               : "",
            "sub_model_part_name"           : "",
            "check_for_flag"                : "",
            "time_frequency"                : 1.00,
            "historical_value"              : true,
            "resultant_solution"            : false
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model_part

        self.params = params

        self.output_file_name = ""
        self.output_variables = []
        self.gauss_points_output_variables = []
        self.frequency = 0.0
        self.time_counter = 0.0
        self.resultant_solution = False

    def ExecuteInitialize(self):
        """ This method is executed at the begining to initialize the process

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # We get the submodelpart
        model_part_name = self.params["model_part_name"].GetString()
        sub_model_part_name = self.params["sub_model_part_name"].GetString()
        if (sub_model_part_name != ""):
            self.sub_model_part = self.model_part[model_part_name].GetSubModelPart(sub_model_part_name)
        else:
            self.sub_model_part = self.model_part[model_part_name]

        # We get the number of subsubmodelparts

        if self.sub_model_part.NumberOfSubModelParts() == 0:
            self.number_of_segments = 1
        else:
            self.number_of_segments = self.sub_model_part.NumberOfSubModelParts()        

        # If we consider any flag
        flag_name = self.params["check_for_flag"].GetString()
        if (flag_name != ""):
            self.flag = globals().get(flag_name)
        else:
            self.flag = None

        self.output_file_name = self.params["output_file_name"].GetString()
        self.output_variables = self.__generate_variable_list_from_input(self.params["output_variables"])
        self.gauss_points_output_variables = self.__generate_variable_list_from_input(self.params["gauss_points_output_variables"])
        self.frequency = self.params["time_frequency"].GetDouble()
        self.resultant_solution = self.params["resultant_solution"].GetBool()
        self.historical_value = self.params["historical_value"].GetBool()

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed before starting the time loop

        This step generates the structure of the dictionary

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        data = {}
        data["TIME"] = []

        # Gauss point values
        for i in range(self.params["gauss_points_output_variables"].size()):
            out = self.params["gauss_points_output_variables"][i]
            variable_name = out.GetString()
            variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
            variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

            data["Mean_Value_of_" + variable_name] = {}

            for seg in range(self.number_of_segments):
                current_segment = self.sub_model_part.GetSubModelPart("Segment_" + str(seg+1))

                for elem in current_segment.Elements:
                    compute = self.__check_flag(elem)

                    if (compute == True):
                        data["Mean_Value_of_" + variable_name]["Segment_" + str(seg+1)]= []

        #write_external_json(self.output_file_name, data)
        with open(self.output_file_name, 'w') as outfile:
            json.dump(data, outfile)

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Here the dictionary containing the solution is filled

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        #data =  read_external_json(self.output_file_name)
        with open(self.output_file_name, 'r') as outfile:
            data = json.load(outfile)
  


        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt

        # Compute the total Area of the Model_part


        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            data["TIME"].append(time)

            #Gauss points values
            for i in range(self.params["gauss_points_output_variables"].size()):
                out = self.params["gauss_points_output_variables"][i]
                variable_name = out.GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                for seg in range(self.number_of_segments):
                    current_segment = self.sub_model_part.GetSubModelPart("Segment_" + str(seg+1))
                    segment_area = 0.0
                    
                    for elem in current_segment.Elements:
                        segment_area += current_segment.GetElement(elem.Id).GetGeometry().Area()


                    total_list_elem_participations = []

                    for elem in current_segment.Elements:
                        compute = self.__check_flag(elem)

                        if (compute == True):
                            value = elem.CalculateOnIntegrationPoints(variable, current_segment.ProcessInfo)

                            # print(variable_name)
                            # print(value[0])

                            gauss_point_number = len(value)
                        
                            if variable_type == "Vector":
                                elem_mean_list = self.__generate_elemental_gp_mean_value(value)
                                elem_area = current_segment.GetElement(elem.Id).GetGeometry().Area()
                                elem_proportion_total_mean = elem_area / segment_area
                                elem_total_mean_participation = [x * elem_proportion_total_mean for x in elem_mean_list]

                            total_list_elem_participations.append(elem_total_mean_participation)

                    list = self.__generate_total_sub_model_part_mean(total_list_elem_participations)

                    data["Mean_Value_of_" + variable_name]["Segment_" + str(seg+1)].append(list)

        #write_external_json(self.output_file_name, data)
        with open(self.output_file_name, 'w') as outfile:
            json.dump(data, outfile) 


    def ExecuteBeforeOutputStep(self):
        """ This method is executed right before the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed right after the ouput process computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def ExecuteFinalize(self):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        pass

    def __kratos_vector_to__python_list(self, value):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        list = []
        for index in range(len(value)):
            list.append(value[index])
        return list
    
    def __generate_elemental_gp_mean_value(self, value):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to build the elemental mean of the gauss point results 
        """
        
        list = []
        length = len(value[0])
        for entry in range(length):
            gp_sum = 0.0
            for gp in range(len(value)):
                gp_sum += value[gp][entry]
            gp_mean = gp_sum / len(value)
            list.append(gp_mean)
        return list

    def __generate_total_sub_model_part_mean(self, item):
        """ This method is executed in order to finalize the current computation

        Keyword arguments:
        self -- It signifies an instance of a class.
        item -- The Total List item to geenrate the total mean over the sub_model_part
        """
        
        list = []
        length = len(item[0])
        for index1 in range(length):
            entry_mean = 0.0
            for index2 in range(len(item)):
                entry_mean += item[index2][index1]
            list.append(entry_mean)
        return list

    def __generate_variable_list_from_input(self, param):
        """ Parse a list of variables from input.

        Keyword arguments:
        self -- It signifies an instance of a class.
        value -- The Kratos vector to transform
        """

        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [KratosMultiphysics.KratosGlobals.GetVariable(param[i].GetString()) for i in range( 0, param.size())]

    def __check_flag(self, component):
        """ Checks the flag over a component

        Keyword arguments:
        self -- It signifies an instance of a class.
        component -- The Kratos node or element to check
        """

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False

        return True
