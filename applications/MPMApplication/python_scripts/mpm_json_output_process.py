# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.json_utilities import read_external_json, write_external_json

# Importing the base class
from KratosMultiphysics.json_output_process import JsonOutputProcess

def Factory(settings, Model):
    if not isinstance(Model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMJsonOutputProcess(Model, settings["Parameters"])

class MPMJsonOutputProcess(JsonOutputProcess):

    def ExecuteBeforeSolutionLoop(self):

        data = {}
        data["TIME"] = []
        count = 0

        # Material point elements values
        for mp in self.sub_model_part.Elements:
            compute = self.__check_flag(mp)

            if compute:
                mp_identifier = "MP_" + str(mp.Id)

                if not self.resultant_solution:
                    data[mp_identifier] = {}
                else:
                    if count == 0:
                        data["RESULTANT"] = {}

                for i in range(self.params["gauss_points_output_variables"].size()):
                    out = self.params["gauss_points_output_variables"][i]
                    variable_name = out.GetString()
                    variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)

                    if variable_type == "Double" or variable_type == "Integer" or variable_type == "Component":
                        if not self.resultant_solution:
                            data[mp_identifier][variable_name] = []
                        else:
                            if count == 0:
                                data["RESULTANT"][variable_name] = []
                    elif variable_type == "Array":
                        if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                            if not self.resultant_solution:
                                data[mp_identifier][variable_name + "_X"] = []
                                data[mp_identifier][variable_name + "_Y"] = []
                                data[mp_identifier][variable_name + "_Z"] = []
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name + "_X"] = []
                                    data["RESULTANT"][variable_name + "_Y"] = []
                                    data["RESULTANT"][variable_name + "_Z"] = []
                        else:
                            if not self.resultant_solution:
                                data[mp_identifier][variable_name] = []
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name] = []
                    elif variable_type == "Vector":
                        if not self.resultant_solution:
                            data[mp_identifier][variable_name] = []
                        else:
                            if count == 0:
                                data["RESULTANT"][variable_name] = []

                count += 1

        write_external_json(self.output_file_name, data)

    def ExecuteFinalizeSolutionStep(self):

        data = read_external_json(self.output_file_name)

        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            data["TIME"].append(time)
            count = 0

            # Material points values
            for mp in self.sub_model_part.Elements:
                compute = self.__check_flag(mp)

                if compute:
                    mp_identifier = "MP_" + str(mp.Id)

                    for i in range(self.params["gauss_points_output_variables"].size()):
                        out = self.params["gauss_points_output_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        values_vector = mp.CalculateOnIntegrationPoints(variable, self.sub_model_part.ProcessInfo)
                        value = values_vector[0]

                        if variable_type == "Double" or variable_type == "Integer" or variable_type == "Component":
                            if not self.resultant_solution:
                                data[mp_identifier][variable_name].append(value)
                            else:
                                if count == 0:
                                    data["RESULTANT"][variable_name].append(value)
                                else:
                                    data["RESULTANT"][variable_name][-1] += value
                        elif variable_type == "Array":
                            if KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double":
                                if not self.resultant_solution:
                                    data[mp_identifier][variable_name + "_X"].append(value[0])
                                    data[mp_identifier][variable_name + "_Y"].append(value[1])
                                    data[mp_identifier][variable_name + "_Z"].append(value[2])
                                else:
                                    if count == 0:
                                        data["RESULTANT"][variable_name + "_X"].append(value[0])
                                        data["RESULTANT"][variable_name + "_Y"].append(value[1])
                                        data["RESULTANT"][variable_name + "_Z"].append(value[2])
                                    else:
                                        data["RESULTANT"][variable_name + "_X"][-1] += value[0]
                                        data["RESULTANT"][variable_name + "_Y"][-1] += value[1]
                                        data["RESULTANT"][variable_name + "_Z"][-1] += value[2]
                            else:
                                if not self.resultant_solution:
                                    list = self.__kratos_vector_to__python_list(value)
                                    data[mp_identifier][variable_name].append(list)
                                else:
                                    aux = 0.0
                                    for index in range(len(value)):
                                        aux += value[index]
                                    if count == 0:
                                        data["RESULTANT"][variable_name].append(aux)
                                    else:
                                        data["RESULTANT"][variable_name][-1] += aux
                        elif variable_type == "Vector":
                            if not self.resultant_solution:
                                list = self.__kratos_vector_to__python_list(value)
                                data[mp_identifier][variable_name].append(list)
                            else:
                                if count == 0:
                                    list = self.__kratos_vector_to__python_list(value)
                                    data["RESULTANT"][variable_name][-1] += list

                count += 1

        write_external_json(self.output_file_name, data)

    def __kratos_vector_to__python_list(self, value):

        list = []
        for index in range(len(value)):
            list.append(value[index])
        return list

    def __check_flag(self, component):

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False
        return True
