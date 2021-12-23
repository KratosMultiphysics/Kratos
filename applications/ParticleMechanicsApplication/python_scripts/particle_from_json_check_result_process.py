from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.from_json_check_result_process import LegacyFromJsonCheckResultProcess # TODO: This must be updated to the new C++ version
from KratosMultiphysics.KratosUnittest import isclose as t_isclose

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ParticleFromJsonCheckResultProcess(Model, settings["Parameters"])

class ParticleFromJsonCheckResultProcess(LegacyFromJsonCheckResultProcess, KratosUnittest.TestCase):

    def __init__(self, model_part, params):
        super(ParticleFromJsonCheckResultProcess, self).__init__(model_part, params)

    def ExecuteFinalizeSolutionStep(self):

        tol = self.params["tolerance"].GetDouble()
        reltol = self.params["relative_tolerance"].GetDouble()
        time = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            input_time_list = self.data["TIME"]

            # Material points values
            for mp in self.model_part.Elements:
                compute = self.__check_flag(mp)

                if (compute == True):
                    for i in range(self.params["check_variables"].size()):
                        out = self.params["check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        value = mp.CalculateOnIntegrationPoints(variable,self.model_part.ProcessInfo)[0]

                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Integer" or variable_type == "Component"):
                            values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name]
                            value_json = self.__linear_interpolation(time, input_time_list, values_json)
                            isclosethis = t_isclose(value, value_json, rel_tol=reltol, abs_tol=tol)
                            self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking particle " + str(mp.Id) + " " + variable_name + " results."))
                        # Array variable
                        elif variable_type == "Array":

                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Double"):
                                for component_index, component in enumerate(["_X", "_Y", "_Z"]):
                                    values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name +component]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = t_isclose(value[component_index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value[component_index]) + " != "+str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking particle " + str(mp.Id) + " " + variable_name + " results."))
                            else:
                                values_json = self.data["PARTICLE_"+str(mp.Id)][variable_name][step - 1]
                                for index in range(len(value)):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    isclosethis = t_isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking particle " + str(mp.Id) + " " + variable_name + " results."))
                        # Vector variable
                        elif variable_type == "Vector":
                            values_json = self.data["PARTICLE_"+str(mp.Id)][variable_name][step - 1]
                            for index in range(len(value)):
                                value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                isclosethis = t_isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking particle " + str(mp.Id) + " " + variable_name + " results."))

    def __linear_interpolation(self, x, x_list, y_list):

        tb = KratosMultiphysics.PiecewiseLinearTable()
        for i in range(len(x_list)):
            tb.AddRow(x_list[i], y_list[i])

        return tb.GetNearestValue(x)

    def __generate_variable_list_from_input(self,param):

        # At least verify that the input is a string
        if not param.IsArray():
            raise Exception("{0} Error: Variable list is unreadable".format(self.__class__.__name__))

        # Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
        return [ KratosMultiphysics.KratosGlobals.GetVariable( param[i].GetString() ) for i in range( 0,param.size() ) ]

    def __check_flag(self, component):

        if self.flag != None:
            if component.Is(self.flag) == False:
                return False

        return True
