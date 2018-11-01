from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics
import json
import math
from json_utilities import *

# Importing the base class
from from_json_check_result_process import FromJsonCheckResultProcess

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ParticleMPMFromJsonCheckResultProcess(Model, settings["Parameters"])

class ParticleMPMFromJsonCheckResultProcess(FromJsonCheckResultProcess, KratosUnittest.TestCase):

    def __init__(self, model_part, params):
        super(ParticleMPMFromJsonCheckResultProcess, self).__init__(model_part, params)

    def ExecuteFinalizeSolutionStep(self):

        tol = self.params["tolerance"].GetDouble()
        reltol = self.params["relative_tolerance"].GetDouble()
        time = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.TIME)
        dt = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.DELTA_TIME)
        step = self.sub_model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)
        self.time_counter += dt
        if self.time_counter > self.frequency:
            self.time_counter = 0.0
            input_time_list = self.data["TIME"]

            # Material points values
            for mp in self.sub_model_part.Elements:
                compute = self.__check_flag(mp)

                if (compute == True):
                    for i in range(self.params["check_variables"].size()):
                        out = self.params["check_variables"][i]
                        variable_name = out.GetString()
                        variable = KratosMultiphysics.KratosGlobals.GetVariable( variable_name )
                        variable_type = KratosMultiphysics.KratosGlobals.GetVariableType(variable_name)
                        value = mp.GetValue(variable)

                        # Scalar variable
                        if (variable_type == "Double" or variable_type == "Component"):
                            values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name]
                            value_json = self.__linear_interpolation(time, input_time_list, values_json)
                            if (self.iscloseavailable == True):
                                isclosethis = math.isclose(value, value_json, rel_tol=reltol, abs_tol=tol)
                                self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking mp " + str(mp.Id) + " " + variable_name + " results."))
                            else:
                                self.assertAlmostEqual(value, value_json, msg=("Error checking mp " + str(mp.Id) + " " + variable_name + " results."), delta=tol)
                        # Array variable
                        elif variable_type == "Array":
                            if (KratosMultiphysics.KratosGlobals.GetVariableType(variable_name + "_X") == "Component"):
                                if (self.iscloseavailable == True):
                                    # X-component
                                    values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name + "_X"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = math.isclose(value[0], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking mp " + str(mp.Id) + " " + variable_name + " results."))
                                    # Y-component
                                    values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name + "_Y"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = math.isclose(value[1], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != "+str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking mp " + str(mp.Id) + " " + variable_name + " results."))
                                    # Z-component
                                    values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name + "_Z"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    isclosethis = math.isclose(value[2], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value)+" != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking mp " + str(mp.Id) + " " + variable_name + " results."))
                                else:
                                    # X-component
                                    values_json = self.data["PARTICLE_"+str(mp.Id)][variable_name  + "_X"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    self.assertAlmostEqual(value[0], value_json, msg=("Error checking mp " + str(mp.Id) + " " + variable_name + " X-component results."), delta=tol)
                                    # Y-component
                                    values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name + "_Y"]
                                    value_json = self.__linear_interpolation(time, input_time_list, values_json)
                                    self.assertAlmostEqual(value[1], value_json, msg=("Error checking mp " + str(mp.Id) + " " + variable_name + " Y-component results."), delta=tol)
                                    # Z-component
                                    values_json = self.data["PARTICLE_" + str(mp.Id)][variable_name + "_Z"]
                                    value_json = self._linear_interpolation(time, input_time_list, values_json)
                                    self.assertAlmostEqual(value[2], value_json, msg=("Error checking mp " + str(mp.Id) + " " + variable_name + " Z-component results."), delta=tol)
                            else:
                                values_json = self.data["PARTICLE_"+str(mp.Id)][variable_name][step - 1]
                                for index in range(len(value)):
                                    value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                    if (self.iscloseavailable == True):
                                        isclosethis = math.isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                        self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking mp " + str(mp.Id) + " " + variable_name + " results."))
                                    else:
                                        self.assertAlmostEqual(value[index], value_json, msg=("Error checking mp " + str(mp.Id) + " " + variable_name + " results."), delta=tol)
                        # Vector variable
                        elif variable_type == "Vector":
                            values_json = self.data["PARTICLE_"+str(mp.Id)][variable_name][step - 1]
                            for index in range(len(value)):
                                value_json = values_json[index] # self.__linear_interpolation(time, input_time_list, values_json[index])
                                if (self.iscloseavailable == True):
                                    isclosethis = math.isclose(value[index], value_json, rel_tol=reltol, abs_tol=tol)
                                    self.assertTrue(isclosethis, msg=(str(value) + " != " + str(value_json) + ", rel_tol = " + str(reltol) + ", abs_tol = " + str(tol) + " : Error checking mp " + str(mp.Id) + " " + variable_name + " results."))
                                else:
                                    self.assertAlmostEqual(value[index], value_json, msg=("Error checking mp " + str(mp.Id) + " " + variable_name + " results."), delta=tol)

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
