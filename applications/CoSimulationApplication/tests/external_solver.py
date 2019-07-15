from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis
import KratosMultiphysics.CoSimulationApplication as KratosCoSim

import os, sys, time

class ExternalSolver(object):
    """This class emulates the behavior of an external solver, hence it can be used for testing
    TFor now it uses directly the EMPIRE_API
    In the future this will be replaced by the CoSimIO (which is included in the external solvers)
    """
    def __init__(self, settings):
        default_settings = KM.Parameters("""{
            "name"                           : "UNSPECIFIED",
            "do_coupling"                    : true,
            "implicit_coupling"              : false,
            "num_steps"                      : 10,
            "coupling_interfaces"            : [],
            "receive_meshes"                 : [],
            "mdpa_file_name"                 : "UNSPECIFIED",
            "api_configuration_file_name"    : "UNSPECIFIED",
            "echo_level"                     : 0,
            "solving_time"                   : 1
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.settings = settings

        self.do_coupling = self.settings["do_coupling"].GetBool()

        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("ExtSolver")
        self.name = self.settings["name"].GetString()
        self.echo_level = self.settings["echo_level"].GetInt()
        self.solving_time = self.settings["solving_time"].GetInt()

        # TODO hardcoded for now, parse it from the variables used in the settings
        self.model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        self.model_part.AddNodalSolutionStepVariable(KM.FORCE)
        self.model_part.AddNodalSolutionStepVariable(KM.MOMENT)

        severity = KM.Logger.GetDefaultOutput().GetSeverity()
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING) # mute MP-IO
        model_part_io = KM.ModelPartIO(self.settings["mdpa_file_name"].GetString())
        model_part_io.ReadModelPart(self.model_part)
        KM.Logger.GetDefaultOutput().SetSeverity(severity)

        KratosCoSim.EMPIRE_API.EMPIRE_API_Connect(self.settings["api_configuration_file_name"].GetString())

    def Run(self):
        num_coupling_interfaces = self.settings["coupling_interfaces"].size()
        if self.do_coupling:
            self.__CustomPrint(1, "Starting export of CouplingInterfaces ...")

            for i in range(num_coupling_interfaces):
                sub_model_part_name = self.settings["coupling_interfaces"][i]["sub_model_part_name"].GetString()
                comm_name = self.settings["coupling_interfaces"][i]["comm_name"].GetString()
                KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(self.model["ExtSolver."+sub_model_part_name], comm_name)

            self.__CustomPrint(1, "Sucessfully finished export of CouplingInterfaces")

        if self.do_coupling:
            if self.settings["receive_meshes"].size() > 0:
                self.__CustomPrint(1, "Starting import of CouplingInterfaces ...")

            for mesh_name in self.settings["receive_meshes"].GetStringArray():
                stopNotImplemented
                model_part_for_receiving = self.model.CreateModelPart(mesh_name)
                KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(model_part_for_receiving)

            if self.settings["receive_meshes"].size() > 0:
                self.__CustomPrint(1, "Sucessfully finished import of CouplingInterfaces")

        # time loop
        self.__CustomPrint(1, "Starting Solution Loop")
        for i in range(self.settings["num_steps"].GetInt()):
            print() # newline
            self.__CustomPrint(1, 'Step: {}/{}'.format(i+1, self.settings["num_steps"].GetInt()))
            self.SolveSolutionStep()
        self.__CustomPrint(1, "Finished")

    def SolveSolutionStep(self):
        num_coupling_interfaces = self.settings["coupling_interfaces"].size()

        if self.do_coupling:
            self.__CustomPrint(1, "Starting import of CouplingInterfaceData ...")

            for i in range(num_coupling_interfaces):
                coupling_interface_settings = self.settings["coupling_interfaces"][i]
                sub_model_part_name = coupling_interface_settings["sub_model_part_name"].GetString()
                data_field_settings = coupling_interface_settings["data_field_send"]

                data_field_name = data_field_settings["data_field_name"].GetString()
                variables = GenerateVariableListFromInput(data_field_settings["variables"])
                KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(self.model["ExtSolver."+sub_model_part_name], data_field_name, *variables)

            self.__CustomPrint(1, "Sucessfully finished import of CouplingInterfaceData")

        self.__CustomPrint(1, "Solving ...")
        time.sleep(self.solving_time)
        self.__CustomPrint(2, "  Solving took {} [sec]".format(self.solving_time))
        # TODO implement random values ... ?

        if self.do_coupling:
            self.__CustomPrint(1, "Starting export of CouplingInterfaceData ...")

            for i in range(num_coupling_interfaces):
                sub_model_part_name = coupling_interface_settings["sub_model_part_name"].GetString()
                data_field_settings = coupling_interface_settings["data_field_recv"]

                data_field_name = data_field_settings["data_field_name"].GetString()
                variables = GenerateVariableListFromInput(data_field_settings["variables"])
                KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(self.model["ExtSolver."+sub_model_part_name], data_field_name, *variables)

            self.__CustomPrint(1, "Sucessfully finished export of CouplingInterfaceData")

    def __CustomPrint(self, echo_level, *args):
        if self.echo_level >= echo_level:
            KM.Logger.PrintInfo(self.name, " ".join(map(str,args)))

def GenerateVariableListFromInput(param):
    """ Retrieve variable name from input (a string) and request the corresponding C++ object to the kernel
    """
    from KratosMultiphysics import KratosGlobals
    return [ KratosGlobals.GetVariable( var_name ) for var_name in param.GetStringArray() ]

if __name__ == '__main__':
    if len(sys.argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python external_solver.py <flower-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = sys.argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    simulation = ExternalSolver(parameters)
    simulation.Run()