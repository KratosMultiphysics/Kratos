from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.colors as colors

# Other imports
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput
import sys, time

class DummyFLOWerSolver(object):
    """This class emulates the behavior of an the CFD solver FLOWer
    It can be used for testing and development of couplings without actually having to couple to FLOWer
    """
    def __init__(self, settings):
        default_settings = KM.Parameters("""{
            "name"                        : "UNSPECIFIED",
            "implicit_coupling"           : false,
            "num_steps"                   : 10,
            "io_settings"                 : {},
            "coupling_interfaces"         : [],
            "receive_meshes"              : [],
            "mdpa_file_name"              : "UNSPECIFIED",
            "debugging_settings" : { }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)
        self.settings = settings

        self.model = KM.Model()
        self.model_part = self.model.CreateModelPart("ExtSolver")
        self.name = self.settings["name"].GetString()

        debugging_settings_defaults = KM.Parameters("""{
            "echo_level"         : 0,
            "api_print_timing"   : false,
            "dry_run"            : true,
            "solving_time"       : 0.0,
            "print_colors"       : false
        }""")

        self.settings["debugging_settings"].ValidateAndAssignDefaults(debugging_settings_defaults)
        debugging_settings = self.settings["debugging_settings"]
        self.echo_level = debugging_settings["echo_level"].GetInt()
        self.dry_run = debugging_settings["dry_run"].GetBool()
        self.solving_time = debugging_settings["solving_time"].GetDouble()
        colors.PRINT_COLORS = debugging_settings["print_colors"].GetBool()

        # TODO hardcoded for now, parse it from the variables used in the settings
        self.model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KM.ROTATION)
        self.model_part.AddNodalSolutionStepVariable(KM.FORCE)
        self.model_part.AddNodalSolutionStepVariable(KM.MOMENT)
        self.model_part.AddNodalSolutionStepVariable(KM.VELOCITY)

        severity = KM.Logger.GetDefaultOutput().GetSeverity()
        KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING) # mute MP-IO
        model_part_io = KM.ModelPartIO(self.settings["mdpa_file_name"].GetString())
        model_part_io.ReadModelPart(self.model_part)
        KM.Logger.GetDefaultOutput().SetSeverity(severity)

        # Note: calling "EMPIRE_API_Connect" is NOT necessary, it is replaced by the next two lines
        KratosCoSim.EMPIRE_API.EMPIRE_API_SetEchoLevel(self.echo_level)
        KratosCoSim.EMPIRE_API.EMPIRE_API_PrintTiming(debugging_settings["api_print_timing"].GetBool())

    def Run(self):
        num_coupling_interfaces = self.settings["coupling_interfaces"].size()
        self.__CustomPrint(1, colors.cyan("Starting export") + " of CouplingInterfaces ...")

        for i in range(num_coupling_interfaces):
            coupling_interface_settings = self.settings["coupling_interfaces"][i]
            sub_model_part_name = coupling_interface_settings["sub_model_part_name"].GetString()
            comm_name = coupling_interface_settings["comm_name"].GetString()

            receive_interface = False
            operation = "Exporting"
            if coupling_interface_settings.Has("receive_interface"):
                receive_interface = coupling_interface_settings["receive_interface"].GetBool()
                operation = "Importing"

            self.__CustomPrint(2, colors.cyan(operation)+' coupling-interface "{}" on ModelPart "{}"'.format(comm_name, sub_model_part_name))

            if not self.dry_run:
                if receive_interface:
                    # use an aux-modelpart, which will not be used again
                    KratosCoSim.EMPIRE_API.EMPIRE_API_recvMesh(self.model.CreateModelPart(sub_model_part_name), comm_name)
                    # TODO maybe restructure and do all receives after all sends
                else:
                    KratosCoSim.EMPIRE_API.EMPIRE_API_sendMesh(self.model["ExtSolver."+sub_model_part_name], comm_name)
            else:
                self.__CustomPrint(2, colors.magenta('... skipped'))

            self.__CustomPrint(1, colors.cyan("Finished "+operation) + " of CouplingInterfaces")

        # time loop
        self.__CustomPrint(1, "Starting Solution Loop")
        for i in range(self.settings["num_steps"].GetInt()):
            if self.echo_level > 0: print() # newline
            self.__CustomPrint(1, colors.bold('Step: {}/{}'.format(i+1, self.settings["num_steps"].GetInt())))
            self.SolveSolutionStep()
        self.__CustomPrint(1, "Finished")

    def SolveSolutionStep(self):
        num_coupling_interfaces = self.settings["coupling_interfaces"].size()

        self.__CustomPrint(1, colors.green("Starting import") + " of CouplingInterfaceData ...")

        for i in range(num_coupling_interfaces):
            coupling_interface_settings = self.settings["coupling_interfaces"][i]
            if not coupling_interface_settings.Has("data_field_recv"):
                continue
            sub_model_part_name = coupling_interface_settings["sub_model_part_name"].GetString()
            data_field_settings = coupling_interface_settings["data_field_recv"]

            data_field_name = data_field_settings["data_field_name"].GetString()

            self.__CustomPrint(2, colors.green('Importing')+' data-field "{}" on ModelPart "{}"'.format(data_field_name, sub_model_part_name))

            variables = GenerateVariableListFromInput(data_field_settings["variables"])

            if not self.dry_run:
                KratosCoSim.EMPIRE_API.EMPIRE_API_recvDataField(self.model["ExtSolver."+sub_model_part_name], data_field_name, *variables) # passing all varibales from the list
            else:
                self.__CustomPrint(2, colors.magenta('... skipped'))

        self.__CustomPrint(1, colors.green("Finished import") + " of CouplingInterfaceData")

        if self.echo_level > 0: print() # newline
        self.__CustomPrint(1, colors.blue("Solving ..."))
        time.sleep(self.solving_time)
        self.__CustomPrint(2, "Solving took {} [sec]".format(self.solving_time))
        if self.echo_level > 0: print() # newline
        # TODO implement random values ... ?

        self.__CustomPrint(1, colors.cyan("Starting export") + " of CouplingInterfaceData ...")

        for i in range(num_coupling_interfaces):
            coupling_interface_settings = self.settings["coupling_interfaces"][i]
            if not coupling_interface_settings.Has("data_field_send"):
                continue
            sub_model_part_name = coupling_interface_settings["sub_model_part_name"].GetString()
            data_field_settings = coupling_interface_settings["data_field_send"]

            data_field_name = data_field_settings["data_field_name"].GetString()

            self.__CustomPrint(2, colors.cyan('Exporting')+' data-field "{}" on ModelPart "{}"'.format(data_field_name, sub_model_part_name))

            variables = GenerateVariableListFromInput(data_field_settings["variables"])

            if not self.dry_run:
                KratosCoSim.EMPIRE_API.EMPIRE_API_sendDataField(self.model["ExtSolver."+sub_model_part_name], data_field_name, *variables) # passing all varibales from the list
            else:
                self.__CustomPrint(2, colors.magenta('... skipped'))

        self.__CustomPrint(1, colors.cyan("Finished export") + " of CouplingInterfaceData")

    def __CustomPrint(self, echo_level, *args):
        if self.echo_level >= echo_level:
            KM.Logger.PrintInfo(self.name, (echo_level-1)*"   " + " ".join(map(str,args)))

if __name__ == '__main__':
    if len(sys.argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python external_solver.py <flower-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = sys.argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    simulation = DummyFLOWerSolver(parameters)
    simulation.Run()