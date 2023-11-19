# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupling_operation import CoSimulationCouplingOperation

# Additional imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

def Create(*args):
    return ComputeBoundaryForce(*args)

class ComputeBoundaryForce(CoSimulationCouplingOperation):
    """This operation is used to compute forces in a boundary, based on the pressure.
    TODO:
    - add messages with different echo-levels
    - add tests
    - more cleanup
    """
    def __init__(self, settings, solver_wrappers, process_info, data_communicator):
        super().__init__(settings, process_info, data_communicator)
        self.model = solver_wrappers[self.settings["solver"].GetString()].model
        self.model_part_name = self.settings["model_part_name"].GetString()
        self.model_part = self.model[self.model_part_name]
        self.write_output_file = self.settings['write_output_file'].GetBool()
        self.format = self.settings["print_format"].GetString()

        self.width = self.settings["width"].GetDouble()
        # If 2D case: width from parameters is used
        # If 3D case: width is not used
        domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
        if domain_size == 3:
            self.width = 1

        self.interval = KM.IntervalUtility(settings)

        if(self.model_part.GetCommunicator().MyPID() == 0):
            if(self.write_output_file):
                output_file_name = self.model_part_name + "_global_force.dat"
                file_handler_settings = KM.Parameters(self.settings["output_file_settings"])
                if file_handler_settings.Has("file_name"):
                    warn_msg = 'Unexpected user-specified entry found in "output_file_settings": {"file_name": '
                    warn_msg += '"' + file_handler_settings["file_name"].GetString() + '"}\n'
                    warn_msg += 'Using this specififed file name instead of the default "' + output_file_name + '"'
                    cs_tools.cs_print_info(self._ClassName(), warn_msg)
                else:
                    file_handler_settings.AddEmptyValue("file_name")
                    file_handler_settings["file_name"].SetString(output_file_name)
                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(self.model_part, file_handler_settings, file_header).file

    def Execute(self):
        current_time = self.model_part.ProcessInfo[KM.TIME]

        if(self.interval.IsInInterval(current_time)):
            results = self._EvaluateGlobalForces()

            if(self.model_part.GetCommunicator().MyPID() == 0):
                output = []
                output.extend(results)
                output_values = [format(val, self.format) for val in output]
                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                output_values.insert(0, str(current_time))

                if(self.echo_level > 2):
                    # print to screen the results at echo level 3 or higher
                    res_labels = ['time: ', 'vel_x: ', 'vel_y: ', 'vel_z: ', 'f_x: ', 'f_y: ', 'f_z: ', 'p: ',]
                    result_msg = 'Boundary Force force evaluation for model part ' + self.model_part_name + '\n'
                    result_msg += ', '.join([a + b for a, b in zip(res_labels, output_values)])
                    cs_tools.cs_print_info(self._ClassName(), result_msg)

                if(self.write_output_file):
                    self.output_file.write(' '.join(output_values) + '\n')

    def _EvaluateGlobalForces(self):
        # vel_x, vel_y, vel_z
        velocity = [0.0, 0.0, 0.0]
        sum_forces = [0.0, 0.0, 0.0]
        pressure_list = [0.0]

        utils = KM.VariableUtils()
        utils.SetVariable(KM.REACTION_X, 0, self.model_part.Nodes)
        utils.SetVariable(KM.REACTION_Y, 0, self.model_part.Nodes)
        utils.SetVariable(KM.REACTION_Z, 0, self.model_part.Nodes)

        for element in self.model_part.Elements:
            geometry = element.GetGeometry()
            nodes = element.GetNodes()
            shape_functions_values = geometry.ShapeFunctionsValues()
            area = geometry.Area()
            unit_normal = geometry.UnitNormal()
            pressure = 0
            c = 0

            for node in nodes:
                pressure_node = node.GetSolutionStepValue(KM.PRESSURE, 0)
                force_node = unit_normal * (-1) * pressure_node * area * shape_functions_values[0, c] * self.width
                pressure += pressure_node * shape_functions_values[0, c]

                node.SetSolutionStepValue(KM.REACTION_X, 0, node.GetSolutionStepValue(KM.REACTION_X, 0) + force_node[0])
                node.SetSolutionStepValue(KM.REACTION_Y, 0, node.GetSolutionStepValue(KM.REACTION_Y, 0) + force_node[1])
                node.SetSolutionStepValue(KM.REACTION_Z, 0, node.GetSolutionStepValue(KM.REACTION_Z, 0) + force_node[2])
                c += 1

            force = unit_normal * pressure * area * self.width

            for i in range(3):
                sum_forces[i] += force[i]

            pressure_list[0] += pressure

        if self.echo_level > 1:
            info_msg = "Computed boundary forces for model part \"" + self.model_part_name  + "\" in solver: \"" + self.settings["solver"].GetString() + "\""
            cs_tools.cs_print_info(self._ClassName(), info_msg)

        return velocity + sum_forces + pressure_list

    def _GetFileHeader(self):
        header = '# Global force for model part ' + self.model_part_name + '\n'
        header += '# Time vel_x vel_y vel_z f_x f_y f_z p\n'
        return header

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "solver"                : "UNSPECIFIED",
            "model_part_name"       : "",
            "interval"              : [0.0, 1e30],
            "print_format"          : ".8f",
            "width"                 : 1.0,
            "write_output_file"     : true,
            "output_file_settings"  : {}
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults
