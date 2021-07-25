import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility


def Factory(params, model):
    if not isinstance(params, KM.Parameters):
        raise Exception('expected input shall be a Parameters object, encapsulating a json string')
    return ComputeBoundaryForceProcess(model, params['Parameters'])


class ComputeBoundaryForceProcess(KM.Process):
    '''
    Get the external accelerations and computes the hydrostatic forces.
    The results are written in a file or printed into screen.
    '''

    def __init__(self, model, params):
        '''Constructor of ComputeBoundaryForceProcess.'''

        super().__init__()

        default_settings = KM.Parameters("""
            {
                "model_part_wall_name"   : "",
                "model_part_bottom_name" : "",
                "interval"               : [0.0, 1e30],                
                "print_to_screen"        : false,
                "print_format"           : ".8f",
                "write_output_file"      : true,
                "output_file_settings"   : {}
            }
            """)

        self.interval = KM.IntervalUtility(params)

        params.ValidateAndAssignDefaults(default_settings)

        self.model_part_wall_name = params['model_part_wall_name'].GetString()
        self.model_part_wall = model[self.model_part_wall_name]
        self.model_part_bottom_name = params['model_part_bottom_name'].GetString()
        self.model_part_bottom = model[self.model_part_bottom_name]

        self.print_to_screen = params['print_to_screen'].GetBool()
        self.write_output_file = params['write_output_file'].GetBool()
        self.print_format = params["print_format"].GetString()

        if (self.model_part_wall.GetCommunicator().MyPID() == 0):
            if (self.write_output_file):

                default_file_name = params["model_part_wall_name"].GetString() + "_global_force.dat"

                file_handler_params = KM.Parameters(
                    params["output_file_settings"])

                if file_handler_params.Has("file_name"):
                    file_name = file_handler_params["file_name"].GetString()
                    msg = 'Unexpected user-specified entry found in "output_file_settings":\n'
                    msg += '\t"file_name" : "{}"\n'
                    msg += 'Using this specified file name instead of the default ("{}")'
                    KM.Logger.PrintWarning("ComputeBoundaryForceProcess", msg.format(file_name, default_file_name))
                else:
                    file_handler_params.AddString("file_name", default_file_name)

                file_header = self._GetFileHeader()
                self.output_file = TimeBasedAsciiFileWriterUtility(
                    self.model_part_wall, file_handler_params, file_header).file


    def ExecuteFinalizeSolutionStep(self):
        '''Print the boundary forces to a file or into screen.'''

        current_time = self.model_part_wall.ProcessInfo[KM.TIME]

        if self.interval.IsInInterval(current_time):
            accelerations, forces = self._EvaluateGlobalForces()

            if self.model_part_wall.GetCommunicator().MyPID() == 0:

                output_values = []
                # not formatting time in order to not lead to problems with time recognition
                # in the file writer when restarting
                output_values.append(str(current_time))
                for val in accelerations : output_values.append(format(val, self.print_format))
                for val in forces : output_values.append(format(val, self.print_format))

                if self.print_to_screen:
                    result_msg = 'Global force evaluation for model part ' + \
                        self.model_part_wall_name + '\n'
                    res_labels = ['time: ','acc_x: ','acc_y: ','acc_z: ','f_x: ', 'f_y: ', 'f_z: ']
                    result_msg += ', '.join([a+b for a,b in zip(res_labels, output_values)])
                    self._PrintToScreen(result_msg)

                if self.write_output_file:
                    self.output_file.write(' '.join(output_values) + '\n')


    def _EvaluateGlobalForces(self):
        for node in self.model_part_wall.Nodes:
            acceleration = node.GetSolutionStepValue(KM.MESH_ACCELERATION)
            break

        process_info = self.model_part_wall.ProcessInfo
        sum_forces = SW.ShallowWaterUtilities().ComputeHydrostaticForces(self.model_part_wall.Conditions, process_info)
        sum_forces += SW.ShallowWaterUtilities().ComputeHydrostaticForces(self.model_part_bottom.Elements, process_info)

        return acceleration, sum_forces


    def _GetFileHeader(self):
        header = '# Global force for model part ' + self.model_part_wall_name + '\n'
        header += '# Time acc_x acc_y acc_z f_x f_y f_z\n'
        return header


    @staticmethod
    def _PrintToScreen(result_msg):
        KM.Logger.PrintInfo('ComputeBoundaryForceProcess', 'Global force results - flow- and body-attached:')
        KM.Logger.PrintInfo('ComputeBoundaryForceProcess', 'Current time: ' + result_msg)
