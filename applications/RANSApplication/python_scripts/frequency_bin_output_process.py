# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS
# other imports
from KratosMultiphysics.RANSApplication import RansAuxiliaryUtilities
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, Model):
    if(type(settings) != Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return FrequencyBinOutputProcess(Model, settings["Parameters"])

class FrequencyBinOutputProcess(Kratos.OutputProcess):
    def __init__(self, model, parameters):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters('''{
            "model_part_name"           : "PLEASE_SPECIFY_MAIN_MODEL_PART_NAME",
            "point_coordinates"         : [0.0, 0.0, 0.0],
            "velocity_direction"        : [0.0, 1.0, 0.0],
            "frequency_bin_indices"     : [1, 2, 3, 4, 5],
            "total_number_of_time_steps": 1000,
            "echo_level"                : 0,
            "output_file_settings"      : {}
        }''')

        self.model = model
        self.parameters = parameters

        if (self.parameters.Has("frequency_bin_indices") and self.parameters["frequency_bin_indices"].IsString()):
            if (self.parameters["frequency_bin_indices"].GetString() == "ALL_FREQUENCY_BIN_INDICES"):
                if (self.parameters.Has("total_number_of_time_steps")):
                    total_steps = self.parameters["total_number_of_time_steps"].GetInt()
                else:
                    total_steps = 1000

                v = Kratos.Vector(total_steps // 2)
                for i in range(total_steps // 2):
                    v[i] = i

                self.parameters["frequency_bin_indices"].SetVector(v)

        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.main_model_part = self.model.GetModelPart(self.parameters["model_part_name"].GetString())

        self.output_file = None
        self.element_id = -1
        self.point_shape_function_values = None

        point_coordinates = self.parameters["point_coordinates"].GetVector()
        if (point_coordinates.Size() != 3):
            raise Exception("Point coordinates must be of a vector with size 3. [ point_coordinates = " + str(point_coordinates) + " ].")
        self.point_coordinates = Kratos.Point(point_coordinates[0], point_coordinates[1], point_coordinates[2])

        velocity_direction = self.parameters["velocity_direction"].GetVector()
        if (velocity_direction.Size() != 3):
            raise Exception("Velocity direction must be of a vector with size 3. [ velocity_direction = " + str(velocity_direction) + " ].")
        self.velocity_direction = Kratos.Array3([velocity_direction[0], velocity_direction[1], velocity_direction[2]])

        self.total_number_of_steps = self.parameters["total_number_of_time_steps"].GetInt()

        self.frequency_bin_indices_list = []
        for v in self.parameters["frequency_bin_indices"].GetVector():
            self.frequency_bin_indices_list.append(int(v))

    def ExecuteInitialize(self):
        self.element_id, self.point_shape_function_values = self._GetElementId()

        # Only rank 0 writes in MPI
        my_rank = 0
        comm = self.main_model_part.GetCommunicator().GetDataCommunicator()
        self.is_writing_rank = my_rank == comm.Rank()
        if self.is_writing_rank:
            file_handler_params = Kratos.Parameters(self.parameters["output_file_settings"])
            file_header = self.GetFileHeader()
            self.output_file =  TimeBasedAsciiFileWriterUtility(self.main_model_part, file_handler_params, file_header).file

    def PrintOutput(self):
        out = str(self.main_model_part.ProcessInfo[Kratos.TIME])

        real_values, imag_values = self._ComputeFrequencyBinComponents()
        for real, imag in zip(real_values, imag_values):
            out += "," + str(real) + "," + str(imag)

        out += "\n"

        if self.is_writing_rank:
            self.output_file.write(out)

    def ExecuteFinalize(self):
        if self.is_writing_rank:
            self.output_file.close()

    def GetFileHeader(self):
        header  = '# Frequency bin results ' + '\n'
        header += '# time'

        for frequency_bin_index in self.frequency_bin_indices:
            header += ", Bin {0:d}-real, Bin {0:d}-imag".format(frequency_bin_index)

        header += "\n"

        return header

    def _ComputeFrequencyBinComponents(self):
        real_values = Kratos.Vector()
        imag_values = Kratos.Vector()

        RansAuxiliaryUtilities.CalculateFrequencyBinValues(
            real_values,
            imag_values,
            self._GetTimeInstantaneousValue(),
            self.main_model_part.ProcessInfo[Kratos.STEP],
            self.total_number_of_steps,
            self.frequency_bin_indices_list)

        return real_values, imag_values

    def _GetTimeInstantaneousValue(self):
        if (self.element_id == -1):
            raise Exception("No element found for give point. [ Point coordinates = " + str(self.point_coordinates) + " ].")

        current_velocity = Kratos.Array3(0.0)
        for index, node in enumerate(self.main_model_part.GetElement(self.element_id).GetGeometry()):
            current_velocity += node.GetSolutionStepValue(Kratos.VELOCITY) * self.point_shape_function_values[index]

        time_instantaneous_value = current_velocity[0] * self.velocity_direction[0] + current_velocity[1] * self.velocity_direction[1] + current_velocity[2] * self.velocity_direction[2]

        return time_instantaneous_value

    def _GetElementData(self):
        point_shape_function_values = Kratos.Vector()
        element_id = Kratos.BruteForcePointLocator(self.main_model_part).FindElement(
            self.point_coordinates,
            point_shape_function_values,
            Kratos.Configuration.Current,
            1e-6)

        return element_id, point_shape_function_values
