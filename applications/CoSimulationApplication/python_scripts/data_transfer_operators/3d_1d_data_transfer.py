# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools
import KratosMultiphysics.CoSimulationApplication.colors as colors

def Create(*args):
    return Kratos3D1DDataTransferOperator(*args)

class Kratos3D1DDataTransferOperator(CoSimulationDataTransferOperator):
    """DataTransferOperator that transfers values from one 3D interface (ModelPart) to another 1D interface (ModelPart).
    The data_transfer_3d_1ds of the Kratos-MappingApplication are used
    """

    def __init__(self, settings, parent_coupled_solver_data_communicator):
        if KM.IsDistributedRun():
            raise Exception("This function can only be called when Kratos is running in MPI!")
        if not settings.Has("3d_1d_data_transfer_settings"):
            raise Exception('No "3d_1d_data_transfer_settings" provided!')
        super().__init__(settings, parent_coupled_solver_data_communicator)
        self.origin_is_3d = None
        self.__data_transfer_process = {}

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        model_part_origin_name = from_solver_data.model_part_name
        variable_origin        = from_solver_data.variable
        identifier_origin      = from_solver_data.solver_name + "." + model_part_origin_name + "." + variable_origin.Name()

        model_part_destination_name = to_solver_data.model_part_name
        variable_destination        = to_solver_data.variable
        identifier_destination      = to_solver_data.solver_name + "." + model_part_destination_name + "." + variable_destination.Name()

        identifier_tuple         = (identifier_origin, identifier_destination)
        inverse_identifier_tuple = (identifier_destination, identifier_origin)

        if identifier_tuple in self.__data_transfer_process:
            if self.origin_is_3d == True:
                self.__data_transfer_process[identifier_tuple].Set(KM.Mapper.USE_TRANSPOSE)
            self.__data_transfer_process[identifier_tuple].Execute()
            if self.origin_is_3d == True:
                self.__data_transfer_process[identifier_tuple].Reset(KM.Mapper.USE_TRANSPOSE)
        elif inverse_identifier_tuple in self.__data_transfer_process:
            if self.origin_is_3d == False:
                self.__data_transfer_process[inverse_identifier_tuple].Set(KM.Mapper.USE_TRANSPOSE)
            self.__data_transfer_process[inverse_identifier_tuple].Execute()
            if self.origin_is_3d == False:
                self.__data_transfer_process[inverse_identifier_tuple].Reset(KM.Mapper.USE_TRANSPOSE)
        else:
            model_part_origin      = from_solver_data.GetModelPart()
            model_part_destination = to_solver_data.GetModelPart()

            if self.echo_level > 0:
                info_msg  = "Creating Mapper:\n"
                info_msg += '    Origin: ModelPart "{}" of solver "{}"\n'.format(model_part_origin_name, from_solver_data.solver_name)
                info_msg += '    Destination: ModelPart "{}" of solver "{}"'.format(model_part_destination_name, to_solver_data.solver_name)

                cs_tools.cs_print_info(colors.bold(self._ClassName()), info_msg)

            # Creating parameters for the data transfer
            # TODO: I should check if the parameters change between different data transfers
            parameters = KM.Parameters(self.settings["3d_1d_data_transfer_settings"].WriteJsonString())
            if not parameters.Has("origin_variables"):
                parameters.AddEmptyArray("origin_variables")
            parameters["origin_variables"].Append(variable_origin.Name())
            if not parameters.Has("destination_variables"):
                parameters.AddEmptyArray("destination_variables")
            parameters["destination_variables"].Append(variable_destination.Name())
            for transfer_option in transfer_options.GetStringArray():
                if transfer_option == "swap_sign":
                    parameters["swap_sign"].SetBool(True)

            # Check if the origin model part is 3D or 1D
            self.origin_is_3d = self.__check_model_part_3D(model_part_origin)
            # Clone is necessary because the settings are validated and defaults assigned, which could influence the creation of other data transfers
            self.__data_transfer_process[identifier_tuple] = KratosCoSim.DataTransfer3D1DProcess(model_part_origin, model_part_destination, parameters.Clone())
            self.__data_transfer_process[inverse_identifier_tuple] = KratosCoSim.DataTransfer3D1DProcess(model_part_origin, model_part_destination, parameters.Clone())

            # Execute the data transfer
            if self.origin_is_3d:
                self.__data_transfer_process[identifier_tuple].Set(KM.Mapper.USE_TRANSPOSE)
            self.__data_transfer_process[identifier_tuple].Execute()
            if self.origin_is_3d:
                self.__data_transfer_process[identifier_tuple].Reset(KM.Mapper.USE_TRANSPOSE)

    def _Check(self, from_solver_data, to_solver_data):
        def CheckData(data_to_check):
            if "node" not in data_to_check.location:
                raise Exception('Transfer only supports nodal values!"{}"\nChecking ModelPart "{}" of solver "{}"'.format(self._ClassName(), data_to_check.model_part_name, data_to_check.solver_name))

        CheckData(from_solver_data)
        CheckData(to_solver_data)

    def __check_model_part_3D(self, model_part):
        is_1d = False
        for elem in model_part.Elements:
            geom = elem.GetGeometry()
            if geom.LocalSpaceDimension() == 1:
                is_1d = True
            break
        return not is_1d

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "3d_1d_data_transfer_settings" : {
                "origin_variables"         : [],
                "destination_variables"    : [],
                "swap_sign"                : false,
                "interpolate_parameters"   : {
                    "data_transfer_3d_1d_type" : "nearest_element",
                    "echo_level"  : 0,
                    "search_settings" : {
                        "max_num_search_iterations"     : 8,
                        "echo_level"                    : 0
                    }
                },
                "search_parameters"        :  {
                    "allocation_size"         : 100,
                    "bucket_size"             : 4,
                    "search_factor"           : 2.0,
                    "search_increment_factor" : 1.5
                }
            }
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return ["swap_sign"]
