# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_data_transfer_operator import CoSimulationDataTransferOperator

# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication as KSM
import KratosMultiphysics.MappingApplication as KratosMapping

def Create(*args):
    return KratosBeamMapping(*args)

class KratosBeamMapping(CoSimulationDataTransferOperator):

    # currently available mapper-flags aka transfer-options
    __mapper_flags_dict = {
        "swap_sign"     : KratosMapping.Mapper.SWAP_SIGN
    }

    def __init__(self, settings, parent_coupled_solver_data_communicator):
        super(KratosBeamMapping, self).__init__(settings, parent_coupled_solver_data_communicator)
        self.beam_mapper = None
        self.model_part_name_beam = self.settings["model_part_name_beam"].GetString()
        self.model_part_name_surface = self.settings["model_part_name_surface"].GetString()
        self.solver_name_beam = self.settings["solver_name_beam"].GetString()
        self.solver_name_surface = self.settings["solver_name_surface"].GetString()

        self.second_variable_displacement = KM.KratosGlobals.GetVariable(self.settings["second_variable_displacement"].GetString())
        self.second_variable_force = KM.KratosGlobals.GetVariable(self.settings["second_variable_force"].GetString())

    def _ExecuteTransferData(self, from_solver_data, to_solver_data, transfer_options):
        # TODO check location of data => should coincide with the one for the mapper
        # or throw if it is not in a suitable location (e.g. on the ProcessInfo)

        self._CheckAvailabilityTransferOptions(transfer_options)
        
        from_solver_name = from_solver_data.solver_name
        to_solver_name = to_solver_data.solver_name

        from_model_part_name = from_solver_data.model_part_name
        to_model_part_name = to_solver_data.model_part_name

        if self.solver_name_beam == from_solver_name and self.solver_name_surface == to_solver_name:
            inverse_map = False

            if from_model_part_name != self.model_part_name_beam:
                raise Exception(
                    f"Invalid 'from_model_part_name' for beam→surface mapping.\n"
                    f"Expected: '{self.model_part_name_beam}', but got: '{from_model_part_name}'."
                )

            if to_model_part_name != self.model_part_name_surface:
                raise Exception(
                    f"Invalid 'to_model_part_name' for beam→surface mapping.\n"
                    f"Expected: '{self.model_part_name_surface}', but got: '{to_model_part_name}'."
                )

            
        elif self.solver_name_surface == from_solver_name and self.solver_name_beam == to_solver_name:
            inverse_map = True

            if to_model_part_name != self.model_part_name_beam:
                raise Exception(
                    f"Invalid 'to_model_part_name' for surface→beam mapping.\n"
                    f"Expected: '{self.model_part_name_beam}', but got: '{to_model_part_name}'."
                )

            if from_model_part_name != self.model_part_name_surface:
                raise Exception(
                    f"Invalid 'from_model_part_name' for surface→beam mapping.\n"
                    f"Expected: '{self.model_part_name_surface}', but got: '{from_model_part_name}'."
                )
        else:
            raise Exception(
                f"Invalid solver combination.\n"
                f"from_solver_name: '{from_solver_name}', to_solver_name: '{to_solver_name}'\n"
                f"Expected either (beam→surface) or (surface→beam)."
            )

        if self.beam_mapper == None:
            if inverse_map:
                origin_model_part = to_solver_data.GetModelPart() # beam
                destination_model_part = from_solver_data.GetModelPart() # surface
            else:
                origin_model_part = from_solver_data.GetModelPart() # beam
                destination_model_part = to_solver_data.GetModelPart() # surface
            
            self.beam_mapper = KratosMapping.MapperFactory.CreateMapper(origin_model_part, destination_model_part, self.settings["mapper_settings"].Clone())

        mapper_flags = self.__GetMapperFlags(transfer_options)

        if inverse_map:
            variable_destination = from_solver_data.variable
            variable_origin = to_solver_data.variable
            if not variable_origin == KSM.POINT_LOAD:
                raise Exception("Wrong variable")
            self.beam_mapper.InverseMap(variable_origin, self.second_variable_force, variable_destination, mapper_flags)
        else:
            variable_origin = from_solver_data.variable
            variable_destination = to_solver_data.variable
            if not variable_origin == KM.DISPLACEMENT:
                raise Exception("Wrong variable")
            self.beam_mapper.Map(variable_origin, self.second_variable_displacement, variable_destination, mapper_flags)

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "mapper_settings"              : { },
            "solver_name_beam"             : "UNSPECIFIED",
            "solver_name_surface"          : "UNSPECIFIED",
            "model_part_name_beam"         : "UNSPECIFIED",
            "model_part_name_surface"      : "UNSPECIFIED",
            "second_variable_displacement" : "ROTATION",
            "second_variable_force"        : "POINT_MOMENT",
            "type": "kratos_beam_mapping",
            "echo_level": 0
        }""")
        this_defaults.AddMissingParameters(super(KratosBeamMapping, cls)._GetDefaultParameterss())
        return this_defaults

    @classmethod
    def _GetListAvailableTransferOptions(cls):
        return cls.__mapper_flags_dict.keys()

    def __GetMapperFlags(self, transfer_options):
        mapper_flags = KM.Flags()
        for flag_name in transfer_options.GetStringArray():
            mapper_flags |= self.__mapper_flags_dict[flag_name]

        return mapper_flags