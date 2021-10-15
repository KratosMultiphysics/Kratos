import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
import KratosMultiphysics.MappingApplication as KratosMapping
default_data_comm = KM.Testing.GetDefaultDataCommunicator()
if default_data_comm.IsDistributed():
    from KratosMultiphysics import mpi as KratosMPI
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

from KratosMultiphysics.testing import utilities as testing_utils
import mapper_test_case
from math import cos
import os
from sys import version_info as py_version_info

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class BladeMappingTests(mapper_test_case.MapperTestCase):
    '''This class contains basic tests for mapping on real geometries
    In this case it is a remodeled NREL Phase VI wind turbine blade
    It also serves as a showcase on how to use the Mapper in FSI
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        structure_mdpa_file_name = "blade_quad"
        fluid_mdpa_file_name     = "blade_tri"
        super().setUpModelParts(structure_mdpa_file_name, fluid_mdpa_file_name)

        # TODO ATTENTION: currently the MapperFactory removes some keys, hence those checks have to be done beforehand => improve this!

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        cls.model_part_structure = cls.model_part_origin
        cls.model_part_fluid = cls.model_part_destination

        if default_data_comm.IsDistributed():
            map_creator = MappingMPIExtension.MPIMapperFactory.CreateMapper
        else:
            map_creator = KratosMapping.MapperFactory.CreateMapper

        cls.mapper_pressure_side = map_creator(
            cls.GetStructureModelPart("pressure_side_quad"),
            cls.model_part_fluid.GetSubModelPart("pressure_side_tri"),
            mapper_parameters.Clone())

        cls.mapper_suction_side = map_creator(
            cls.GetStructureModelPart("suction_side_quad"),
            cls.model_part_fluid.GetSubModelPart("suction_side_tri"),
            mapper_parameters.Clone())

        cls.print_output = False # this can be overridden in derived classes to print the output

    @classmethod
    def GetStructureModelPart(cls, sub_model_part_name):
        return cls.model_part_structure.GetSubModelPart(sub_model_part_name)

    def test_map_displacements(self):
        SetDisplacements(self.model_part_structure)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.DISPLACEMENT, "Blade_" + self.mapper_type + "_Structure_prescr_disp")

        self.mapper_pressure_side.Map(KM.DISPLACEMENT, KM.MESH_DISPLACEMENT)
        self.mapper_suction_side.Map(KM.DISPLACEMENT, KM.MESH_DISPLACEMENT)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.MESH_DISPLACEMENT, "Blade_" + self.mapper_type + "_Fluid_mapped_disp")

        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_fluid, KM.MESH_DISPLACEMENT, GetFilePath(self._GetFileName("balde_map_disp")))

    def test_map_forces(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_" + self.mapper_type + "_Fluid_prescr_force")

        # this would be POINT_LOAD in regular StructuralMechanics (using FORCE to avoid the StructuralMechanics import)
        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.REACTION, KratosMapping.Mapper.SWAP_SIGN)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.REACTION, KratosMapping.Mapper.SWAP_SIGN)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_" + self.mapper_type + "_Structure_mapped_force")

        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.FORCE, GetFilePath(self._GetFileName("balde_map_force")))

    def test_map_forces_conservative(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_" + self.mapper_type + "_Fluid_prescr_force")

         # this would be POINT_LOAD in regular StructuralMechanics (using FORCE to avoid the StructuralMechanics import)
        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.REACTION, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.USE_TRANSPOSE)
        self.__CheckValuesSum(self.model_part_fluid.GetSubModelPart("pressure_side_tri"), self.GetStructureModelPart("pressure_side_quad"), KM.REACTION, KM.FORCE)

        # Note: Setting the solution again because in this case some nodes are shared and hence
        # would slightly influence the computation
        SetReactions(self.model_part_fluid)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.REACTION, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.USE_TRANSPOSE)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_" + self.mapper_type + "_Structure_mapped_force_conserv")

        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.FORCE, GetFilePath(self._GetFileName("balde_map_force_conserv")))

        self.__CheckValuesSum(self.model_part_fluid.GetSubModelPart("suction_side_tri"), self.GetStructureModelPart("suction_side_quad"), KM.REACTION, KM.FORCE)

    def __CheckValuesSum(self, mp1, mp2, var1, var2):
        if mp1.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank():
            val_1 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var1, mp1, 0)

        if mp2.GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank():
            val_2 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var2, mp2, 0)

        if mp1.GetCommunicator().GetDataCommunicator().IsNullOnThisRank() or mp2.GetCommunicator().GetDataCommunicator().IsNullOnThisRank():
            return # ranks that don't have either ModelPart don't do the check

        # minus because SWAP_SIGN is used
        self.assertVectorAlmostEqual(val_1, (-1)*val_2)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, file_appendix)

def SetDisplacements(model_part_structure):
    for node in model_part_structure.Nodes:
        disp_x = 0.0 # edgewise
        disp_y = 0.008*(node.Z)**3 # flapwise
        disp_z = 0.0 # spanwise
        node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([disp_x, disp_y, disp_z]))

def SetReactions(model_part_fluid):
    for node in model_part_fluid.Nodes:
        react_x = cos(node.X*5)
        react_y = 0.0
        react_z = cos(node.Z*2)
        node.SetSolutionStepValue(KM.REACTION, KM.Vector([react_x, react_y, react_z]))


@KratosUnittest.skipUnless(KM.IsDistributedRun(), "this test requires MPI")
@KratosUnittest.skipIf((py_version_info[0], py_version_info[1]) < (3,8), 'this test requires at least python 3.8 (bcs of "addClassCleanup"')
class BladeMappingTestsSerialModelPart(BladeMappingTests):
    '''In these tests the structural ModelPart is serial and the fluid ModelPart is distributed
    '''

    @classmethod
    def ReadModelParts(cls):
        data_comm_name_rank_0 = "OnlyRank0"
        cls.sub_comm_rank_0 = KratosMPI.DataCommunicatorFactory.CreateFromRanksAndRegister(default_data_comm, [0], data_comm_name_rank_0)
        cls.addClassCleanup(KM.ParallelEnvironment.UnregisterDataCommunicator, data_comm_name_rank_0)

        if default_data_comm.Rank() == 0:
            testing_utils.ReadSerialModelPart(cls.input_file_origin, cls.model_part_origin) # structure
        else:
            # set the DataComm that only has rank 0 in the "dummy" ModelPart to be used in the mapper
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(cls.model_part_origin, cls.sub_comm_rank_0)

        testing_utils.ReadDistributedModelPart(cls.input_file_destination, cls.model_part_destination) # fluid

    @classmethod
    def GetStructureModelPart(cls, sub_model_part_name):
        if default_data_comm.Rank() == 0:
            return cls.model_part_structure.GetSubModelPart(sub_model_part_name)
        else:
            # other ranks return a dummy ModelPart (here using the MainModelPart of the structure as example)
            return cls.model_part_structure

@KratosUnittest.skipUnless(default_data_comm.Size() > 2, "This test needs at least 3 mpi processes")
@KratosUnittest.skipIf((py_version_info[0], py_version_info[1]) < (3,8), 'this test requires at least python 3.8 (bcs of "addClassCleanup"')
class BladeMappingTestsLessRanksModelPart(BladeMappingTests):
    '''In these tests the structural ModelPart is distributed over less ranks
    and the fluid ModelPart is distributed over all ranks
    '''

    @classmethod
    def GetRanksForStructure(cls):
        raise NotImplementedError("This function mus be implemented in the derived classes!")

    @classmethod
    def ReadModelParts(cls):
        ranks_structure = cls.GetRanksForStructure()
        data_comm_name = "structure_comm"
        cls.sub_comm_structure = KratosMPI.DataCommunicatorFactory.CreateFromRanksAndRegister(default_data_comm, ranks_structure, data_comm_name)
        cls.addClassCleanup(KM.ParallelEnvironment.UnregisterDataCommunicator, data_comm_name)

        importer_settings = KM.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + cls.input_file_origin + """\",
                "partition_in_memory" : true,
                "data_communicator_name": \"""" + data_comm_name + """\"
            },
            "echo_level" : 0
        }""")
        if cls.sub_comm_structure.IsDefinedOnThisRank():
            testing_utils.ReadDistributedModelPart(cls.input_file_origin, cls.model_part_origin, importer_settings) # structure
        else:
            # set the DataComm of the structure also in the other ranks, where it is not defined!
            KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(cls.model_part_origin, cls.sub_comm_structure)

        testing_utils.ReadDistributedModelPart(cls.input_file_destination, cls.model_part_destination) # fluid

    @classmethod
    def GetStructureModelPart(cls, sub_model_part_name):
        if cls.sub_comm_structure.IsDefinedOnThisRank():
            return cls.model_part_structure.GetSubModelPart(sub_model_part_name)
        else:
            # other ranks return a dummy ModelPart (here using the MainModelPart of the structure as example)
            return cls.model_part_structure

class BladeMappingTestsAllRanksExceptLast(BladeMappingTestsLessRanksModelPart):
    @classmethod
    def GetRanksForStructure(cls):
        num_ranks = default_data_comm.Size()
        return list(range(num_ranks-1))

class BladeMappingTestsAllRanksExceptFirst(BladeMappingTestsLessRanksModelPart):
    @classmethod
    def GetRanksForStructure(cls):
        num_ranks = default_data_comm.Size()
        return list(range(1, num_ranks))

class BladeMappingTestsUnevenRanks(BladeMappingTestsLessRanksModelPart):
    @classmethod
    def GetRanksForStructure(cls):
        num_ranks = default_data_comm.Size()
        return list(range(1, num_ranks, 2))

