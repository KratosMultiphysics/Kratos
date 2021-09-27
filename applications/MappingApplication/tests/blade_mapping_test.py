import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
data_comm = KM.Testing.GetDefaultDataCommunicator()
if data_comm.IsDistributed():
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

import mapper_test_case
from math import cos
import os

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
        super(BladeMappingTests, cls).setUpModelParts(structure_mdpa_file_name, fluid_mdpa_file_name)

        # TODO ATTENTION: currently the MapperFactory removes some keys, hence those checks have to be done beforehand => improve this!

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        mapper_parameters_p = mapper_parameters.Clone()
        mapper_parameters_s = mapper_parameters.Clone()

        mapper_parameters_p.AddEmptyValue("interface_submodel_part_origin").SetString("pressure_side_quad")
        mapper_parameters_p.AddEmptyValue("interface_submodel_part_destination").SetString("pressure_side_tri")

        mapper_parameters_s.AddEmptyValue("interface_submodel_part_origin").SetString("suction_side_quad")
        mapper_parameters_s.AddEmptyValue("interface_submodel_part_destination").SetString("suction_side_tri")

        cls.model_part_structure = cls.model_part_origin
        cls.model_part_fluid = cls.model_part_destination

        if data_comm.IsDistributed():
            cls.mapper_pressure_side = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
            cls.mapper_suction_side = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)
        else:
            cls.mapper_pressure_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_p)
            cls.mapper_suction_side = KratosMapping.MapperFactory.CreateMapper(cls.model_part_structure, cls.model_part_fluid, mapper_parameters_s)

        cls.print_output = False # this can be overridden in derived classes to print the output

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
        self.__CheckValuesSum(self.model_part_fluid.GetSubModelPart("pressure_side_tri"), self.model_part_structure.GetSubModelPart("pressure_side_quad"), KM.REACTION, KM.FORCE)

        # Note: Setting the solution again because in this case some nodes are shared and hence
        # would slightly influence the computation
        SetReactions(self.model_part_fluid)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.REACTION, KratosMapping.Mapper.SWAP_SIGN | KratosMapping.Mapper.USE_TRANSPOSE)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_" + self.mapper_type + "_Structure_mapped_force_conserv")

        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_structure, KM.FORCE, GetFilePath(self._GetFileName("balde_map_force_conserv")))

        self.__CheckValuesSum(self.model_part_fluid.GetSubModelPart("suction_side_tri"), self.model_part_structure.GetSubModelPart("suction_side_quad"), KM.REACTION, KM.FORCE)

    def __CheckValuesSum(self, mp1, mp2, var1, var2):
        val_1 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var1, mp1, 0)
        val_2 = KM.VariableUtils().SumHistoricalNodeVectorVariable(var2, mp2, 0)
        # minus because SWAP_SIGN is used
        self.assertAlmostEqual(val_1[0], -val_2[0])
        self.assertAlmostEqual(val_1[1], -val_2[1])
        self.assertAlmostEqual(val_1[2], -val_2[2])

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
