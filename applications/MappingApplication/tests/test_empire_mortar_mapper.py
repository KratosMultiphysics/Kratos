import KratosMultiphysics as KM
from KratosMultiphysics import KratosUnittest
from KratosMultiphysics.MappingApplication import Mapper
from KratosMultiphysics.MappingApplication import empire_mortar_mapper
import mapper_test_case
from math import cos
import os


def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "mdpa_files", file_name)


class TestEmpireMortarMapper(KratosUnittest.TestCase):
    """This test depends on EMPIRE and is only used for debugging"""

    @classmethod
    def setUpClass(self):
        self.model = KM.Model()
        self.model_part_structure = self.model.CreateModelPart("structure")
        self.model_part_fluid = self.model.CreateModelPart("fluid")

        self.model_part_structure.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        self.model_part_structure.AddNodalSolutionStepVariable(KM.FORCE)

        self.model_part_fluid.AddNodalSolutionStepVariable(KM.MESH_DISPLACEMENT)
        self.model_part_fluid.AddNodalSolutionStepVariable(KM.REACTION)

        structure_mdpa_file_name = "blade_quad"
        fluid_mdpa_file_name     = "blade_tri"

        ReadModelPart(self.model_part_structure, structure_mdpa_file_name)
        ReadModelPart(self.model_part_fluid, fluid_mdpa_file_name)

    def setUp(self):
        mapper_settings = KM.Parameters("""{
            "path_mapper_lib" : "/home/philippb/software/empire/EMPIRE-Core/lib/libEMPIRE_MapperLib.so",
            "echo_level"      : 1
        }""")

        self.mapper_pressure_side = empire_mortar_mapper.Create(self.model_part_structure, self.model_part_fluid, mapper_settings)
        self.mapper_suction_side  = empire_mortar_mapper.Create(self.model_part_structure, self.model_part_fluid, mapper_settings)

        self.print_output = True

    def test_map_displacements(self):
        SetDisplacements(self.model_part_structure)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.DISPLACEMENT, "Blade_empire_mortar_Structure_prescr_disp")

        self.mapper_pressure_side.Map(KM.DISPLACEMENT, KM.MESH_DISPLACEMENT)
        self.mapper_suction_side.Map(KM.DISPLACEMENT, KM.MESH_DISPLACEMENT)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.MESH_DISPLACEMENT, "Blade_empire_mortar_Fluid_mapped_disp")

    def test_map_forces(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_empire_mortar_Fluid_prescr_force")

        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.REACTION, Mapper.SWAP_SIGN)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.REACTION, Mapper.SWAP_SIGN)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_empire_mortar_Structure_mapped_force")

    def test_map_forces_conservative(self):
        SetReactions(self.model_part_fluid)
        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_fluid, KM.REACTION, "Blade_empire_mortar_Fluid_prescr_force")

        self.mapper_pressure_side.InverseMap(KM.FORCE, KM.REACTION, Mapper.SWAP_SIGN | Mapper.USE_TRANSPOSE)
        self.mapper_suction_side.InverseMap(KM.FORCE, KM.REACTION, Mapper.SWAP_SIGN | Mapper.USE_TRANSPOSE)

        if self.print_output:
            mapper_test_case.VtkOutputNodesHistorical(self.model_part_structure, KM.FORCE, "Blade_empire_mortar_Structure_mapped_force_conserv")


def ReadModelPart(model_part, mdpa_file_name):
    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
    KM.ModelPartIO(GetFilePath(mdpa_file_name), import_flags).ReadModelPart(model_part)

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


if __name__ == '__main__':
    KratosUnittest.main()
