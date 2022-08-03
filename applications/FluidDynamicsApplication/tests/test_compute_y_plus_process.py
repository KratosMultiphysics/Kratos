import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.process_factory import KratosProcessFactory

class ComputeYPlusProcessTest(UnitTest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 2
        cls.model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        ReadModelPart("AdjointVMSSensitivity2DTest/cylinder_test", cls.model_part)

        cls.density = 2.0
        cls.kinematic_viscosity = 3.0

        for element in cls.model_part.Elements:
            element.Properties[Kratos.DENSITY] = cls.density
            element.Properties[Kratos.DYNAMIC_VISCOSITY] = cls.kinematic_viscosity * cls.density
            break

        tmoc = Kratos.TetrahedralMeshOrientationCheck
        flags = (tmoc.COMPUTE_NODAL_NORMALS).AsFalse() | (tmoc.COMPUTE_CONDITION_NORMALS).AsFalse() | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS
        Kratos.TetrahedralMeshOrientationCheck(cls.model_part, False, flags).Execute()

        # now populate the reactions
        for node in cls.model_part.Nodes:
            node.SetSolutionStepValue(Kratos.REACTION, Kratos.Array3([node.Id, node.Id * 2, node.Id * 3]))

    def testComputeYPlusProcess(self):
        parameters = Kratos.Parameters("""{
            "model_part_name"                  : "test",
            "output_variable_name"             : "Y_PLUS",
            "output_to_elements"               : true,
            "calculate_normals_every_time_step": false,
            "echo_level"                       : 1
        }""")
        compute_y_plus_process = KratosCFD.ComputeYPlusProcess(self.model, parameters)
        compute_y_plus_process.Check()
        compute_y_plus_process.ExecuteInitializeSolutionStep()
        compute_y_plus_process.ExecuteFinalizeSolutionStep()

        ref_y_values = {
            1:0.000255058,
            2:0.000277998,
            3:0.000262711,
            4:0.000262711,
            5:0.000277998,
            6:0.000255058,
            7:0.000300881,
            8:0.000300881,
            9:0.000331303,
            10:0.000331303,
        }

        # now check the shear stresses
        for condition_id, y in ref_y_values.items():
            condition = self.model_part.GetCondition(condition_id)
            geometry = condition.GetGeometry()
            y_plus = condition.GetValue(KratosCFD.Y_PLUS)
            u_tau = y_plus * self.kinematic_viscosity / y
            shear_stress = self.density * (u_tau ** 2)
            y_plus_based_tangential_reaction_magnitude = shear_stress * geometry.Area()

            reaction = Kratos.Array3([0.0, 0.0, 0.0])
            for node in geometry:
                nodal_normal = node.GetSolutionStepValue(Kratos.NORMAL)
                nodal_normal_magnitude = math.sqrt(ComputeYPlusProcessTest.__InnerProd(nodal_normal, nodal_normal))
                reaction += node.GetSolutionStepValue(Kratos.REACTION) / nodal_normal_magnitude
            reaction *= geometry.Area() / len(geometry)

            condition_unit_normal = condition.GetValue(Kratos.NORMAL)
            condition_unit_normal /= math.sqrt(ComputeYPlusProcessTest.__InnerProd(condition_unit_normal, condition_unit_normal))
            perpendicular_reaction = condition_unit_normal * ComputeYPlusProcessTest.__InnerProd(reaction, condition_unit_normal)
            analytical_tangential_reaction = reaction - perpendicular_reaction
            analytical_tangential_reaction_magnitude = math.sqrt(ComputeYPlusProcessTest.__InnerProd(analytical_tangential_reaction, analytical_tangential_reaction))

            self.assertAlmostEqual(y_plus_based_tangential_reaction_magnitude/analytical_tangential_reaction_magnitude - 1.0, 0.0, 5)

    def testYPlusOutputProcess(self):
        settings = Kratos.Parameters(r'''
        [
            {
                "kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
                "python_module" : "y_plus_output_process",
                "Parameters" : {
                    "model_part_name"                  : "test",
                    "interval"                         : [0.0, 1e30],
                    "y_plus_output_limit"              : 0.1,
                    "print_to_screen"                  : false,
                    "print_format"                     : ".8f",
                    "write_output_file"                : true,
                    "output_step"                      : 2,
                    "output_to_elements"               : false,
                    "calculate_normals_every_time_step": false,
                    "echo_level"                       : 0,
                    "output_file_settings"             : {}
                }
            }
        ]''')

        factory = KratosProcessFactory(self.model)
        process_list = factory.ConstructListOfProcesses(settings)

        for process in process_list:
            process.Check()
            
        time_steps = [1.0, 1.5, 2.0, 2.5]
        for step, time_step in enumerate(time_steps, 1):
            self.model_part.ProcessInfo[Kratos.TIME] = time_step
            self.model_part.ProcessInfo[Kratos.STEP] = step

            for process in process_list:
                process.ExecuteInitializeSolutionStep()
                process.ExecuteFinalizeSolutionStep()
                if process.IsOutputStep():
                    process.PrintOutput()

        for process in process_list:
            process.ExecuteFinalize()

    def tearDown(self):
        kratos_utilities.DeleteFileIfExisting("test_y_plus.dat")

    @staticmethod
    def __InnerProd(v1: Kratos.Array3, v2: Kratos.Array3):
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

if __name__ == '__main__':
    UnitTest.main()
