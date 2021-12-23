import KratosMultiphysics
from  KratosMultiphysics import FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication import navier_stokes_compressible_explicit_solver

from KratosMultiphysics import KratosUnittest
from KratosMultiphysics import kratos_utilities

class NavierStokesCompressibleExplicitSolverTest(KratosUnittest.TestCase):
    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = True

    def test_DimensionNumber(self):
        model = KratosMultiphysics.Model()

        # 3D is supported: no exception should be raised
        settings = self.GetSettings(3)
        _ = navier_stokes_compressible_explicit_solver.CreateSolver(model, settings)

        # 4D is unsupported: an exception should be raised
        settings = self.GetSettings(4)
        with self.assertRaises(Exception) as context:
            _ = navier_stokes_compressible_explicit_solver.CreateSolver(model, settings)

        self.assertIn("Wrong domain size", str(context.exception))

    def test_ChoiceRungeKutta(self):
        model = KratosMultiphysics.Model()
        model.CreateModelPart("FluidModelPart").CreateSubModelPart("fluid_computational_model_part")

        settings = self.GetSettings(3)
        settings["time_scheme"].SetString("WrongValue")

        # Ensuring helpful error
        solver = navier_stokes_compressible_explicit_solver.CreateSolver(model, settings)
        solver.ImportModelPart()
        with self.assertRaises(Exception) as context:
            solver.Initialize()
        self.assertIn("WrongValue", str(context.exception))
        self.assertIn("RK4", str(context.exception))
        self.assertIn("RK3-TVD", str(context.exception))
        self.assertIn("forward_euler", str(context.exception))

    @classmethod
    def GetSettings(cls, n_dimensions):
        return KratosMultiphysics.Parameters("""
            {
                "solver_type": "compressible_solver_from_defaults",
                "model_part_name": "FluidModelPart",
                "domain_size": %d,
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                    "materials_filename": "FluidMaterials.json"
                },
                "echo_level": 1,
                "time_order": 2,
                "move_mesh_flag": false,
                "time_scheme" : "forward_euler",
                "compute_reactions": false,
                "reform_dofs_at_each_step" : false,
                "assign_neighbour_elements_to_conditions": true,
                "volume_model_part_name" : "volume_model_part",
                "skin_parts": [""],
                "no_skin_parts":[""],
                "time_stepping"                : {
                    "automatic_time_step" : true,
                    "CFL_number"          : 1.0,
                    "minimum_delta_time"  : 1.0e-8,
                    "maximum_delta_time"  : 1.0e-2
                },
                "use_oss" : true
            }""" % n_dimensions)

if __name__ == '__main__':
    KratosUnittest.main()
