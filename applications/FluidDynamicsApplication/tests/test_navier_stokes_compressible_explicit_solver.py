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
        solver = navier_stokes_compressible_explicit_solver.CreateSolver(model, settings)

        # 4D is unsupported: an exception should be raised
        settings = self.GetSettings(4)
        with self.assertRaises(Exception) as context:
            solver = navier_stokes_compressible_explicit_solver.CreateSolver(model, settings)
        
        self.assertIn("Wrong domain size", str(context.exception))

    def test_CFL_number(self):
        solver_cfl_incompressible = self.CreateSolverAndModelPart(False)
        solver_cfl_compressible = self.CreateSolverAndModelPart(True)

        vel = 300
        dt  = 0.1
        c   = 0.0

        # Case 1: ensuring both results are the same with c=0.0
        v = KratosMultiphysics.Array3()
        v[0] = vel
        v[1] = 0.0
        v[2] = 0.0

        for node in solver_cfl_incompressible.GetComputingModelPart().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v)
            node.SetValue(KratosMultiphysics.SOUND_VELOCITY, c)
        
        for node in solver_cfl_compressible.GetComputingModelPart().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, v)
            node.SetValue(KratosMultiphysics.SOUND_VELOCITY, c)

        solver_cfl_compressible.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, dt)
        solver_cfl_incompressible.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, dt)

        dt_comp = solver_cfl_compressible._ComputeDeltaTime()
        dt_inc = solver_cfl_incompressible._ComputeDeltaTime()

        self.assertEqual(dt_comp, dt_inc, "Compressible and incompressible DELTA_TIME do not match with c=0")

        solver_cfl_compressible.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, dt)
        solver_cfl_incompressible.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, dt)

        # Case 2: Ensuring compressible dt is more restricted that incompressible with c > 0
        c = 340.0

        for node in solver_cfl_incompressible.GetComputingModelPart().Nodes:
            node.SetValue(KratosMultiphysics.SOUND_VELOCITY, c)
        
        for node in solver_cfl_compressible.GetComputingModelPart().Nodes:
            node.SetValue(KratosMultiphysics.SOUND_VELOCITY, c)

        dt_comp = solver_cfl_compressible._ComputeDeltaTime()
        dt_inc = solver_cfl_incompressible._ComputeDeltaTime()

        self.assertLess(dt_comp, dt_inc, "Compressible DELTA_TIME is not smaller than incompressible's despite having c=%f" % c)


    @classmethod
    def CreateSolverAndModelPart(cls, consider_compressibility):
        model = KratosMultiphysics.Model()
        mpart = model.CreateModelPart("FluidModelPart")
        mpart = mpart.CreateSubModelPart("volume_model_part")

        solver = navier_stokes_compressible_explicit_solver.CreateSolver(model, cls.GetSettings(2, consider_compressibility))
        solver.AddVariables()

        mpart.CreateNewNode(1, 0.0, 0.0, 0.0)
        mpart.CreateNewNode(2, 1e-3, 0.0, 0.0)
        mpart.CreateNewNode(3, 0.0, 1e-3, 0.0)
        mpart.CreateNewElement("Element2D3N", 1, [1,2,3], mpart.Properties[0])

        solver.ImportModelPart()
        solver.PrepareModelPart()
        solver.AddDofs()
        solver.Initialize()

        return solver


    @classmethod
    def GetSettings(cls, n_dimensions = 2, consider_compressibility = True):
        return KratosMultiphysics.Parameters("""
            {
                "solver_type": "compressible_solver_from_defaults",
                "model_part_name": "FluidModelPart",
                "domain_size": %d,
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                    "materials_filename": ""
                },
                "echo_level": 1,
                "time_order": 2,
                "move_mesh_flag": false,
                "shock_capturing": true,
                "compute_reactions": false,
                "reform_dofs_at_each_step" : false,
                "assign_neighbour_elements_to_conditions": true,
                "volume_model_part_name" : "volume_model_part",
                "skin_parts": [],
                "no_skin_parts":[],
                "time_stepping"                : {
                    "automatic_time_step" : true,
                    "CFL_number"          : 1.0,
                    "minimum_delta_time"  : 1.0e-8,
                    "maximum_delta_time"  : 1.0e-2,
                    "consider_compressibility" : %s
                },
                "use_oss" : true
            }""" % (n_dimensions, str(consider_compressibility).lower()))

if __name__ == '__main__':
    KratosUnittest.main()