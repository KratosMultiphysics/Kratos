# Import Kratos
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Additional imports
import constitutive_law_python_utility as constitutive_law_utils

def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName

class SprismTests(KratosUnittest.TestCase):

    def setUp(self):
        # Modelpart for the solid
        model_part = ModelPart("StructuralPart")
        model_part.AddNodalSolutionStepVariable(ALPHA_EAS)
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(REACTION)

        # Initialize GiD  I/O
        input_file_name = GetFilePath("SPRISM3D6N/patch_test")

        # Reading the fluid part
        model_part_io_fluid = ModelPartIO(input_file_name)
        model_part_io_fluid.ReadModelPart(model_part)
        
        # Find neighbours if required
        sprism_neighbour_search = SprismNeighbours(model_part)
        sprism_neighbour_search.Execute()

        # Setting up the buffer size
        model_part.SetBufferSize(3)

        # Adding dofs
        for node in model_part.Nodes:
          node.AddDof(DISPLACEMENT_X, REACTION_X)
          node.AddDof(DISPLACEMENT_Y, REACTION_Y)
          node.AddDof(DISPLACEMENT_Z, REACTION_Z)
          
        # Set the constitutive law
        constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(model_part, 3);
        constitutive_law.Initialize();

        max_iters = 30

        linear_solver = SkylineLUFactorizationSolver()

        # Create solver
        builder_and_solver = ResidualBasedBuilderAndSolver(linear_solver)
        mechanical_scheme = ResidualBasedStaticScheme()
        mechanical_convergence_criterion = ResidualCriteria(1e-4, 1e-4)

        self.mechanical_solver = ResidualBasedNewtonRaphsonStrategy(
            model_part,
            mechanical_scheme,
            linear_solver,
            mechanical_convergence_criterion,
            builder_and_solver,
            max_iters,
            False ,
            True,
            True)

        self.mechanical_solver.SetEchoLevel(0)
        self.mechanical_solver.Initialize()

        self.model_part = model_part

    def test_MembranePacth(self):
        step = 0
        time = 0.0
        Dt = 1.0
        final_time = 1.0

        # Add BC
        for node in self.model_part.Nodes:
            if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(
                    DISPLACEMENT_X,
                    0,
                    1.0e-7 * (node.X + node.Y / 2))
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(
                    DISPLACEMENT_Y,
                    0,
                    1.0e-7 * (node.Y + node.X / 2))

        while(time <= final_time):

            self.model_part.CloneTimeStep(time)

            time += Dt
            step += 1

            self.mechanical_solver.Solve()

        for node in self.model_part.Nodes:
            value = node.GetSolutionStepValue(DISPLACEMENT_X,0)
            self.assertAlmostEqual(value, 1.0e-7 * (node.X + node.Y / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            self.assertAlmostEqual(value, 1.0e-7 * (node.Y + node.X / 2))

    def test_BendingPacth(self):
        step = 0
        time = 0.0
        Dt = 1.0
        final_time = 1.0

        # Add BC
        for node in self.model_part.Nodes:
            if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(
                    DISPLACEMENT_X, 0,
                    -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(
                    DISPLACEMENT_Y, 0,
                    -1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
                node.Fix(DISPLACEMENT_Z)
                node.SetSolutionStepValue(
                    DISPLACEMENT_Z, 0,
                    0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y ** 2))

        while(time <= final_time):

            self.model_part.CloneTimeStep(time)

            time += Dt
            step += 1
            
            self.mechanical_solver.Solve()

        for node in self.model_part.Nodes:
            value = node.GetSolutionStepValue(DISPLACEMENT_X,0)
            self.assertAlmostEqual(
                value,
                -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            self.assertAlmostEqual(
                value,
                -1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Z,0)
            self.assertAlmostEqual(value,
                0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y **2 ))

    def tearDown(self):
        pass