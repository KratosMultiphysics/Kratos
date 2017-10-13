from KratosMultiphysics import *
from KratosMultiphysics.FDApplication import *

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class Bfecc(KratosUnittest.TestCase):

    def setUp(self):
        # Modelpart for the fluid
        model_part = ModelPart("FluidPart")

        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part.AddNodalSolutionStepVariable(DENSITY)

        # initialize GiD  I/O
        input_file_name = GetFilePath("IsoCubeFluid")

        # reading the fluid part
        model_part_io_fluid = ModelPartIO(input_file_name)
        model_part_io_fluid.ReadModelPart(model_part)

        # setting up the buffer size: SHOULD BE DONE AFTER READING!!!
        model_part.SetBufferSize(1)

        self.model_part = model_part

    def tearDown(self):
        pass

    def testCreateSolver(self):
        solver = BfeccSolverStrategy(self.model_part)
        solver.Initialize()

    def testReadIncorrectModelpart(self):
        model_part = ModelPart("FluidPart")

        # Density is missing
        model_part.AddNodalSolutionStepVariable(PRESSURE)
        model_part.AddNodalSolutionStepVariable(VELOCITY)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)

        # initialize GiD  I/O
        input_file_name = GetFilePath("IsoCubeFluid")

        # reading the fluid part
        with self.assertRaisesRegex(RuntimeError, "kratos"):
            model_part_io_fluid = ModelPartIO(input_file_name)
            model_part_io_fluid.ReadModelPart(model_part)
