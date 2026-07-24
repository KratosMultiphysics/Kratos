import KratosMultiphysics
from KratosMultiphysics import IsDistributedRun

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtils


def SelectAndVerifyLinearSolver(settings, skiptest):
    # The mechanical solver selects automatically the fastest linear-solver available
    # this might not be appropriate for a test, therefore in case nothing is specified,
    # the previous default linear-solver is set
    if not settings["solver_settings"].Has("linear_solver_settings"):
        # check if running in MPI because there we use a different default linear solver
        if IsDistributedRun():
            default_lin_solver_settings = KratosMultiphysics.Parameters(
                """{
                "solver_type" : "amesos",
                "amesos_solver_type" : "Amesos_Klu"
            }"""
            )

        else:
            default_lin_solver_settings = KratosMultiphysics.Parameters(
                """{
                "solver_type": "LinearSolversApplication.sparse_lu"
            }"""
            )
        settings["solver_settings"].AddValue(
            "linear_solver_settings", default_lin_solver_settings
        )

    solver_type = settings["solver_settings"]["linear_solver_settings"][
        "solver_type"
    ].GetString()
    solver_type_splitted = solver_type.split(".")
    if len(solver_type_splitted) == 2:
        # this means that we use a solver from an application
        # hence we have to check if it exists, otherwise skip the test
        app_name = solver_type_splitted[0]
        solver_name = solver_type_splitted[1]
        if not kratos_utils.CheckIfApplicationsAvailable(app_name):
            skiptest(
                'Application "{}" is needed for the specified solver "{}" but is not available'.format(
                    app_name, solver_name
                )
            )


class TestDataGenerationFactory(KratosUnittest.TestCase):
    def _check_data_dimensions(self):
        input_data = DataLoadingUtils.ImportDataFromFile(
            self.output_file_in, "InputData"
        )
        self.assertEqual(input_data[0].data.size, self.dimension_in)
        # os.remove(self.output_file_in)
        output_data = DataLoadingUtils.ImportDataFromFile(
            self.output_file_out, "OutputData"
        )
        self.assertEqual(output_data[0].data.size, self.dimension_out)
        # os.remove(self.output_file_out)

    def setUp(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            # Reading the ProjectParameters
            with open(self.file_name + "_parameters.json", "r") as parameter_file:
                ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            SelectAndVerifyLinearSolver(ProjectParameters, self.skipTest)

            self.modify_parameters(ProjectParameters)

            # To avoid many prints
            if ProjectParameters["problem_data"]["echo_level"].GetInt() == 0:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
                    KratosMultiphysics.Logger.Severity.WARNING
                )
            else:
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(
                    KratosMultiphysics.Logger.Severity.INFO
                )

            # Creating the test
            model = KratosMultiphysics.Model()
            self.test = StructuralMechanicsAnalysis(model, ProjectParameters)
            self.test.Initialize()

    def modify_parameters(self, project_parameters):
        """This function can be used in derived classes to modify existing parameters
        before the execution of the test (e.g. switch to MPI)
        """
        pass

    def test_execution(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.RunSolutionLoop()
            self._check_data_dimensions()

    def tearDown(self):
        # Within this location context:
        with KratosUnittest.WorkFolderScope(".", __file__):
            self.test.Finalize()
            # self._check_data_dimensions()


class TestDataGeneration1I1OVariablesText(TestDataGenerationFactory):
    file_name = "beam_cantilever_nonlinear/1I1OVariablesText"
    output_file_in = "data/training_in_raw.dat"
    output_file_out = "data/training_out_raw.dat"
    dimension_in = 1
    dimension_out = 1


class TestDataGeneration1I2OVariablesText(TestDataGenerationFactory):
    file_name = "beam_cantilever_nonlinear/1I2OVariablesText"
    output_file_in = "data/training_in_raw.dat"
    output_file_out = "data/training_out_raw.dat"
    dimension_in = 1
    dimension_out = 2


class TestDataGeneration1IVectorOVariablesText(TestDataGenerationFactory):
    file_name = "beam_cantilever_nonlinear/1IVectorOVariablesText"
    output_file_in = "data/training_in_raw.dat"
    output_file_out = "data/training_out_raw.dat"
    dimension_in = 1
    dimension_out = 3


class TestDataGeneration1I1OVariablesHDF5(TestDataGenerationFactory):
    file_name = "beam_cantilever_nonlinear/1I1OVariablesHDF5"
    output_file_in = "data/training_in_raw_1o.h5"
    output_file_out = "data/training_out_raw_1o.h5"
    dimension_in = 1
    dimension_out = 1


class TestDataGeneration1I2OVariablesHDF5(TestDataGenerationFactory):
    file_name = "beam_cantilever_nonlinear/1I2OVariablesHDF5"
    output_file_in = "data/training_in_raw_2o.h5"
    output_file_out = "data/training_out_raw_2o.h5"
    dimension_in = 1
    dimension_out = 2


class TestDataGeneration1IVectorOVariablesHDF5(TestDataGenerationFactory):
    file_name = "beam_cantilever_nonlinear/1IVectorOVariablesHDF5"
    output_file_in = "data/training_in_raw_vector.h5"
    output_file_out = "data/training_out_raw_vector.h5"
    dimension_in = 1
    dimension_out = 3


if __name__ == "__main__":
    KratosUnittest.main()
