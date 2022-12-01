from KratosMultiphysics import *
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.testing.utilities as TestingUtilities
import KratosMultiphysics.gid_output_process

class TestVariableRedistributionUtility(UnitTest.TestCase):

    def setUp(self):
        self.input_file = "redistribution_test"
        self.work_folder = "auxiliar_files_for_python_unittest/variable_redistribution_utility"

        self.domain_size = 3
        self.redistribution_iterations = 100
        self.redistribution_tolerance = 1e-8

        self.check_tolerance = 1e-6
        self.print_output = False

        self.mapped_function = None

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            KratosUtilities.DeleteTimeFiles('.')

    def test_linear_function(self):
        def Flag1Check(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) == 1.0

        def ReferenceSolution(node,variable):
            node.SetSolutionStepValue(variable,node.X)

        self.RunTestCase(
            Flag1Check,
            ReferenceSolution,
            PRESSURE,
            NODAL_PAUX,
            TEMPERATURE)

    def test_sharp_corners(self):
        def FlagDefinedCheck(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) > 0.0

        def ReferenceSolution(node,variable):
            node.SetSolutionStepValue(variable, 10.*node.X + node.Y)

        self.RunTestCase(
            FlagDefinedCheck,
            ReferenceSolution,
            PRESSURE,
            NODAL_PAUX,
            TEMPERATURE)

    def test_quadratic(self):
        def FlagDefinedCheck(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) > 0.0

        def ReferenceSolution(node,variable):
            node.SetSolutionStepValue(variable, node.X*node.X + node.Y*node.Z)

        self.RunTestCase(
            FlagDefinedCheck,
            ReferenceSolution,
            PRESSURE,
            NODAL_PAUX,
            TEMPERATURE)

    def test_vector(self):
        def FlagDefinedCheck(node):
            return node.GetSolutionStepValue(FLAG_VARIABLE) > 0.0

        def ReferenceSolution(node,variable):
            value = Array3()
            value[0] = node.Y
            value[1] = 10.*node.X + node.Z
            value[2] = 500
            node.SetSolutionStepValue(variable, value)

        self.RunTestCase(
            FlagDefinedCheck,
            ReferenceSolution,
            VELOCITY,
            VORTICITY,
            ACCELERATION)

    def test_quad_elements_quadratic_scalar_field(self):
        # Read mdpa and set test model part
        self.domain_size = 2
        self.input_file = "two_dim_unstructured_square_quads"
        self.work_folder = "auxiliar_files_for_python_unittest/mdpa_files"
        self._SetUpProblem()

        # Set scalar origin field
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(PRESSURE, 0, node.X**2 + node.Y**2)

        # Perform the variable redistribution test (first make point and then distribute again)
        reference_variable = PRESSURE
        intermediate_variable = NODAL_PAUX
        result_variable = TEMPERATURE

        VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self.model_part,
            self.model_part.Elements,
            reference_variable,
            intermediate_variable)

        VariableRedistributionUtility.DistributePointValues(
            self.model_part,
            self.model_part.Elements,
            intermediate_variable,
            result_variable,
            self.redistribution_tolerance,
            self.redistribution_iterations)

        if self.print_output:
            self._PrintOutput()

        # Check results
        self._CheckDoubleResults(self.model_part, reference_variable, result_variable)

    def test_quad_elements_quadratic_vector_field_non_historical(self):
        # Read mdpa and set test model part
        self.domain_size = 2
        self.input_file = "two_dim_unstructured_square_quads"
        self.work_folder = "auxiliar_files_for_python_unittest/mdpa_files"
        self._SetUpProblem()

        # Set scalar origin field
        for node in self.model_part.Nodes:
            node.SetValue(VELOCITY, [node.X**2 + node.Y**2, node.X + node.Y, 0.0])
        self.model_part.GetCommunicator().SynchronizeNonHistoricalVariable(VELOCITY)

        # Perform the variable redistribution test (first make point and then distribute again)
        reference_variable = VELOCITY
        intermediate_variable = VORTICITY
        result_variable = ACCELERATION

        VariableRedistributionUtility.ConvertDistributedValuesToPointNonHistorical(
            self.model_part,
            self.model_part.Elements,
            reference_variable,
            intermediate_variable)

        VariableRedistributionUtility.DistributePointValuesNonHistorical(
            self.model_part,
            self.model_part.Elements,
            intermediate_variable,
            result_variable,
            self.redistribution_tolerance,
            self.redistribution_iterations)

        if self.print_output:
            self._PrintOutput()

        # Check results
        self._CheckArrayResults(self.model_part, reference_variable, result_variable)

    def RunTestCase(self,inteface_check,set_reference,reference_variable,intermediate_variable,result_variable):
        # Read mdpa and set test model part
        self._SetUpProblem()

        for node in self.model_part.Nodes:
            if inteface_check(node):
                set_reference(node,reference_variable)
        self.model_part.GetCommunicator().SynchronizeVariable(reference_variable)

        VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self.model_part,
            self.model_part.Conditions,
            reference_variable,
            intermediate_variable)

        VariableRedistributionUtility.DistributePointValues(
            self.model_part,
            self.model_part.Conditions,
            intermediate_variable,
            result_variable,
            self.redistribution_tolerance,
            self.redistribution_iterations)

        if self.print_output:
            self._PrintOutput()

        if KratosGlobals.Kernel.HasDoubleVariable(reference_variable.Name()):
            self._CheckDoubleResults(self.model_part, reference_variable,result_variable)
        elif KratosGlobals.Kernel.HasArrayVariable(reference_variable.Name()):
            self._CheckArrayResults(self.model_part, reference_variable,result_variable)
        else:
            self.fail("Failing due to incorrect test definition: Wrong variable type")

    def test_nodal_area(self):
        self.input_file = "square10"

        self._SetUpProblem()

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(PRESSURE,1.0)
        self.model_part.GetCommunicator().SynchronizeVariable(PRESSURE)

        VariableRedistributionUtility.ConvertDistributedValuesToPoint(
            self.model_part,
            self.model_part.Conditions,
            PRESSURE,
            TEMPERATURE)

        for cond in self.model_part.Conditions:
            area = cond.GetGeometry().Area()
            for node in cond.GetNodes():
                nodal_area = node.GetSolutionStepValue(NODAL_PAUX)
                node.SetSolutionStepValue(NODAL_PAUX,nodal_area+area/3.0)
        self.model_part.GetCommunicator().AssembleCurrentData(NODAL_PAUX)

        if self.print_output:
            self._PrintOutput()

        self._CheckDoubleResults(self.model_part, TEMPERATURE, NODAL_PAUX)

    def _SetUpProblem(self):
        with UnitTest.WorkFolderScope(self.work_folder,__file__):
            self.current_model = Model()
            self.model_part = self.current_model.CreateModelPart("Interface")

            self.model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
            self.model_part.AddNodalSolutionStepVariable(PRESSURE)
            self.model_part.AddNodalSolutionStepVariable(NODAL_PAUX)
            self.model_part.AddNodalSolutionStepVariable(TEMPERATURE)
            self.model_part.AddNodalSolutionStepVariable(VELOCITY)
            self.model_part.AddNodalSolutionStepVariable(VORTICITY)
            self.model_part.AddNodalSolutionStepVariable(ACCELERATION)

            self.model_part.ProcessInfo[DOMAIN_SIZE] = self.domain_size
            self.model_part.SetBufferSize(1)

            TestingUtilities.ReadModelPart(self.input_file, self.model_part)

    def _CheckDoubleResults(self, model_part, reference_variable,result_variable):
        for node in model_part.Nodes:
            reference = node.GetSolutionStepValue(reference_variable)
            result = node.GetSolutionStepValue(result_variable)
            self.assertAlmostEqual(reference, result, delta=self.check_tolerance)

    def _CheckArrayResults(self, model_part, reference_variable,result_variable):
        for node in model_part.Nodes:
            reference = node.GetSolutionStepValue(reference_variable)
            result = node.GetSolutionStepValue(result_variable)
            self.assertAlmostEqual(reference[0], result[0], delta=self.check_tolerance)
            self.assertAlmostEqual(reference[1], result[1], delta=self.check_tolerance)
            self.assertAlmostEqual(reference[2], result[2], delta=self.check_tolerance)

    def _PrintOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters(r'''{
            "Parameters" : {
                "model_part_name" : "Interface",
                "output_name" : "test_variable_redistribution",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags": {
                            "GiDPostMode": "GiD_PostBinary",
                            "WriteDeformedMeshFlag": "WriteUndeformed",
                            "WriteConditionsFlag": "WriteConditions",
                            "MultiFileFlag": "SingleFile"
                        },
                        "nodal_results" : ["PRESSURE", "NODAL_PAUX", "TEMPERATURE", "VELOCITY", "VORTICITY", "ACCELERATION"]
                    }
                }
            }
        }''')
        gid_output_process = KratosMultiphysics.gid_output_process.Factory(gid_output_settings ,self.current_model)
        gid_output_process.ExecuteInitialize()
        gid_output_process.ExecuteBeforeSolutionLoop()
        gid_output_process.ExecuteInitializeSolutionStep()
        gid_output_process.PrintOutput()
        gid_output_process.ExecuteFinalizeSolutionStep()
        gid_output_process.ExecuteFinalize()

if __name__ == '__main__':
    UnitTest.main()
