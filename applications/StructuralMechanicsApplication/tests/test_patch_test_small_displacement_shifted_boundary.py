import KratosMultiphysics

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestPatchTestSmallDisplacementShiftedBoundary(KratosUnittest.TestCase):

    def setUp(self):
        self.print_output = False

    def testPatchTestSmallDisplacementShiftedBoundary2D3N(self):
        self.work_folder = "patch_test/small_disp_shifted_boundary"
        self.file_name = "patch_test_2D_tension_tri_shifted_boundary"
        self.__RunTest()

    def testSmallDisplacementShiftedBoundaryLocalSystemConsistency2D3N(self):
        self.work_folder = "patch_test/small_disp_shifted_boundary"
        self.file_name = "patch_test_2D_tension_tri_shifted_boundary"

        class SmallDisplacementShiftedBoundaryPatchTestAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):

            def ModifyInitialGeometry(self):
                super().ModifyInitialGeometry()
                level_set_x_position = 0
                for node in self._GetSolver().GetComputingModelPart().Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, node.X - level_set_x_position)

            def ApplyBoundaryConditions(self):
                super().ApplyBoundaryConditions()

                penalty_factor = 1.0e3
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.PENALTY_COEFFICIENT] = penalty_factor
                for condition in self._GetSolver().GetComputingModelPart().GetSubModelPart("shifted_boundary").Conditions:
                    condition.SetValue(KratosMultiphysics.DISPLACEMENT, [0.0,0.0,0.0])

        with open(f"{self.work_folder}/{self.file_name}_parameters.json",'r') as parameter_file:
            self.settings = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = SmallDisplacementShiftedBoundaryPatchTestAnalysis(model, self.settings)
            simulation.Run()

            computing_model_part = simulation._GetSolver().GetComputingModelPart()
            process_info = computing_model_part.ProcessInfo

            # Check every DisplacementShiftedBoundaryCondition on the immersed boundary
            n_checked_conditions = 0
            for condition in computing_model_part.GetSubModelPart("shifted_boundary").Conditions:
                self.__AssertLocalSystemConsistency(condition, process_info)
                n_checked_conditions += 1
            self.assertGreater(n_checked_conditions, 0)

            # Check every SmallDisplacementShiftedBoundaryElement touching the surrogate boundary
            n_checked_elements = 0
            for element in computing_model_part.Elements:
                if element.Is(KratosMultiphysics.INTERFACE):
                    self.__AssertLocalSystemConsistency(element, process_info)
                    n_checked_elements += 1
            self.assertGreater(n_checked_elements, 0)

    def __AssertLocalSystemConsistency(self, entity, process_info):
        lhs_combined = KratosMultiphysics.Matrix(0, 0)
        rhs_combined = KratosMultiphysics.Vector(0)
        entity.CalculateLocalSystem(lhs_combined, rhs_combined, process_info)

        lhs_split = KratosMultiphysics.Matrix(0, 0)
        entity.CalculateLeftHandSide(lhs_split, process_info)

        rhs_split = KratosMultiphysics.Vector(0)
        entity.CalculateRightHandSide(rhs_split, process_info)

        self.assertEqual(lhs_combined.Size1(), lhs_split.Size1())
        self.assertEqual(lhs_combined.Size2(), lhs_split.Size2())
        for i in range(lhs_combined.Size1()):
            for j in range(lhs_combined.Size2()):
                a = lhs_combined[i, j]
                b = lhs_split[i, j]
                self.assertAlmostEqual(a, b, delta=1e-9 * (abs(a) + abs(b) + 1.0))

        self.assertEqual(rhs_combined.Size(), rhs_split.Size())
        for i in range(rhs_combined.Size()):
            a = rhs_combined[i]
            b = rhs_split[i]
            self.assertAlmostEqual(a, b, delta=1e-9 * (abs(a) + abs(b) + 1.0))

    def __RunTest(self):
        class SmallDisplacementShiftedBoundaryPatchTestAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):

            def ModifyInitialGeometry(self):
                super().ModifyInitialGeometry()
                level_set_x_position = 0
                for node in self._GetSolver().GetComputingModelPart().Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, node.X - level_set_x_position)

            def ApplyBoundaryConditions(self):
                super().ApplyBoundaryConditions()

                penalty_factor = 1.0e3
                self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.PENALTY_COEFFICIENT] = penalty_factor
                for condition in self._GetSolver().GetComputingModelPart().GetSubModelPart("shifted_boundary").Conditions:
                    condition.SetValue(KratosMultiphysics.DISPLACEMENT, [0.0,0.0,0.0])

        with open(f"{self.work_folder}/{self.file_name}_parameters.json",'r') as parameter_file:
            # Read and customize settings
            self.settings = KratosMultiphysics.Parameters(parameter_file.read())
            if self.print_output:
                self.__AddOutput()

            # Creating the test
            model = KratosMultiphysics.Model()
            simulation = SmallDisplacementShiftedBoundaryPatchTestAnalysis(model, self.settings)
            simulation.Run()

    def __AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "Parameters"    : {
                "model_part_name"        : "Structure",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_interval"             : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["DISPLACEMENT","DISTANCE"],
                        "nodal_flags_results"         : ["INTERFACE","BOUNDARY","ACTIVE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        gid_output_settings["Parameters"]["output_name"].SetString(f"{self.work_folder}/{self.file_name}")
        self.settings["output_processes"]["gid_output"].Append(gid_output_settings)

class TestPatchTestShellShiftedBoundary(KratosUnittest.TestCase):
    E = 200000.0
    NU = 0.0
    L = 4.0
    Q_LOAD = 1.0
    THICKNESS = 1.0

    @classmethod
    def _ExactMidsurfaceDeflection(cls, x, t):
        # Guided (x=0) - simply-supported (x=L) beam strip under a uniformly distributed
        # transverse load q: w(x,t) = q/(2*kappa*G*t)*(L^2-x^2) + q/(2*E*t^3)*(5*L^4-6*L^2*x^2+x^4)
        G = cls.E / (2.0 * (1.0 + cls.NU))
        kappa = 5.0 / 6.0
        q = cls.Q_LOAD
        return (q / (2.0 * kappa * G * t)) * (cls.L**2 - x**2) \
            + (q / (2.0 * cls.E * t**3)) * (5.0 * cls.L**4 - 6.0 * cls.L**2 * x**2 + x**4)

    def __GetGuidedEndDisplacement(self, model_part):
        guided_end_node = min(
            (node for node in model_part.Nodes if node.X0 > 1e-9),
            key=lambda node: node.X0)
        uz = guided_end_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[2]
        exact_uz = self._ExactMidsurfaceDeflection(guided_end_node.X0, self.THICKNESS)
        return uz, exact_uz

    def __RunSBM(self, formulation_type):
        work_folder = "patch_test/shell_shifted_boundary"
        with open(f"{work_folder}/plate_guided_ss_sbm_parameters.json", 'r') as parameter_file:
            settings = KratosMultiphysics.Parameters(parameter_file.read())

        class ShellShiftedBoundaryAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):
            def ModifyInitialGeometry(self):
                super().ModifyInitialGeometry()
                for node in self._GetSolver().GetComputingModelPart().Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, node.X)

            def ApplyBoundaryConditions(self):
                super().ApplyBoundaryConditions()
                model_part = self._GetSolver().GetComputingModelPart()
                for condition in model_part.GetSubModelPart("shifted_boundary").Conditions:
                    condition.SetValue(KratosMultiphysics.DISPLACEMENT, [0.0, 0.0, 0.0])
                    condition.SetValue(KratosMultiphysics.ROTATION, [0.0, 0.0, 0.0])
                model_part.GetProperties(1).SetValue(SMA.SBM_FORMULATION_TYPE, formulation_type)
                model_part.GetProperties(1).SetValue(
                    SMA.SBM_CONSTRAINED_DOFS, KratosMultiphysics.Vector([1.0, 1.0, 0.0, 1.0, 1.0, 1.0]))

            def Initialize(self):
                super().Initialize()
                model_part = self._GetSolver().GetComputingModelPart()
                density = model_part.GetProperties(1).GetValue(KratosMultiphysics.DENSITY)
                accel_z = TestPatchTestShellShiftedBoundary.Q_LOAD / (density * TestPatchTestShellShiftedBoundary.THICKNESS)
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION, [0.0, 0.0, accel_z])

        model = KratosMultiphysics.Model()
        simulation = ShellShiftedBoundaryAnalysis(model, settings)
        simulation.Run()
        return self.__GetGuidedEndDisplacement(simulation._GetSolver().GetComputingModelPart())

    def __RunBodyFitted(self):
        work_folder = "patch_test/shell_shifted_boundary"
        with open(f"{work_folder}/plate_guided_ss_bodyfitted_parameters.json", 'r') as parameter_file:
            settings = KratosMultiphysics.Parameters(parameter_file.read())

        class ShellBodyFittedAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):
            def Initialize(self):
                super().Initialize()
                model_part = self._GetSolver().GetComputingModelPart()
                density = model_part.GetProperties(1).GetValue(KratosMultiphysics.DENSITY)
                accel_z = TestPatchTestShellShiftedBoundary.Q_LOAD / (density * TestPatchTestShellShiftedBoundary.THICKNESS)
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION, [0.0, 0.0, accel_z])

        model = KratosMultiphysics.Model()
        simulation = ShellBodyFittedAnalysis(model, settings)
        simulation.Run()
        # No shifted-boundary/ghost machinery here -- x=0 is a genuine mesh node.
        model_part = simulation._GetSolver().GetComputingModelPart()
        guided_end_node = min(model_part.Nodes, key=lambda node: node.X0)
        uz = guided_end_node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[2]
        exact_uz = self._ExactMidsurfaceDeflection(guided_end_node.X0, self.THICKNESS)
        return uz, exact_uz

    def __RunNitsche(self):
        work_folder = "patch_test/shell_shifted_boundary"
        with open(f"{work_folder}/plate_guided_ss_nitsche_parameters.json", 'r') as parameter_file:
            settings = KratosMultiphysics.Parameters(parameter_file.read())

        class ShellNitscheBoundaryAnalysis(structural_mechanics_analysis.StructuralMechanicsAnalysis):
            def ModifyInitialGeometry(self):
                super().ModifyInitialGeometry()
                for node in self._GetSolver().GetComputingModelPart().Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, node.X - 0.001)

            def ApplyBoundaryConditions(self):
                super().ApplyBoundaryConditions()
                model_part = self._GetSolver().GetComputingModelPart()
                model_part.GetProperties(1).SetValue(
                    SMA.SBM_CONSTRAINED_DOFS, KratosMultiphysics.Vector([1.0, 1.0, 0.0, 1.0, 1.0, 1.0]))

            def Initialize(self):
                super().Initialize()
                model_part = self._GetSolver().GetComputingModelPart()
                density = model_part.GetProperties(1).GetValue(KratosMultiphysics.DENSITY)
                accel_z = TestPatchTestShellShiftedBoundary.Q_LOAD / (density * TestPatchTestShellShiftedBoundary.THICKNESS)
                for node in model_part.Nodes:
                    node.SetSolutionStepValue(KratosMultiphysics.VOLUME_ACCELERATION, [0.0, 0.0, accel_z])

        model = KratosMultiphysics.Model()
        simulation = ShellNitscheBoundaryAnalysis(model, settings)
        simulation.Run()
        return self.__GetGuidedEndDisplacement(simulation._GetSolver().GetComputingModelPart())

    def __AssertRegression(self, uz, exact_uz, regression_value, sanity_rel_tol):
        self.assertAlmostEqual(uz, regression_value, delta=1e-7)
        self.assertLess(abs(uz - exact_uz), sanity_rel_tol * abs(exact_uz))

    def testShellShiftedBoundaryGuidedSimplySupportedBodyFitted(self):
        uz, exact_uz = self.__RunBodyFitted()
        self.__AssertRegression(uz, exact_uz, regression_value=0.0032283949, sanity_rel_tol=0.10)

    def testShellShiftedBoundaryGuidedSimplySupportedMLS(self):
        uz, exact_uz = self.__RunSBM(0)
        self.__AssertRegression(uz, exact_uz, regression_value=0.0027105547, sanity_rel_tol=0.30)

    def testShellShiftedBoundaryGuidedSimplySupportedTaylor(self):
        uz, exact_uz = self.__RunSBM(1)
        self.__AssertRegression(uz, exact_uz, regression_value=0.0026299507, sanity_rel_tol=0.30)

    def testShellShiftedBoundaryGuidedSimplySupportedGapSBM(self):
        uz, exact_uz = self.__RunSBM(2)
        self.__AssertRegression(uz, exact_uz, regression_value=0.0026257076, sanity_rel_tol=0.30)

    def testShellShiftedBoundaryGuidedSimplySupportedNitsche(self):
        uz, exact_uz = self.__RunNitsche()
        self.__AssertRegression(uz, exact_uz, regression_value=0.0029875002, sanity_rel_tol=0.10)

if __name__ == '__main__':
    # KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()