from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_analysis

external_solvers_application_available = kratos_utils.CheckIfApplicationsAvailable("ExternalSolversApplication")
eigen_solvers_application_available = kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication")
hdf5_application_available = kratos_utils.CheckIfApplicationsAvailable("HDF5Application")

from math import sqrt
from cmath import phase

class HarmonicAnalysisTests(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)

    def _solve_eigen(self,mp,echo=0):
        feast_system_solver_settings = KratosMultiphysics.Parameters("""{ }""")

        eigensolver_settings = KratosMultiphysics.Parameters("""
        {
                "perform_stochastic_estimate": false,
                "lambda_min": 1.0e-2,
                "lambda_max": 25.0,
                "search_dimension": 7
        }
        """)
        if not (hasattr(KratosMultiphysics.ExternalSolversApplication,"PastixComplexSolver")):
            self.skipTest('"PastixComplexSolver" not available')

        feast_system_solver = ExternalSolversApplication.PastixComplexSolver(feast_system_solver_settings)
        eigen_solver = ExternalSolversApplication.FEASTSolver(eigensolver_settings, feast_system_solver)
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(eigen_solver)

        eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        eig_strategy = StructuralMechanicsApplication.EigensolverStrategy(mp,
                                                                    eigen_scheme,
                                                                    builder_and_solver)

        eig_strategy.SetEchoLevel(echo)
        eig_strategy.Solve()

    def _setup_harmonic_solver(self,mp,echo=0):
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
        eigen_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        harmonic_strategy = StructuralMechanicsApplication.HarmonicAnalysisStrategy(mp, eigen_scheme, builder_and_solver, False)
        harmonic_strategy.SetEchoLevel(echo)

        return harmonic_strategy


    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)

    def _create_2dof_geometry(self, current_model, stiffness, mass, damping=0):
        mp = current_model.CreateModelPart("mdof")
        self._add_variables(mp)

        base = mp.CreateNewNode(3,0.0,0.0,0.0)
        node1 = mp.CreateNewNode(1,10.0,0.0,0.0)
        node2 = mp.CreateNewNode(2,20.0,0.0,0.0)

        self._add_dofs(mp)
        for no in mp.Nodes:
            no.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            no.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        base.Fix(KratosMultiphysics.DISPLACEMENT_X)

        #create elements and conditions
        spring1 = mp.CreateNewElement("SpringDamperElement3D2N", 1, [3,1], mp.GetProperties()[1])
        spring2 = mp.CreateNewElement("SpringDamperElement3D2N", 2, [2,1], mp.GetProperties()[1])
        mass1 = mp.CreateNewElement("NodalConcentratedElement3D1N", 3, [1], mp.GetProperties()[1])
        mass2 = mp.CreateNewElement("NodalConcentratedElement3D1N", 4, [2], mp.GetProperties()[1])
        mp.CreateNewCondition("PointLoadCondition3D1N",1,[1],mp.GetProperties()[1])
        mp.CreateNewCondition("PointLoadCondition3D1N",2,[2],mp.GetProperties()[1])

        mass1.SetValue(KratosMultiphysics.NODAL_MASS,mass)
        mass2.SetValue(KratosMultiphysics.NODAL_MASS,mass/2)
        spring1.SetValue(StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS,[stiffness,0,0])
        spring2.SetValue(StructuralMechanicsApplication.NODAL_DISPLACEMENT_STIFFNESS,[stiffness/2,0,0])
        node1.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[1,0,0])
        mp.GetProperties()[1].SetValue(StructuralMechanicsApplication.SYSTEM_DAMPING_RATIO, damping)

        return mp

    @KratosUnittest.skipUnless(external_solvers_application_available,"Missing required application: ExternalSolversApplication")
    def test_undamped_mdof_harmonic(self):
        import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
        current_model = KratosMultiphysics.Model()
        #analytic solution taken from Humar - Dynamics of Structures p. 675

        #material properties
        stiffness = 10.0
        mass = 2.0

        #create the model
        mp = self._create_2dof_geometry(current_model, stiffness, mass)

        #solve the eigenproblem
        self._solve_eigen(mp)

        #perform the harmonic analysis
        harmonic_solver = self._setup_harmonic_solver(mp)

        exfreq = 1.0
        max_exfreq = 2.0
        df = 0.05

        while(exfreq <= max_exfreq):
            mp.CloneTimeStep(exfreq)
            harmonic_solver.Solve()

            disp_x1_expected = ((1.0/3.0 * stiffness) / (0.5 - (exfreq / sqrt(stiffness/mass))**2) \
                + (1.0/1.5 * stiffness) / (2.0 - (exfreq / sqrt(stiffness/mass))**2)) / stiffness**2
            disp_x2_expected = ((2.0/3.0 * stiffness) / (0.5 - (exfreq / sqrt(stiffness/mass))**2) \
                - (2.0/3.0 * stiffness) / (2.0 - (exfreq / sqrt(stiffness/mass))**2)) / stiffness**2

            self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0), \
                disp_x1_expected,delta=1e-5)
            self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0), \
                disp_x2_expected,delta=1e-5)

            exfreq = exfreq + df

    @KratosUnittest.skipUnless(external_solvers_application_available,"Missing required application: ExternalSolversApplication")
    def test_damped_mdof_harmonic(self):
        import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
        current_model = KratosMultiphysics.Model()

        #analytic solution taken from Humar - Dynamics of Structures p. 677

        #material properties
        stiffness = 10.0
        mass = 2.0
        damping = 0.1

        #create the model
        mp = self._create_2dof_geometry(current_model,stiffness, mass, damping)

        #solve the eigenproblem
        self._solve_eigen(mp)

        #perform the harmonic analysis
        harmonic_solver = self._setup_harmonic_solver(mp)

        exfreq = 1.0
        max_exfreq = 2.0
        df = 0.05

        while(exfreq <= max_exfreq):
            mp.CloneTimeStep(exfreq)
            harmonic_solver.Solve()

            exfreq_d = exfreq / sqrt(stiffness/mass)

            disp_x1_expected_complex = 1 / (3*stiffness * complex(0.5 - exfreq_d**2, sqrt(2) * damping * exfreq_d)) \
                + 2 / (3*stiffness * complex(2 - exfreq_d**2, 2 * sqrt(2) * damping * exfreq_d))
            if disp_x1_expected_complex.real < 0:
                disp_x1_expected = - abs(disp_x1_expected_complex)
            else:
                disp_x1_expected = abs(disp_x1_expected_complex)

            disp_x2_expected_complex = 2 / (3*stiffness * complex(0.5 - exfreq_d**2, sqrt(2) * damping * exfreq_d)) \
                - 2 / (3*stiffness * complex(2 - exfreq_d**2, 2 * sqrt(2) * damping * exfreq_d))
            if disp_x2_expected_complex.real < 0:
                disp_x2_expected = - abs(disp_x2_expected_complex)
            else:
                disp_x2_expected = abs(disp_x2_expected_complex)

            #test displacement
            self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0), \
                disp_x1_expected,delta=1e-5)
            self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0), \
                disp_x2_expected,delta=1e-5)

            #test phase angle
            self.assertAlmostEqual(mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.REACTION_X,0), \
                phase(disp_x1_expected_complex),delta=1e-5)
            self.assertAlmostEqual(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.REACTION_X,0), \
                phase(disp_x2_expected_complex),delta=1e-5)

            exfreq = exfreq + df

class HarmonicAnalysisTestsWithHDF5(KratosUnittest.TestCase):
    @KratosUnittest.skipUnless(hdf5_application_available and eigen_solvers_application_available,"Missing required application: HDF5Application, EigenSolversApplication")
    def test_harmonic_mdpa_input(self):

        with KratosUnittest.WorkFolderScope(".",__file__):
            #run simulation and write to hdf5 file
            model = KratosMultiphysics.Model()
            project_parameter_file_name = "harmonic_analysis_test/harmonic_analysis_test_eigenproblem_parameters.json"
            with open(project_parameter_file_name,'r') as parameter_file:
                project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
            test = structural_mechanics_analysis.StructuralMechanicsAnalysis(model,project_parameters)
            test.Run()

            #start new simulation and read from hdf5 file
            model = KratosMultiphysics.Model()
            project_parameter_file_name = "harmonic_analysis_test/harmonic_analysis_test_parameters.json"
            with open(project_parameter_file_name,'r') as parameter_file:
                project_parameters = KratosMultiphysics.Parameters(parameter_file.read())
            test = structural_mechanics_analysis.StructuralMechanicsAnalysis(model,project_parameters)
            test.Run()

            # remove hdf5 file
            kratos_utils.DeleteFileIfExisting("harmonic_analysis_test/eigen_results.h5")
            kratos_utils.DeleteFileIfExisting("harmonic_analysis_test/harmonic_analysis_test.time")

if __name__ == '__main__':
    KratosUnittest.main()
