from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.MORApplication as MOR
import KratosMultiphysics.LinearSolversApplication as LSA
from KratosMultiphysics import eigen_solver_factory

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestAcousticElements(KratosUnittest.TestCase):

    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication","StructuralMechanicsApplication")
    @KratosUnittest.skipUnless(LSA.HasFEAST(),"FEAST not found in LinearSolversApplication, skipping.")
    def test_acoustic_2d4n_element(self):
        Model = KratosMultiphysics.Model()
        mp = Model.CreateModelPart('mp')
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        input_file = 'test_input/channel_2d'
        material_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "materials_filename" : "test_input/AcousticMaterials.json"
            }
        }
        """)

        KratosMultiphysics.ModelPartIO(input_file).ReadModelPart(mp)
        KratosMultiphysics.ReadMaterialsUtility(material_settings, Model)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE,mp)

        mp.SetBufferSize(2)
        mp.Check(mp.ProcessInfo)

        eigensolver_settings = KratosMultiphysics.Parameters("""{
                "solver_type": "eigen_feast",
                "echo_level": 0,
                "search_lowest_eigenvalues": false,
                "e_min" : 1.0e3,
                "e_max" : 5.0e7,
                "subspace_size" : 20
            }""")

        linear_solver = LSA.FEASTSymmetricEigensystemSolver(eigensolver_settings)
        scheme = SMA.EigensolverDynamicScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        strategy = SMA.EigensolverStrategy(mp, scheme, builder_and_solver, 1.0, -1.0)

        strategy.Check()
        strategy.Initialize()

        strategy.Solve()

        eigen_vals = mp.ProcessInfo[SMA.EIGENVALUE_VECTOR]
        eigen_vals_exp = [1141020.1092768596,4565206.671350432,10275939.498756096,18278855.305498056,28581853.224728845,
            28581853.22472898,29722873.334005963,33147059.896079566,38857792.72348527,41195102.48913678,46860708.53022728]
        for ev, ev_exp in zip(eigen_vals, eigen_vals_exp):
            self.assertAlmostEqual(ev, ev_exp)


    @KratosUnittest.skipIfApplicationsNotAvailable("LinearSolversApplication","StructuralMechanicsApplication")
    @KratosUnittest.skipUnless(LSA.HasFEAST(),"FEAST not found in LinearSolversApplication, skipping.")
    def test_acoustic_3d8n_element(self):
        Model = KratosMultiphysics.Model()
        mp = Model.CreateModelPart('mp')
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        input_file = 'test_input/channel_3d'
        material_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "materials_filename" : "test_input/AcousticMaterials.json"
            }
        }
        """)

        KratosMultiphysics.ModelPartIO(input_file).ReadModelPart(mp)
        KratosMultiphysics.ReadMaterialsUtility(material_settings, Model)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE,mp)

        mp.SetBufferSize(2)
        mp.Check(mp.ProcessInfo)

        eigensolver_settings = KratosMultiphysics.Parameters("""{
                "solver_type": "eigen_feast",
                "echo_level": 0,
                "search_lowest_eigenvalues": false,
                "e_min" : 1.0e3,
                "e_max" : 1.0e7,
                "subspace_size": 5
            }""")

        linear_solver = LSA.FEASTSymmetricEigensystemSolver(eigensolver_settings)
        scheme = SMA.EigensolverDynamicScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        strategy = SMA.EigensolverStrategy(mp, scheme, builder_and_solver, 1.0, -1.0)

        strategy.Check()
        strategy.Initialize()

        strategy.Solve()

        eigen_vals = mp.ProcessInfo[SMA.EIGENVALUE_VECTOR]
        eigen_vals_exp = [1141301.6678376605,4569713.826374584]
        for ev, ev_exp in zip(eigen_vals, eigen_vals_exp):
            self.assertAlmostEqual(ev, ev_exp)

if __name__ == '__main__':
    KratosUnittest.main()
