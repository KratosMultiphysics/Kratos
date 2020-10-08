from __future__ import print_function, absolute_import, division
import csv

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.MORApplication as MOR

from KratosMultiphysics.python_linear_solver_factory import ConstructSolver

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
linear_solvers_application_available = kratos_utils.CheckIfApplicationsAvailable("LinearSolversApplication")

class TestIrkaStrategies(KratosUnittest.TestCase):

    def _add_variables(self,mp):
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        mp.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        mp.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        mp.AddNodalSolutionStepVariable(MOR.COMPONENT_OUTPUT)

    def _add_dofs(self,mp):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,mp)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,mp)

    def _create_plate_model(self, material_settings):
        model = KratosMultiphysics.Model()
        mp = model.CreateModelPart('mp')
        self._add_variables(mp)

        input_file = 'test_input/cantilever_plate'

        KratosMultiphysics.ModelPartIO(input_file).ReadModelPart(mp)
        KratosMultiphysics.ReadMaterialsUtility(material_settings, model)
        self._add_dofs(mp)

        mp.SetBufferSize(2)
        mp.Check(mp.ProcessInfo)

        # add load
        for node in model.GetModelPart('mp.PointLoad3D_tipright').Nodes:
            node.SetSolutionStepValue(StructuralMechanicsApplication.POINT_LOAD,0,[0,0,-1000])

        for node in model.GetModelPart('mp.DISPLACEMENT_left').Nodes:
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
            node.Fix(KratosMultiphysics.ROTATION_X)
            node.Fix(KratosMultiphysics.ROTATION_Y)
            node.Fix(KratosMultiphysics.ROTATION_Z)

        return model

    @KratosUnittest.skipUnless(linear_solvers_application_available,"Missing required application: LinearSolversApplication")
    def test_complex_irka(self):
        import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication

        material_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "materials_filename" : "test_input/StructuralMaterials_damped_modal.json"
            }
        }
        """)
        model = self._create_plate_model(material_settings)
        mp = model.GetModelPart("mp")

        output_process_settings = KratosMultiphysics.Parameters("""
        {
            "build_output_structure": true,
            "output_structure_type": "vector",
            "model_part_names" : ["GENERIC_right"],
            "output_variable_names" : ["DISPLACEMENT_Z"]
        }
        """)

        solver_type = 'skyline_lu_complex'
        cplx_linear_solver_settings = KratosMultiphysics.Parameters('{ "solver_type" : "' + solver_type + '" }')
        cplx_linear_solver = ConstructSolver(cplx_linear_solver_settings)

        solver_type = 'complex_dense_col_piv_householder_qr'
        cplx_dense_linear_solver_settings = KratosMultiphysics.Parameters('{ "solver_type" : "LinearSolversApplication.' + solver_type + '" }')
        cplx_dense_linear_solver = LinearSolversApplication.ComplexDenseLinearSolverFactory().Create(cplx_dense_linear_solver_settings)

        scheme = MOR.MatrixBuilderScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
        off_strategy = MOR.MorSecondOrderComplexIrkaStrategy(mp,
            scheme,
            builder_and_solver,
            cplx_linear_solver,
            [15j, 25j, 45j],
            20,
            1e-10)
        off_strategy.SetEchoLevel(0)

        strategy = MOR.MorComplexOnlineStrategy(mp,
            cplx_dense_linear_solver,
            off_strategy,
            False)

        freq = 10
        max_freq = 30
        df = 10

        result = []
        n_nodes = len(model.GetModelPart('mp.GENERIC_right').Nodes)
        for cond in model.GetModelPart('mp.GENERIC_right').Conditions:
            cond.SetValue(MOR.COMPONENT_OUTPUT,[0,0,1])

        while freq <= max_freq:
            mp.CloneTimeStep(freq)
            mp.ProcessInfo[MOR.FREQUENCY] = freq

            strategy.Solve()

            disp_z = strategy.GetScalarResult()
            result.append([freq, disp_z/n_nodes])

            freq = freq + df

        # with open('test_results/irka_complex_result.csv', 'w', newline='') as csvfile:
        #     writer = csv.writer(csvfile)
        #     writer.writerows(result)
        with open('test_results/irka_complex_result.csv', newline='') as csvfile:
            expected_result = csv.reader(csvfile)
            for r, xr in zip(result, expected_result):
                self.assertAlmostEqual(r[1], complex(xr[1]))

    @KratosUnittest.skipUnless(linear_solvers_application_available,"Missing required application: LinearSolversApplication")
    def test_real_irka(self):
        import KratosMultiphysics.LinearSolversApplication as LinearSolversApplication
        material_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                "materials_filename" : "test_input/StructuralMaterials_damped.json"
            }
        }
        """)
        model = self._create_plate_model(material_settings)
        mp = model.GetModelPart("mp")

        output_process_settings = KratosMultiphysics.Parameters("""
        {
            "build_output_structure": true,
            "output_structure_type": "vector",
            "model_part_names" : ["GENERIC_right"],
            "output_variable_names" : ["DISPLACEMENT_Z"]
        }
        """)

        solver_type = 'skyline_lu_complex'
        cplx_linear_solver_settings = KratosMultiphysics.Parameters('{ "solver_type" : "' + solver_type + '" }')
        cplx_linear_solver = ConstructSolver(cplx_linear_solver_settings)

        solver_type = 'complex_dense_col_piv_householder_qr'
        cplx_dense_linear_solver_settings = KratosMultiphysics.Parameters('{ "solver_type" : "LinearSolversApplication.' + solver_type + '" }')
        cplx_dense_linear_solver = LinearSolversApplication.ComplexDenseLinearSolverFactory().Create(cplx_dense_linear_solver_settings)

        scheme = MOR.MatrixBuilderScheme()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(KratosMultiphysics.LinearSolver())
        off_strategy = MOR.MorSecondOrderRealIrkaStrategy(mp,
            scheme,
            builder_and_solver,
            cplx_linear_solver,
            [15j, 25j, 45j],
            20,
            1e-12)
        off_strategy.SetEchoLevel(0)

        strategy = MOR.MorRealOnlineStrategy(mp,
            cplx_dense_linear_solver,
            off_strategy,
            False)

        freq = 10
        max_freq = 30
        df = 10

        result = []
        n_nodes = len(model.GetModelPart('mp.GENERIC_right').Nodes)
        for cond in model.GetModelPart('mp.GENERIC_right').Conditions:
            cond.SetValue(MOR.COMPONENT_OUTPUT,[0,0,1])

        while freq <= max_freq:
            mp.CloneTimeStep(freq)
            mp.ProcessInfo[MOR.FREQUENCY] = freq

            strategy.Solve()

            disp_z = strategy.GetScalarResult()
            result.append([freq, disp_z/n_nodes])

            freq = freq + df

        # with open('test_results/irka_real_result.csv', 'w', newline='') as csvfile:
        #     writer = csv.writer(csvfile)
        #     writer.writerows(result)
        with open('test_results/irka_real_result.csv', newline='') as csvfile:
            expected_result = csv.reader(csvfile)
            for r, xr in zip(result, expected_result):
                self.assertAlmostEqual(r[1], complex(xr[1]))

if __name__ == '__main__':
    KratosUnittest.main()
