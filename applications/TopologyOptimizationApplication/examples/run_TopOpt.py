#import kratos core and applications
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
import KratosMultiphysics.TopologyOptimizationApplication as kto
import KratosMultiphysics.LinearSolversApplication as kls
from KratosMultiphysics import process_factory
import os
import OptimizationParameters as opt_parameters
from importlib import import_module
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import KratosMultiphysics.KratosUnittest as KratosUnittest

class RunTopOPt(KratosUnittest.TestCase):

    def _apply_material_properties(self, ModelPart, Dimesion):
         # Define material properties
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.YOUNG_MODULUS, 200.0e9)
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.POISSON_RATIO, 0.4)
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.DENSITY, 1.0)

        # Define body force
        g = [0,0,0]
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.VOLUME_ACCELERATION, g)

        # Define constitutive law
        ModelPart.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, cons_law)

    def _solve(self, ModelPart):
        # Define a linear strategy to solve the problem
        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True

        strategy = KratosMultiphysics.ResidualBasedLinearStrategy(
            ModelPart,
            scheme,
            builder_and_solver,
            compute_reactions,
            reform_step_dofs,
            calculate_norm_dx,
            move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Check()

        # Solve the problem
        strategy.Solve()

    def RunTopOPt(self):
        parameter_file = open("/home/philipp/opt/kratosDev/applications/TopologyOptimizationApplication/examples/01_Small_Cantilever_Hexahedra/ProjectParameters.json",'r')
        ProjectParameters = km.Parameters(parameter_file.read())
        echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
        dimension = 2
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Structure")
        self._add_variables(model_part)
        self.AddVariables(model_part)
        self.AddDofs(model_part)
        self.ImportModelPart(model_part)
        self._apply_material_properties(model_part, dimension)

        self._solve(model_part)

    def __post_process(self, main_model_part, post_type = "gid"):
        if post_type == "gid":
            self.gid_output = GiDOutputProcess(
                main_model_part,
                main_model_part.Name,
                KratosMultiphysics.Parameters(r"""
                {
                    "result_file_configuration" : {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "nodal_results"       : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
                    "gauss_point_results" : ["CAUCHY_STRESS_VECTOR"]
                    }
                }"""))

            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()
            self.gid_output.ExecuteInitializeSolutionStep()
            self.gid_output.PrintOutput()
            self.gid_output.ExecuteFinalizeSolutionStep()
            self.gid_output.ExecuteFinalize()
if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()