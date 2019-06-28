import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.CoSimulationApplication.co_simulation_interface import CoSimulationInterface
from KratosMultiphysics.CoSimulationApplication.co_simulation_tools import ImportDataStructure
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools


class TestConvergenceCriterionRelativeNorm(KratosUnittest.TestCase):
    def test_convergence_criterion_relative_norm(self):
        parameter_file_name = "test_parameters.json"
        cs_data_structure = ImportDataStructure(parameter_file_name)

        m = 10
        dz = 2.0
        a0 = 10.0
        a1 = 1.0e-4
        a2 = 1.0e-6
        interface_settings = {"wall": "AREA"}

        # Create interface
        variable = "AREA"
        model = cs_data_structure.Model()
        model_part = model.CreateModelPart("wall")
        model_part.AddNodalSolutionStepVariable(variable)
        for i in range(m):
            model_part.CreateNewNode(i, 0.0, 0.0, i * dz)
        step = 0
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, step, a0)
        interface = CoSimulationInterface(model, interface_settings)

        parameter_file_name = "convergence_criteria/test_relative_norm.json"
        with open(parameter_file_name, 'r') as parameter_file:
            settings = cs_data_structure.Parameters(parameter_file.read())

        convergence_criterion_relative_norm = cs_tools.CreateInstance(settings)
        convergence_criterion_relative_norm.Initialize()
        for i in range(3):
            convergence_criterion_relative_norm.InitializeSolutionStep()
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a0)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a1)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a2)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertTrue(is_satisfied)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(variable, step, a1)
            convergence_criterion_relative_norm.Update(interface)
            is_satisfied = convergence_criterion_relative_norm.IsSatisfied()
            self.assertFalse(is_satisfied)
            convergence_criterion_relative_norm.FinalizeSolutionStep()


if __name__ == '__main__':
    KratosUnittest.main()
