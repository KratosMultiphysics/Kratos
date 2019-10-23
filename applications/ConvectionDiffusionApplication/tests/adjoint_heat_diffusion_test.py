import KratosMultiphysics as kratos
import KratosMultiphysics.ConvectionDiffusionApplication as convdiff

import KratosMultiphysics.KratosUnittest as unittest

from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting

class AdjointHeatDiffusionTest(unittest.TestCase):

    def setUp(self):
        self.input_file = "adjoint_test"
        self.work_folder = "adjoint_diffusion_test"
        self.primal_parameter_file_name = "ProjectParametersPrimal.json"
        self.adjoint_parameter_file_name = "ProjectParametersAdjoint.json"

        self.check_tolerance = 1e-4
        self.perturbation_magnitude = 1e-7
        self.print_output = False
        self.node_ids_to_test = range(1,10)

        self.volume_source = 0.0

    def tearDown(self):

        with unittest.WorkFolderScope(self.work_folder, __file__):
            DeleteFileIfExisting(self.input_file+"_x.mdpa")
            DeleteFileIfExisting(self.input_file+"_y.mdpa")
            if not self.print_output:
                DeleteFileIfExisting("adjoint_diffusion_test.post.lst")
                DeleteFileIfExisting("diffusion_test_primal.post.bin")
                DeleteFileIfExisting("diffusion_test_adjoint.post.bin")

    def testAdjointHeatDiffusion(self):
        self.directSensitivityCheck()

    def testAdjointHeatDiffusionWithSourceTerm(self):
        self.volume_source = 200.0
        self.directSensitivityCheck()

    def testAdjointHeatDiffusionWithConvectionCondition(self):
        self.primal_parameter_file_name = "ProjectParametersConvectionConditionPrimal.json"
        self.adjoint_parameter_file_name = "ProjectParametersConvectionConditionAdjoint.json"
        self.directSensitivityCheck()

    def testAdjointHeatDiffusionWithRadiation(self):
        self.primal_parameter_file_name = "ProjectParametersRadiationPrimal.json"
        self.adjoint_parameter_file_name = "ProjectParametersRadiationAdjoint.json"
        self.directSensitivityCheck()

    def directSensitivityCheck(self):
        model = kratos.Model()
        settings = kratos.Parameters(r'''{}''')

        with unittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.primal_parameter_file_name,'r') as primal_parameter_file:
                settings.AddValue("primal_settings", kratos.Parameters(primal_parameter_file.read()))

            # enable volume source if needed
            other_processes_list = settings["primal_settings"]["processes"]["list_other_processes"]
            for i in range(other_processes_list.size()):
                if other_processes_list[i]["process_name"].GetString() == "AssignScalarVariableProcess" and other_processes_list[i]["Parameters"]["variable_name"].GetString() == "HEAT_FLUX":
                    other_processes_list[i]["Parameters"]["value"].SetDouble(self.volume_source)

            self.primalSolution(model, settings["primal_settings"])

            with open(self.adjoint_parameter_file_name,'r') as adjoint_parameter_file:
                settings.AddValue("adjoint_settings", kratos.Parameters(adjoint_parameter_file.read()))

            adjoint_analysis = ConvectionDiffusionAnalysis(model,settings["adjoint_settings"])
            adjoint_analysis.Run()

            objective_model_part_name = settings["adjoint_settings"]["solver_settings"]["response_function_settings"]["custom_settings"]["model_part_name"].GetString()
            finite_difference_results = self.perturbedSolution(model, settings["primal_settings"], objective_model_part_name, self.node_ids_to_test)

            adjoint_model_part = model.GetModelPart(settings["adjoint_settings"]["solver_settings"]["model_part_name"].GetString())
            for node_id,fd_sensitivity in zip(self.node_ids_to_test, finite_difference_results):
                node = adjoint_model_part.Nodes[node_id]
                adjoint_sensitivity = node.GetSolutionStepValue(kratos.SHAPE_SENSITIVITY,0)
                #diff = adjoint_sensitivity-fd_sensitivity
                #if (diff[0]**2 + diff[1]**2 + diff[2]**2)**0.5 < self.check_tolerance:
                #    status = "PASS"
                #else:
                #    status = "FAIL"
                #print(status, node_id, adjoint_sensitivity, fd_sensitivity, diff)
                self.assertAlmostEqual(adjoint_sensitivity[0], fd_sensitivity[0], delta=self.check_tolerance)
                self.assertAlmostEqual(adjoint_sensitivity[1], fd_sensitivity[1], delta=self.check_tolerance)

    def primalSolution(self, model, settings):

        primal_analysis = ConvectionDiffusionAnalysis(model,settings)
        primal_analysis.Run()

    def perturbedSolution(self, model, settings, objective_model_part_name, perturbed_node_ids):
        model_part_name = settings["solver_settings"]["model_part_name"].GetString()
        input_file_name = settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

        base_temp = self.average_temperature_objective(model.GetModelPart(model_part_name+"."+objective_model_part_name))
        results = []

        for node_id in perturbed_node_ids:
            coord = self._readNodalCoordinates(node_id,input_file_name)
            # X perturbation
            perturbed_coord = [coord[0]+self.perturbation_magnitude, coord[1], coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,input_file_name,input_file_name+"_x")
            model_x = kratos.Model()
            settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(input_file_name+"_x")
            self.primalSolution(model_x,settings)
            x_temp = self.average_temperature_objective(model_x.GetModelPart(model_part_name+"."+objective_model_part_name))

			# Y perturbation
            perturbed_coord = [coord[0], coord[1]+self.perturbation_magnitude, coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,input_file_name,input_file_name+"_y")
            model_y = kratos.Model()
            settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(input_file_name+"_y")
            self.primalSolution(model_y,settings)
            y_temp = self.average_temperature_objective(model_y.GetModelPart(model_part_name+"."+objective_model_part_name))

            result = kratos.Array3()
            result[0] = (x_temp-base_temp)/self.perturbation_magnitude
            result[1] = (y_temp-base_temp)/self.perturbation_magnitude

            results.append(result)

        return results

    def average_temperature_objective(self, model_part):
        temp = 0.0
        for node in model_part.Nodes:
            temp += node.GetSolutionStepValue(kratos.TEMPERATURE)
        return temp / len(model_part.Nodes)

    # the following two functions are taken from applications/FluidDynamicsApplication/tests/adjoint_vms_sensitivity_2d.py
    # Maybe they can be moved to a common utility?
    def _readNodalCoordinates(self,node_id,model_part_file_name):
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        line = lines[node_id] # assumes consecutive node numbering starting with 1
        components = line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        return [float(components[i]) for i in range(1,4)]

    def _writeNodalCoordinates(self,node_id,coords,model_part_file_name,output_file_name=None):
        if output_file_name is None:
            output_file_name = model_part_file_name
        with open(model_part_file_name + '.mdpa', 'r') as model_part_file:
            lines = model_part_file.readlines()
        node_lines = lines[lines.index('Begin Nodes\n'):lines.index('End Nodes\n')]
        old_line = node_lines[node_id] # assumes consecutive node numbering starting with 1
        components = old_line.split()
        if int(components[0]) != node_id:
            raise RuntimeError('Error parsing file ' + model_part_file_name)
        new_line = '{:5d}'.format(node_id) + ' ' \
             + '{:19.10f}'.format(coords[0]) + ' ' \
             + '{:19.10f}'.format(coords[1]) + ' ' \
             + '{:19.10f}'.format(coords[2]) + '\n'
        lines[lines.index(old_line)] = new_line
        with open(output_file_name + '.mdpa', 'w') as model_part_file:
            model_part_file.writelines(lines)


if __name__ == '__main__':
    unittest.main()