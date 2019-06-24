import KratosMultiphysics as kratos
import KratosMultiphysics.ConvectionDiffusionApplication as convdiff

import KratosMultiphysics.KratosUnittest as unittest

from convection_diffusion_analysis import ConvectionDiffusionAnalysis

class AdjointHeatDiffusionTest(unittest.TestCase):

    def setUp(self):
        self.input_file = "adjoint_test"
        #self.reference_file = "reference10_qasgs"
        self.work_folder = "adjoint_diffusion_test"
        self.primal_parameter_file_name = "ProjectParametersPrimal.json"
        self.adjoint_parameter_file_name = "ProjectParametersAdjoint.json"

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

    def testAdjointHeatDiffusion(self):
        model = kratos.Model()
        settings = kratos.Parameters(r'''{}''')

        with unittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.primal_parameter_file_name,'r') as primal_parameter_file:
                settings.AddValue("primal_settings", kratos.Parameters(primal_parameter_file.read()))

            self.primalSolution(model, settings["primal_settings"])

            with open(self.adjoint_parameter_file_name,'r') as adjoint_parameter_file:
                settings.AddValue("adjoint_settings", kratos.Parameters(adjoint_parameter_file.read()))

            adjoint_analysis = ConvectionDiffusionAnalysis(model,settings["adjoint_settings"])
            adjoint_analysis.Run()

            self.perturbedSolution(model, settings["primal_settings"],[8],1e-6)

    def primalSolution(self, model, settings):

        primal_analysis = ConvectionDiffusionAnalysis(model,settings)
        primal_analysis.Run()

    def perturbedSolution(self, model, settings, perturbed_node_ids, perturbation_magnitude):
        model_part_name = settings["solver_settings"]["model_part_name"].GetString()
        input_file_name = settings["solver_settings"]["model_import_settings"]["input_filename"].GetString()

        base_temp = model.GetModelPart(model_part_name).Nodes[perturbed_node_ids[0]].GetSolutionStepValue(kratos.TEMPERATURE)

        for node_id in perturbed_node_ids:
            coord = self._readNodalCoordinates(node_id,input_file_name)
            # X perturbation
            perturbed_coord = [coord[0]+perturbation_magnitude, coord[1], coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,input_file_name,input_file_name+"_x")
            model_x = kratos.Model()
            settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(input_file_name+"_x")
            self.primalSolution(model_x,settings)
            x_temp = model_x.GetModelPart(model_part_name).Nodes[perturbed_node_ids[0]].GetSolutionStepValue(kratos.TEMPERATURE)
            # Y perturbation
            perturbed_coord = [coord[0], coord[1]+perturbation_magnitude, coord[2]]
            self._writeNodalCoordinates(node_id,perturbed_coord,input_file_name,input_file_name+"_y")
            model_y = kratos.Model()
            settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(input_file_name+"_y")
            self.primalSolution(model_y,settings)
            y_temp = model_y.GetModelPart(model_part_name).Nodes[perturbed_node_ids[0]].GetSolutionStepValue(kratos.TEMPERATURE)

        print(base_temp/perturbation_magnitude)
        print(x_temp/perturbation_magnitude)
        print(y_temp/perturbation_magnitude)
        print([(x_temp-base_temp)/perturbation_magnitude, (y_temp-base_temp)/perturbation_magnitude])


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