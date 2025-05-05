import os
import tempfile
import shutil
import numpy as np

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as kratos_geo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.test_helper as test_helper

class prob_analysis:
    def __init__(self, template_project_path, input_material_parameters, output_parameters):
        """
        Create and initialize model - all construction stages

        :param template_project_path:
        :param input_parameters:
        :param output_parameters:

        :return:
        """

        self.input_material_parameters = input_material_parameters
        self.output_parameters = output_parameters

        # record working current folder
        self.starting_cwd = os.getcwd()

        # create a tmp folder
        self.tmp_folder = tempfile.mkdtemp()
        os.chdir(self.tmp_folder)

        # copy contents of template_project_path to tmp_folder and change working directory
        shutil.copytree(template_project_path, self.tmp_folder, dirs_exist_ok=True)

        # get names of files with template 'ProjectParameters_stageX.json' where X is the stage number and is automatically found in a sorted list
        parameter_file_names = sorted([f for f in os.listdir(template_project_path)
                                       if f.startswith('ProjectParameters_stage') and f.endswith('.json')])

        # set stage parameters
        parameters_stages = []
        for parameter_file_name in parameter_file_names:
            with open(parameter_file_name, 'r') as parameter_file:
                parameters_stages.append(Kratos.Parameters(parameter_file.read()))

        self.model = Kratos.Model()
        self.parameters_stages = parameters_stages

    def update_material_parameters(self, input_parameter, input_value):
        # update the parameters in the materials on the model
        # input_parameter structure is stage number, model part, material parameter name
        # e.g. input_parameter = [0, "Soil_1", "YOUNG_MODULUS"]
        model_part = self.model.GetModelPart(input_parameter[1])
        list_expected = False

        if len(input_parameter) == 4:
            list_expected = True

        for ind, element in enumerate(model_part.Elements):
            properties = element.Properties
            if ind == 0 and list_expected:
                vector_values = properties[input_parameter[2]]
                vector_values[input_parameter[3]] = input_value
                properties[input_parameter[2]] = vector_values
            elif ind==0:
                properties[input_parameter[2]] = input_value
            element.Properties = properties


    def calculate(self, input_values):
        """
        Calculate the output values for the given input values
        :param input_values:
        :return:
        """

        output_values = []
        for stage_no, stage_parameters in enumerate(self.parameters_stages):
            stage_obj = analysis.GeoMechanicsAnalysis(self.model, stage_parameters)
            stage_obj.Initialize()

            for i, input_parameter in enumerate(self.input_material_parameters):
                if stage_no < input_parameter[0]:
                    print("Skipping input parameter for stage", stage_no)
                    continue
                value = input_values[i]
                self.update_material_parameters(input_parameter, value)

            stage_obj.RunSolutionLoop()
            stage_obj.Finalize()

            for output_parameter in self.output_parameters:
                if stage_no != output_parameter[0]:
                    continue
                model_part = self.model.GetModelPart(output_parameter[3])
                values = output_parameter[2](model_part, *output_parameter[4:])
                if output_parameter[1] is not None:
                    output_values.append(output_parameter[1](values))
                else:
                    output_values.append(values)
        return output_values

    def finalize(self):
        """
        Finalize the analysis and clean up
        :return:
        """
        # remove the tmp folder
        os.chdir(self.starting_cwd)
        shutil.rmtree(self.tmp_folder)

if __name__ == "__main__":

    # Example usage
    input_young_moduli = [100E5, 500E5, 900.0E05, 1000.0E05, 2000.0E05]
    min_displacements = []
    displacements_at_node_13 = []
    for young_modulus in input_young_moduli:
        template_project_path = r"C:\tmp\FEA-Tools\Quay_Wall\Quay_Wall_4Stage_quadratic_master_slave_interface_with_anchor_reset_displacement_and_prestress_truss_anchor"
        input_parameters = [[0, "PorousDomain.Parts_Solid_layer_1|1", Kratos.YOUNG_MODULUS],
                            [0, "PorousDomain.Parts_Solid_layer_1|2", Kratos.YOUNG_MODULUS],
                            [0, "PorousDomain.Parts_Solid_layer_2|1", Kratos.YOUNG_MODULUS],
                            [0, "PorousDomain.Parts_Solid_layer_2|2", Kratos.YOUNG_MODULUS],
                            [0, "PorousDomain.Parts_Solid_layer_3|1", Kratos.YOUNG_MODULUS],
                            [0, "PorousDomain.Parts_Solid_layer_3|2", Kratos.YOUNG_MODULUS],
                            [0, "PorousDomain.Parts_Solid_layer_4|1", Kratos.YOUNG_MODULUS]]
        input_values = [young_modulus] * len(input_parameters)
        output_parameters = [[0, np.min, test_helper.get_nodal_variable, "PorousDomain.porous_computational_model_part", Kratos.DISPLACEMENT_Y],
                             [1, None, test_helper.get_nodal_variable, "PorousDomain.porous_computational_model_part", Kratos.DISPLACEMENT_Y, 13]]

        try:
            prob_analysis_instance = prob_analysis(template_project_path, input_parameters, output_parameters)
            output_values = prob_analysis_instance.calculate(input_values)
            min_displacements.append(output_values[0])
            displacements_at_node_13.append(output_values[1])
        except Exception as e:
            print("An error occurred:", e)
        finally:
            # Clean up and remove the temporary folder
            prob_analysis_instance.finalize()

    print("Minimum displacements:", min_displacements)
    print("Displacements at node 13:", displacements_at_node_13)