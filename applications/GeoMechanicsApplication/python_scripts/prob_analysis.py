import os
import tempfile
import shutil

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as kratos_geo
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
import applications.GeoMechanicsApplication.tests.test_helper as test_helper

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
                print(type(input_parameter[0]), type(stage_no))
                print(stage_no, "<", input_parameter[0], "=", stage_no < input_parameter[0])
                input()
                if stage_no < input_parameter[0]:
                    print("Skipping input parameter for stage", stage_no)
                    continue
                input(f"Stage {stage_no} > {input_parameter[0]} ..... Press Enter to continue")
                value = input_values[i]
                self.update_material_parameters(input_parameter, value)

            stage_obj.RunSolutionLoop()
            stage_obj.Finalize()

            for output_parameter in self.output_parameters:
                if stage_no != output_parameter[0]:
                    continue
                model_part = self.model.GetModelPart(output_parameter[2])
                output_values.append(output_parameter[1](model_part, *output_parameter[3:]))

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
    template_project_path = r"c:\Users\nuttall\OneDrive - Stichting Deltares\Desktop\dummy 2 stage test"
    input_parameters = [[1, "PorousDomain.Soil", kratos_geo.UMAT_PARAMETERS, 0],
                        [1, "PorousDomain.Soil", kratos_geo.UMAT_PARAMETERS, 1]]
    input_values = [5000000.0, 0.2]
    output_parameters = [[0, test_helper.get_nodal_variable, "PorousDomain.Soil", Kratos.DISPLACEMENT_Y, 1],
                         [1, test_helper.get_nodal_variable, "PorousDomain.Soil", Kratos.DISPLACEMENT_Y, 1]]


    try:
        prob_analysis_instance = prob_analysis(template_project_path, input_parameters, output_parameters)
        output_values = prob_analysis_instance.calculate(input_values)
        print(output_values)
    except Exception as e:
        print("An error occurred:", e)
    finally:
        # Clean up and remove the temporary folder
        prob_analysis_instance.finalize()