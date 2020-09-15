import KratosMultiphysics as km
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis

import os

def RunParametricTestCase(
    settings_file_name,
    work_folder,
    parameters_dict,
    print_output = False):

    test_folder = os.path.join("..", "tests", work_folder)

    with UnitTest.WorkFolderScope(test_folder, __file__):
        model = km.Model()
        with open(settings_file_name, 'r') as settings_file:
            file_data = settings_file.read()

        for key, value in parameters_dict.items():
            file_data = file_data.replace(key, value)

        settings = km.Parameters(file_data)

        # to check the results: add output settings block if needed
        if print_output:
            settings.AddValue(
                "output_processes",
                km.Parameters(r'''{
                "gid_output" : [{
                    "python_module" : "gid_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "GiDOutputProcess",
                    "help"          : "This process writes postprocessing files for GiD",
                    "Parameters"    : {
                        "model_part_name"        : "fluid_computational_model_part",
                        "output_name"            : "interface_test",
                        "postprocess_parameters" : {
                            "result_file_configuration" : {
                                "gidpost_flags" : {
                                    "GiDPostMode"           : "GiD_PostBinary",
                                    "WriteDeformedMeshFlag" : "WriteUndeformed",
                                    "WriteConditionsFlag"   : "WriteElementsOnly",
                                    "MultiFileFlag"         : "SingleFile"
                                },
                                "file_label"          : "time",
                                "output_control_type" : "step",
                                "output_frequency"    : 1,
                                "body_output"         : true,
                                "node_output"         : false,
                                "skin_output"         : false,
                                "plane_output"        : [],
                                "nodal_results"       : ["VELOCITY","PRESSURE"],
                                "gauss_point_results" : []
                            },
                            "point_data_configuration"  : []
                        }
                    }
                }]
            }'''))

        analysis = RANSAnalysis(model, settings)
        analysis.Run()