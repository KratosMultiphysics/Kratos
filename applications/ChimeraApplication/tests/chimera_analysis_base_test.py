import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError:
    have_external_solvers = False

from fluid_chimera_analysis import FluidChimeraAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class ChimeraAnalysisBaseTest(UnitTest.TestCase):

    def setUp(self):
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        self.print_output = False

    def _run_test(self,settings_file_name):
        model = km.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = km.Parameters(settings_file.read())
        # to check the results: add output settings block if needed
        if self.print_output:
            settings.AddValue("output_processes", km.Parameters(r'''{
                "vtk_output" : [{
                    "python_module" : "vtk_output_process",
                    "kratos_module" : "KratosMultiphysics",
                    "process_name"  : "VtkOutputProcess",
                    "help"          : "This process writes postprocessing files for Paraview",
                    "Parameters"    : {
                        "model_part_name"                    : "FluidModelPart",
                        "output_control_type"                : "step",
                        "output_frequency"                   : 1,
                        "file_format"                        : "ascii",
                        "output_precision"                   : 7,
                        "output_sub_model_parts"             : false,
                        "folder_name"                        : "vtk_output",
                        "save_output_files_in_folder"        : true,
                        "nodal_solution_step_data_variables" : ["VELOCITY","PRESSURE"],
                        "nodal_data_value_variables"         : ["DISTANCE"],
                        "element_flags"                      : ["ACTIVE"],
                        "element_data_value_variables"       : [],
                        "condition_data_value_variables"     : []
                    }
                }]
            }'''))

        analysis = FluidChimeraAnalysis(model,settings)
        analysis.Run()