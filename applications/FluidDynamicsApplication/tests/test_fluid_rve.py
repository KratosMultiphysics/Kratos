import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication # Might be needed for testing
import KratosMultiphysics.KratosUnittest as KratosUnittest

import KratosMultiphysics.kratos_utilities as kratos_utilities # Might be needed for testing

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis_rve import FluidDynamicsAnalysisRVE

# Previous line will import the Fluid RVE analysis class

class TestFluidRVETest(KratosUnittest.TestCase):
    
    def setUp(self):
        self.print_output = False
    
    def test_fluid_rve_computation_2d(self):
        #Within location context:
        with KratosUnittest.WorkFolderScope(".",__file__):            
            with open("FluidRVETest/fluid_rve_test_parameters.json", 'r') as parameter_file:
                parameters =  KratosMultiphysics.Parameters(parameter_file.read())
                
            parameters["solver_settings"]["domain_size"].SetInt(2)
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("FluidRVETest/fluid_rve_test_2D")
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("FluidRVETest/fluid_rve_test_materials_2D.json")
            parameters["rve_settings"]["boundary_mp_name"].SetString("FluidModelPart.Slip2D.Boundaries")
            parameters["solver_settings"]["skin_parts"][0].SetString("Slip2D")
            parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][0].SetInt(0)
            parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][1].SetInt(0)
                        
            self._aux_fluid_rve_computation(parameters)
        
    def test_fluid_rve_computation_3d(self):
        #Within location context:
        with KratosUnittest.WorkFolderScope(".",__file__):
            with open("FluidRVETest/fluid_rve_test_parameters.json", 'r') as parameter_file:
                parameters =  KratosMultiphysics.Parameters(parameter_file.read())
                
            parameters["solver_settings"]["domain_size"].SetInt(3)
            parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("FluidRVETest/fluid_rve_test_3D")
            parameters["solver_settings"]["material_import_settings"]["materials_filename"].SetString("FluidRVETest/fluid_rve_test_materials_3D.json")
            parameters["rve_settings"]["boundary_mp_name"].SetString("FluidModelPart.Slip3D.Boundaries")
            parameters["solver_settings"]["skin_parts"][0].SetString("Slip3D")
            parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][0].SetInt(0)
            parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][1].SetInt(0)
            parameters["processes"]["initial_conditions_process_list"][0]["Parameters"]["value"][2].SetInt(0)
                       
            self._aux_fluid_rve_computation(parameters)
        
    def _aux_fluid_rve_computation(self, parameters):
        
        domain_size = parameters["solver_settings"]["domain_size"].GetInt()
        if self.print_output :
            output_settings = KratosMultiphysics.Parameters(R'''[{
                "python_module" : "gid_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "GiDOutputProcess",
                "help"          : "This process writes postprocessing files for GiD",
                "Parameters"    : {
                    "model_part_name"        : "FluidModelPart",
                    "output_name"            : "TO_BE_SET",
                    "postprocess_parameters" : {
                        "result_file_configuration" : {
                            "gidpost_flags"         : {
                                "GiDPostMode"           : "GiD_PostBinary",
                                "WriteDeformedMeshFlag" : "WriteDeformed",
                                "WriteConditionsFlag"   : "WriteConditions",
                                "MultiFileFlag"         : "SingleFile"
                            },
                            "file_label"          : "time",
                            "output_control_type" : "step",
                            "output_interval"     : 1,
                            "body_output"         : true,
                            "node_output"         : false,
                            "skin_output"         : false,
                            "plane_output"        : [],
                            "nodal_results"       : ["VELOCITY","PRESSURE"],
                            "gauss_point_results" : [],
                            "nodal_flags_results": ["MASTER","SLAVE"]
                        },
                        "point_data_configuration"  : []
                    }

                }
            }]''')
            output_settings[0]["Parameters"]["output_name"].SetString("FluidRVETest/fluid_rve_test_{}D".format(domain_size))
            parameters["output_processes"].AddValue("gid_output", output_settings)
        
        model = KratosMultiphysics.Model()
        simulation = FluidDynamicsAnalysisRVE(model, parameters)
        simulation.Run()
        
        # Space for testing


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()  
        
            