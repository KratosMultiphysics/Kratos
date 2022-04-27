import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis_rve import FluidDynamicsAnalysisRVE

# Previous line will import the Fluid RVE analysis class

class TestFluidRVETest(KratosUnittest.TestCase):
    
    def setUp(self):
        # Change to True to generate GiD Output 
        self.print_output = False
        # Change to True if no refence value file exists
        self.check_tolerance = 1e-4
        self.print_reference_values = False
    
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
            
            self.reference_file = "reference_rve_2D"
                        
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
            
            self.reference_file = "reference_rve_3D"
            
            self._aux_fluid_rve_computation(parameters)
        
    def _aux_fluid_rve_computation(self, parameters):
        
        self.domain_size = parameters["solver_settings"]["domain_size"].GetInt()
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
            output_settings[0]["Parameters"]["output_name"].SetString("FluidRVETest/fluid_rve_test_{}D".format(self.domain_size))
            parameters["output_processes"].AddValue("gid_output", output_settings)
        
        model = KratosMultiphysics.Model()
        self.simulation = FluidDynamicsAnalysisRVE(model, parameters)
        self.simulation.Run()
               
        self._CheckResults()
        
    def _CheckResults(self):
        model_part = self.simulation._GetSolver().GetComputingModelPart()
        
        if self.print_reference_values:
            with open('FluidRVETest/' + self.reference_file + '.csv','w') as ref_file:
                if self.domain_size == 2:
                    ref_file.write("#ID, VELOCITY_X, VELOCITY_Y\n")
                else:
                    ref_file.write("#Id, VELOCITY_X, VELOCITY_Y, VELOCITY_Z\n")
                for node in model_part.Nodes:
                    vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                    if self.domain_size == 2:
                        ref_file.write("{0}, {1}, {2}\n".format(node.Id, vel[0], vel[1]))
                    else:
                        ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], vel[2]))
        else:
            with open('FluidRVETest/' + self.reference_file + '.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()

                for node in model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    reference_vel_x = values[1]
                    reference_vel_y = values[2]
                    if self.domain_size == 3:
                        reference_vel_z = values[3]

                    velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                    self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                    self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)
                    if self.domain_size == 3:
                        self.assertAlmostEqual(reference_vel_z, velocity[2], delta = self.check_tolerance)
                   

                    line = reference_file.readline()
                if line != '': 
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")
            


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()  
        
            