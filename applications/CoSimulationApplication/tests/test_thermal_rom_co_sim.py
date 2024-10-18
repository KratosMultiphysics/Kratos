import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    import KratosMultiphysics.ConvectionDiffusionApplication

if kratos_utilities.CheckIfApplicationsAvailable("CoSimulationApplication"):
    from KratosMultiphysics.CoSimulationApplication.co_simulation_analysis import CoSimulationAnalysis

@KratosUnittest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication", "RomApplication", "MappingApplication")
class TestThermalRomCoSim(KratosUnittest.TestCase):

    def testConvDiffStationaryRom2DCoSim(self):
        self.work_folder = "thermal_rom_co_sim_test_files/"
        parameter_file_name = "ProjectParameters_CoSimulation_rom.json"

        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            # Set up simulation
            with open(parameter_file_name, 'r') as parameter_file:
                cosim_parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Generate output.dat for center domain
            parameter_file_name_center = "ProjectParameters_center.json"
            with open(parameter_file_name_center,'r') as f:
                center_parameters = KratosMultiphysics.Parameters(f.read())
                center_parameters["processes"].AddValue("testing_processes", KratosMultiphysics.Parameters("""[{
                        "kratos_module"   : "KratosMultiphysics",
                        "python_module"   : "point_output_process",
                        "process_name"    : "PointOutputProcess",
                        "Parameters" : {
                            "position"         : [0.5, 0.8, 0.0],
                            "entity_type"      : "node",
                            "model_part_name"  : "ThermalModelPart",
                            "output_file_settings": {
                                "file_name"  : "center_output.dat"
                            },
                            "output_variables" : [
                                "TEMPERATURE",
                                "REACTION_FLUX"]
                            }
                        },{
                        "python_module"   : "compare_two_files_check_process",
                        "kratos_module"   : "KratosMultiphysics",
                        "process_name"    : "CompareTwoFilesCheckProcess",
                        "Parameters" :{
                            "output_file_name"    : "center_output.dat",
                            "reference_file_name" : "center_reference.dat",
                            "comparison_type"     : "dat_file_variables_time_history",
                            "tolerance"      : 1e-8,
                            "relative_tolerance"    : 1e-8
                            }
                        }]"""))

            with open(parameter_file_name_center,'w') as f:
                f.write(center_parameters.PrettyPrintJsonString())

            # Generate output.dat for outside domain
            parameter_file_name_outside = "ProjectParameters_outside.json"
            with open(parameter_file_name_outside,'r') as f:
                outside_parameters = KratosMultiphysics.Parameters(f.read())
                outside_parameters["processes"].AddValue("testing_processes", KratosMultiphysics.Parameters("""[{
                        "kratos_module"   : "KratosMultiphysics",
                        "python_module"   : "point_output_process",
                        "process_name"    : "PointOutputProcess",
                        "Parameters" : {
                            "position"         : [0.5, 0.8, 0.0],
                            "entity_type"      : "node",
                            "model_part_name"  : "ThermalModelPart",
                            "output_file_settings": {
                                "file_name"  : "outside_output.dat"
                            },
                            "output_variables" : [
                                "TEMPERATURE",
                                "REACTION_FLUX"]
                            }
                        },{
                        "python_module"   : "compare_two_files_check_process",
                        "kratos_module"   : "KratosMultiphysics",
                        "process_name"    : "CompareTwoFilesCheckProcess",
                        "Parameters" :{
                            "output_file_name"    : "outside_output.dat",
                            "reference_file_name" : "outside_reference.dat",
                            "comparison_type"     : "dat_file_variables_time_history",
                            "tolerance"      : 1e-8,
                            "relative_tolerance"    : 1e-8
                            }
                        }]"""))

            with open(parameter_file_name_outside,'w') as f:
                f.write(outside_parameters.PrettyPrintJsonString())

            simulation = CoSimulationAnalysis(cosim_parameters)
            simulation.Run()

            # After the simulation, clean up the center JSON file
            with open(parameter_file_name_center,'r') as f:
                center_parameters = KratosMultiphysics.Parameters(f.read())
                if center_parameters["processes"].Has("testing_processes"):
                    center_parameters["processes"].RemoveValue("testing_processes")

            with open(parameter_file_name_center,'w') as f:
                f.write(center_parameters.PrettyPrintJsonString())

            # After the simulation, clean up the outside JSON file
            with open(parameter_file_name_outside,'r') as f:
                outside_parameters = KratosMultiphysics.Parameters(f.read())
                if outside_parameters["processes"].Has("testing_processes"):
                    outside_parameters["processes"].RemoveValue("testing_processes")

            with open(parameter_file_name_outside,'w') as f:
                f.write(outside_parameters.PrettyPrintJsonString())

##########################################################################################

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()