import KratosMultiphysics as km

from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

class MonolithicKOmegaSSTTest(UnitTest.TestCase):

    def setUp(self):
        # Set to true to get post-process files for the test
        self.print_output = False

    @UnitTest.skipUnless(not km.IsDistributedRun(), "Running with MPI")
    def testMonolithicKOmegaSSTAfcTke(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "OpenMP", "algebraic_flux_corrected", "turbulent_kinetic_energy_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(not km.IsDistributedRun(), "Running with MPI")
    def testMonolithicKOmegaSSTAfcVelocity(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "OpenMP", "algebraic_flux_corrected", "velocity_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(not km.IsDistributedRun(), "Running with MPI")
    def testMonolithicKOmegaSSTRfcTke(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "OpenMP", "residual_based_flux_corrected", "turbulent_kinetic_energy_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(not km.IsDistributedRun(), "Running with MPI")
    def testMonolithicKOmegaSSTRfcVelocity(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "OpenMP", "residual_based_flux_corrected", "velocity_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(km.IsDistributedRun(), "Running without MPI")
    def testMonolithicKOmegaSSTAfcTkeMPI(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "MPI", "algebraic_flux_corrected", "turbulent_kinetic_energy_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(km.IsDistributedRun(), "Running without MPI")
    def testMonolithicKOmegaSSTAfcVelocityMPI(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "MPI", "algebraic_flux_corrected", "velocity_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(km.IsDistributedRun(), "Running without MPI")
    def testMonolithicKOmegaSSTRfcTkeMPI(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "MPI", "residual_based_flux_corrected", "turbulent_kinetic_energy_based")
            kratos_utilities.DeleteTimeFiles(".")

    @UnitTest.skipUnless(km.IsDistributedRun(), "Running without MPI")
    def testMonolithicKOmegaSSTRfcVelocityMPI(self):
        work_folder = "BackwardFacingStepTest"
        settings_file_name = "backward_facing_step_monolithic_k_omega_sst_parameters.json"

        with UnitTest.WorkFolderScope(work_folder, __file__):
            self._runTest(settings_file_name, "MPI", "residual_based_flux_corrected", "velocity_based")
            kratos_utilities.DeleteTimeFiles(".")

    def _runTest(self,settings_file_name, parallel_type, stabilization_method, wall_friction_calculation_method):
        model = km.Model()
        with open(settings_file_name,'r') as settings_file:
            file_data = settings_file.read()

        file_data = file_data.replace("<PARALLEL_TYPE>", parallel_type)
        file_data = file_data.replace("<STABILIZATION_METHOD>", stabilization_method)
        file_data = file_data.replace("<WALL_FRICTION_VELOCITY_CALCULATION_METHOD>", wall_friction_calculation_method)
        settings = km.Parameters(file_data)

        # to check the results: add output settings block if needed
        if self.print_output:
            settings.AddValue("output_processes", km.Parameters(r'''{
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

        analysis = RANSAnalysis(model,settings)
        analysis.Run()

if __name__ == '__main__':
    UnitTest.main()

