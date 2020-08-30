# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class CouetteFlowTest(KratosUnittest.TestCase):
    def testCouetteFlow2DSymbolicStokes(self):
        self.solver_type = "MonolithicStokes"
        self.element_type = "symbolic_stokes"
        self.time_scheme = "bdf2"
        self._CustomizeSimulationSettings()

    def testCouetteFlow2DSymbolicNavierStokes(self):
        self.solver_type = "Monolithic"
        self.element_type = "symbolic"
        self.time_scheme = "bdf2"
        self._CustomizeSimulationSettings()

    def setUp(self):
        self.print_output = False
        self.check_tolerance = 1.0e-6
        self.print_reference_values = False
        self.work_folder = "CouetteFlowTest"
        settings_filename = "ProjectParameters.json"

        # Read the simulation settings
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            with open(settings_filename,'r') as parameter_file:
                self.parameters = KratosMultiphysics.Parameters(parameter_file.read())

        # If required, add the output process
        if self.print_output:
            self._AddOutput()

    def runTest(self):
        # Create the test simulation
        with KratosUnittest.WorkFolderScope(self.work_folder,__file__):
            self.model = KratosMultiphysics.Model()
            simulation = FluidDynamicsAnalysis(self.model, self.parameters)
            simulation.Run()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting('couette_flow_test.time')

    def checkResults(self):
        fluid_model_part = self.model.GetModelPart("FluidModelPart")
        self.reference_file = "reference_couette_flow_test_" + self.element_type
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            if self.print_reference_values:
                with open(self.reference_file + '.csv','w') as ref_file:
                    ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, PRESSURE\n")
                    for node in fluid_model_part.Nodes:
                        vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                        temp = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE,0)
                        ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], temp))
            else:
                with open(self.reference_file + '.csv','r') as reference_file:
                    reference_file.readline() # skip header
                    line = reference_file.readline()

                    for node in fluid_model_part.Nodes:
                        values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                        #node_id = values[0]
                        reference_vel_x = values[1]
                        reference_vel_y = values[2]
                        reference_press = values[3]

                        velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                        self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                        self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)
                        pressure = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                        self.assertAlmostEqual(reference_press, pressure, delta = self.check_tolerance)

                        line = reference_file.readline()
                    if line != '': # If we did not reach the end of the reference file
                        self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def _CustomizeSimulationSettings(self):
        # Customize simulation settings
        self.parameters["solver_settings"]["solver_type"].SetString(self.solver_type)
        self.parameters["solver_settings"]["formulation"]["element_type"].SetString(self.element_type)
        self.parameters["solver_settings"]["time_scheme"].SetString(self.time_scheme)

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                "output_name"            : "couette_flow_test",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "time",
                        "output_control_type"         : "step",
                        "output_interval"             : 1,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["VELOCITY","PRESSURE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        output_name = gid_output_settings["Parameters"]["output_name"].GetString()
        output_name += "_" + self.element_type
        gid_output_settings["Parameters"]["output_name"].SetString(output_name)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

if __name__ == '__main__':
    test = CouetteFlowTest()
    test.setUp()
    # test.testCouetteFlow2DSymbolicStokes()
    test.testCouetteFlow2DSymbolicNavierStokes()
    test.runTest()
    test.checkResults()
    test.tearDown()
