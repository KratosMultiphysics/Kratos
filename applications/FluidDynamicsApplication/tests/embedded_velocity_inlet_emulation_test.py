import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.kratos_utilities as KratosUtilities
have_external_solvers = KratosUtilities.CheckIfApplicationsAvailable("ExternalSolversApplication")

import sys
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class EmbeddedVelocityInletEmulationTest(UnitTest.TestCase):
    def testEmbeddedVelocityInletEmulationSymbolic2D(self):
        self.print_output = False
        self.check_tolerance = 1.0e-10
        self.print_reference_values = False
        self.work_folder = "EmbeddedVelocityInletEmulationTest"
        self.reference_file = "reference_embedded_symbolic_navier_stokes"
        self.settings = "EmbeddedVelocityInletEmulationTest.json"
        self.formulation_settings = KratosMultiphysics.Parameters(r'''{
            "element_type"        : "embedded_symbolic_navier_stokes",
            "is_slip"             : false,
            "slip_length"         : 1.0e8,
            "penalty_coefficient" : 10.0,
            "dynamic_tau"         : 1.0
        }''')
        self.ExecuteTest()

    def testEmbeddedVelocityInletEmulationEmbedded2D(self):
        self.print_output = False
        self.check_tolerance = 1.0e-10
        self.print_reference_values = False
        self.work_folder = "EmbeddedVelocityInletEmulationTest"
        self.reference_file = "reference_embedded_navier_stokes"
        self.settings = "EmbeddedVelocityInletEmulationTest.json"
        self.formulation_settings = KratosMultiphysics.Parameters(r'''{
            "element_type": "embedded_navier_stokes",
            "is_slip": false,
            "slip_length": 1.0e8,
            "penalty_coefficient": 10.0,
            "dynamic_tau": 1.0
        }''')
        self.ExecuteTest()

    def ExecuteTest(self):
        self.buildSimulation()
        self.runTest()
        self.tearDown()
        self.checkResults()

    def buildSimulation(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
                self.ProjectParameters["solver_settings"]["formulation"] = self.formulation_settings
                if self.print_output:
                    self.ProjectParameters["output_processes"] = self._get_output_process_settings()

        self.model = KratosMultiphysics.Model()
        self.simulation = FluidDynamicsAnalysis(self.model,self.ProjectParameters)

    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            self.simulation.Initialize()
            self._set_inlet_emulation_simulation()
            self.simulation.RunSolutionLoop()
            self.simulation.Finalize()

    def checkResults(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            ## 2D results check
            main_model_part = self.simulation._GetSolver().main_model_part
            if (main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                if self.print_reference_values:
                    with open(self.reference_file+'.csv','w') as ref_file:
                        ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, PRESSURE\n")
                        for node in main_model_part.Nodes:
                            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                            pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE,0)
                            ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], pres))
                else:
                    with open(self.reference_file+'.csv','r') as reference_file:
                        reference_file.readline() # skip header
                        line = reference_file.readline()

                        for node in main_model_part.Nodes:
                            values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                            reference_vel_x = values[1]
                            reference_vel_y = values[2]
                            reference_pres = values[3]

                            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                            pressure = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                            self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                            self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)
                            self.assertAlmostEqual(reference_pres, pressure, delta = self.check_tolerance)

                            line = reference_file.readline()
                        if line != '': # If we did not reach the end of the reference file
                            self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            KratosUtilities.DeleteFileIfExisting(
                self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')

    def _set_inlet_emulation_simulation(self):
        main_model_part = self.simulation._GetSolver().main_model_part
        # Set distance field
        x_level_set = 2.5
        for node in main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.X - x_level_set)

        # Set embedded velocity in cut elements and deactivate negative distance elements
        embedded_velocity = KratosMultiphysics.Vector(3)
        embedded_velocity[0] = 1.0
        embedded_velocity[1] = 0.0
        embedded_velocity[2] = 0.0
        deactivate_negative_elems = True
        for elem in main_model_part.Elements:
            n_pos = 0
            n_neg = 0
            for node in elem.GetNodes():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                    n_neg += 1
                else:
                    n_pos += 1
            if (n_neg != 0 and n_pos != 0):
                for node in elem.GetNodes():
                    node.SetValue(KratosMultiphysics.EMBEDDED_VELOCITY, embedded_velocity)
            elif (n_neg == len(elem.GetNodes()) and deactivate_negative_elems):
                elem.Set(KratosMultiphysics.ACTIVE, False)

    def _get_output_process_settings(self):
        return KratosMultiphysics.Parameters(r'''{
            "gid_output" : [{
                "python_module" : "gid_output_process",
                "kratos_module" : "KratosMultiphysics",
                "process_name"  : "GiDOutputProcess",
                "help"          : "This process writes postprocessing files for GiD",
                "Parameters"    : {
                    "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                    "output_name"            : "embedded_velocity_inlet_emulation_test",
                    "postprocess_parameters" : {
                        "result_file_configuration" : {
                            "gidpost_flags"       : {
                                "GiDPostMode"           : "GiD_PostBinary",
                                "WriteDeformedMeshFlag" : "WriteDeformed",
                                "WriteConditionsFlag"   : "WriteConditions",
                                "MultiFileFlag"         : "SingleFile"
                            },
                            "file_label"          : "time",
                            "output_control_type" : "step",
                            "output_interval"     : 1.0,
                            "body_output"         : true,
                            "node_output"         : false,
                            "skin_output"         : false,
                            "plane_output"        : [],
                            "nodal_results"       : ["VELOCITY","PRESSURE","DISTANCE"],
                            "gauss_point_results" : []
                        },
                        "point_data_configuration"  : []
                    }
                }
            }]
        }''')

if __name__ == '__main__':
    test = EmbeddedVelocityInletEmulationTest()
    test.testEmbeddedVelocityInletEmulationSymbolic2D()
    test.testEmbeddedVelocityInletEmulationEmbedded2D()
