from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities

have_convection_diffusion = KratosUtilities.CheckIfApplicationsAvailable("ConvectionDiffusionApplication")
if have_convection_diffusion:
    from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

    class BuoyancyTestConvectionDiffusionAnalysis(ConvectionDiffusionAnalysis):
        def ModifyAfterSolverInitialize(self):
            super().ModifyAfterSolverInitialize()

            xmin = 0.0
            xmax = 1.0
            ymin = 0.0
            ymax = 1.0

            ## For Ra~1e6
            g = 9.81     # acceleration of gravity m/s2
            T1 = 293.15  # Cold (reference) temperature K
            T2 = 303.15  # Hot temperature K
            rho = 1.2039 # (reference) density kg/m3
            c = 1004.84  # Specific heat J/kg K
            k = ( (rho*c)**2*g*(T2-T1)*(xmax-xmin)**3 / (1e6*T1*0.71) )**0.5 # Given Ra=1e6 & Pr=0.71
            mu = 0.71*k/c # For Prandlt = 0.71
            nu = mu/rho

            ## Set initial and boundary conditions
            model_part = self._GetSolver().GetComputingModelPart()
            model_part.ProcessInfo.SetValue(AMBIENT_TEMPERATURE,T1)
            for node in model_part.Nodes:
                node.SetSolutionStepValue(DENSITY, rho)
                node.SetSolutionStepValue(VISCOSITY, nu)
                node.SetSolutionStepValue(CONDUCTIVITY, k)
                node.SetSolutionStepValue(SPECIFIC_HEAT, c)

                if node.X == xmin or node.X == xmax or node.Y == ymin or node.Y == ymax:
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    if node.X == xmin and node.Y == ymin:
                        node.Fix(PRESSURE)

                if node.X == xmin:
                    node.Fix(TEMPERATURE)
                    node.SetSolutionStepValue(TEMPERATURE, T2)
                elif node.X == xmax:
                    node.Fix(TEMPERATURE)
                    node.SetSolutionStepValue(TEMPERATURE, T1)
                else:
                    T = T2 + (T1-T2)*(node.X-xmin)/(xmax-xmin)
                    node.SetSolutionStepValue(TEMPERATURE, T)

@UnitTest.skipIfApplicationsNotAvailable("ConvectionDiffusionApplication")
class BuoyancyTest(UnitTest.TestCase):

    def testEulerian(self):
        self.convection_diffusion_solver = "eulerian"
        self.reference_file = "reference10_eulerian"
        self.thermal_expansion_coefficient = None # If set, it will be used instead of 1./AmbientTemperature

    def testThermalExpansionCoefficient(self):
        self.convection_diffusion_solver = "eulerian"
        self.reference_file = "reference10_eulerian"
        self.thermal_expansion_coefficient = 1.0/293.15

    def testBFECC(self):
        self.convection_diffusion_solver = "bfecc"
        self.reference_file = "reference10_bfecc"
        self.thermal_expansion_coefficient = None # If set, it will be used instead of 1./AmbientTemperature
        self.check_tolerance = 1e-1 # The bfecc solver shows some variation between runs, we cannot be too strict here

    def validationEulerian(self):
        self.convection_diffusion_solver = "eulerian"
        self.reference_file = "reference80_eulerian"
        self.thermal_expansion_coefficient = None # If set, it will be used instead of 1./AmbientTemperature

        #TODO: MODIFY THIS IN THE JSON SETTINGS
        self.nsteps = 200
        self.input_file = "cavity80"

    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

        # Open parameters
        with UnitTest.WorkFolderScope("BuoyancyTest", __file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                self.parameters = Parameters(parameter_file.read())

    def runTest(self):
        with UnitTest.WorkFolderScope("BuoyancyTest", __file__):
            # If required, add the output process to the test settings
            if self.print_output:
                self._AddOutput()

            # Set thermal solver
            if self.convection_diffusion_solver == 'bfecc':
                self.parameters["solver_settings"]["thermal_solver_settings"]["time_integration_method"].SetString("semi_eulerian")
            elif self.convection_diffusion_solver == 'eulerian':
                self.parameters["solver_settings"]["thermal_solver_settings"]["time_integration_method"].SetString("implicit")
            else:
                raise Exception("Unsupported convection-diffusion solver option: {0}".format(self.convection_diffusion_solver))

            # Set buoyancy process
            gravity = Vector(3)
            gravity[0] = 0.0
            gravity[1] = 9.81
            gravity[2] = 0.0
            self.parameters["processes"]["boussinesq_process_list"][0]["Parameters"].AddVector("gravity", gravity)
            if self.thermal_expansion_coefficient is not None:
                self.parameters["processes"]["boussinesq_process_list"][0]["Parameters"].AddDouble("thermal_expansion_coefficient", self.thermal_expansion_coefficient)

            # Set up the modified analysis stage with initial and boundary conditions and run
            model = Model()
            self.simulation = BuoyancyTestConvectionDiffusionAnalysis(model, self.parameters)
            self.simulation.Run()

            # Check results
            self.checkResults()

    def checkResults(self):
        model_part = self.simulation._GetSolver().GetComputingModelPart()
        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, TEMPERATURE\n")
                for node in model_part.Nodes:
                    vel = node.GetSolutionStepValue(VELOCITY,0)
                    temp = node.GetSolutionStepValue(TEMPERATURE,0)
                    ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], temp))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()
                for node in model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_vel_x = values[1]
                    reference_vel_y = values[2]
                    reference_temp = values[3]

                    velocity = node.GetSolutionStepValue(VELOCITY)
                    self.assertAlmostEqual(reference_vel_x, velocity[0], delta=self.check_tolerance)
                    self.assertAlmostEqual(reference_vel_y, velocity[1], delta=self.check_tolerance)
                    temperature = node.GetSolutionStepValue(TEMPERATURE)
                    self.assertAlmostEqual(reference_temp, temperature, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def tearDown(self):
        with UnitTest.WorkFolderScope("BuoyancyTest", __file__):
            filename = self.parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]["input_filename"].GetString()
            KratosUtilities.DeleteFileIfExisting(filename + '.time')

    def _AddOutput(self):
        gid_output_settings = Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart",
                "output_name"            : "TO_BE_DEFINED",
                "postprocess_parameters" : {
                    "result_file_configuration" : {
                        "gidpost_flags"               : {
                            "GiDPostMode"           : "GiD_PostBinary",
                            "WriteDeformedMeshFlag" : "WriteDeformed",
                            "WriteConditionsFlag"   : "WriteConditions",
                            "MultiFileFlag"         : "SingleFile"
                        },
                        "file_label"                  : "step",
                        "output_control_type"         : "step",
                        "output_frequency"            : 1.0,
                        "body_output"                 : true,
                        "node_output"                 : false,
                        "skin_output"                 : false,
                        "plane_output"                : [],
                        "nodal_results"               : ["VELOCITY","PRESSURE","TEMPERATURE","DENSITY","VISCOSITY","CONDUCTIVITY","SPECIFIC_HEAT","BODY_FORCE"]
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        gid_output_settings["Parameters"]["output_name"].SetString(self.input_file)
        self.parameters["output_processes"]["gid_output"].Append(gid_output_settings)

if __name__ == '__main__':
    test = BuoyancyTest()
    test.setUp()
    #test.print_reference_values = True
    test.print_output = True
    test.testEulerian()
    #test.testThermalExpansionCoefficient()
    #test.testBFECC()
    test.runTest()
    test.tearDown()
