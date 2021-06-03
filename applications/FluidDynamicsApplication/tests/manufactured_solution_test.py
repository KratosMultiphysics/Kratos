from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import Python modules
import math

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

have_external_solvers = KratosUtilities.CheckIfApplicationsAvailable("LinearSolversApplication")

@KratosUnittest.skipUnless(have_external_solvers, "Missing required application: LinearSolversApplication")
class ManufacturedSolutionTest(KratosUnittest.TestCase):
    def testManufacturedSolution(self):
        self.runTest()

    def setUp(self):
        self.print_output = False
        self.print_convergence_plot = False
        self.problem_type = "manufactured_solution" # Available problem types: "manufactured_solution" "analytical_solution"
        self.analytical_solution_type = "sinusoidal_transient_field" # Available fields: "nonlinear_transient_field" "sinusoidal_transient_field" "nonlinear_stationary_field"

        self.work_folder = "manufactured_solution_test"
        self.settings = "ManufacturedSolutionTestParameters.json"

        self.meshes_list = ["manufactured_solution_ref0",
                            "manufactured_solution_ref1",
                            "manufactured_solution_ref2",
                            "manufactured_solution_ref3"]
                            # "manufactured_solution_ref4"]

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            for filename in self.meshes_list:
                KratosUtilities.DeleteFileIfExisting(filename + '.time')

    def runTest(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            with open(self.settings, 'r') as parameter_file:
                self.OriginalProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            h = []
            err_p = []
            err_v = []

            den = 1
            mesh_0_characteristic_size = 0.2

            # Solve the manufactured solution problem for each one of the refinements
            for mesh_name in self.meshes_list:
                # Solve the problem imposing the previously obtained values
                CaseProjectParameters = self.OriginalProjectParameters.Clone()
                FluidProblem = ManufacturedSolutionProblem(CaseProjectParameters, mesh_name, self.print_output, self.problem_type, self.analytical_solution_type)
                FluidProblem.Run()

                # Compute the obtained solution error
                h.append(mesh_0_characteristic_size/den)
                err_p.append(FluidProblem.ComputePressureErrorNorm())
                err_v.append(FluidProblem.ComputeVelocityErrorNorm())
                den *= 2

            # Compute average convergence slopes
            average_slope_pressure = (math.log(err_p[0])-math.log(err_p[-1]))/(math.log(h[0])-math.log(h[-1]))
            average_slope_velocity = (math.log(err_v[0])-math.log(err_v[-1]))/(math.log(h[0])-math.log(h[-1]))

            # Convergence plot print
            if (self.print_convergence_plot == True):
                # Plot the convergence graphs
                import matplotlib.pyplot as plt

                h_1_v = []
                h_2_v = []
                h_3_v = []
                h_1_p = []
                h_2_p = []
                h_3_p = []

                den = 1.0
                for i in range(0,len(err_p)):
                    h_1_v.append(err_v[0]/den)
                    h_2_v.append(err_v[0]/den**2)
                    h_3_v.append(err_v[0]/den**3)
                    h_1_p.append(err_p[0]/den)
                    h_2_p.append(err_p[0]/den**2)
                    h_3_p.append(err_p[0]/den**3)
                    den *= 2

                plt.rc('text', usetex=True)
                plt.rc('font', family='serif')
                plt.loglog(h, err_p, '-+', color = 'r', label = 'Obtained pressure')
                plt.loglog(h, h_1_p, '--', color = 'r', label = 'Linear pressure')
                plt.loglog(h, h_2_p, ':' , color = 'r', label = 'Quadratic pressure')
                plt.loglog(h, h_3_p, '-.', color = 'r', label = 'Cubic pressure')
                plt.loglog(h, err_v, '-x', color = 'k', label = 'Obtained velocity')
                plt.loglog(h, h_1_v, '--', color = 'k', label = 'Linear velocity')
                plt.loglog(h, h_2_v, ':' , color = 'k', label = 'Quadratic velocity')
                plt.loglog(h, h_3_v, ':' , color = 'k', label = 'Cubic velocity')

                plt.title('L2 norm convergence')
                plt.ylabel(r'$\displaystyle\sum_{i=1}^{n_{n}}\frac{A_{i}\Vert\mathbf{u}_{i}-\mathbf{\bar{u}}_{i}\Vert}{A_{T}}$')
                plt.xlabel('h')
                plt.xlim([0.01, 0.25])
                # plt.ylim([1e-4, 1])
                plt.legend(loc=4, ncol=2)
                plt.tight_layout()
                plt.savefig('l2_norm_convergence.png')

            # Check obtained solution
            expected_velocity_errors = [0.01620828837402704, 0.005181479513472467, 0.0012944061647766365, 0.0003215070966125507, 8.577306473469384e-05]
            expected_pressure_errors = [43.89258161784579, 5.237699459045566, 0.9428016719803961, 0.24671881227746087, 0.09607645760134441]

            for i in range(len(self.meshes_list)):
                self.assertAlmostEqual(err_v[i], expected_velocity_errors[i])
                self.assertAlmostEqual(err_p[i], expected_pressure_errors[i])

class ManufacturedSolutionProblem(FluidDynamicsAnalysis):

    def __init__(self, project_parameters, input_file_name, print_output, problem_type, analytical_solution_type):
        self.problem_type = problem_type
        self.print_output = print_output
        self.input_file_name = input_file_name
        self.project_parameters = project_parameters
        self.analytical_solution_type = analytical_solution_type
        self.model = KratosMultiphysics.Model()

        ## Set the current mesh case problem info
        if (self.problem_type == "analytical_solution"):
            self.project_parameters["problem_data"]["problem_name"].SetString(self.input_file_name+"_manufactured")
        else:
            self.project_parameters["problem_data"]["problem_name"].SetString(self.input_file_name)
        self.project_parameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.input_file_name)

        ## If required, set up the GiD I/O
        if self.print_output:
            self._AddOutput()

        ## Note that the base fluid analysis constructor is called after the creation of the model and parameters customization
        super().__init__(self.model, self.project_parameters)

    def ModifyAfterSolverInitialize(self):
        super().ModifyAfterSolverInitialize()

        ## Calculate NODAL_AREA
        nodal_area_process = KratosMultiphysics.CalculateNonHistoricalNodalAreaProcess(self._GetSolver().GetComputingModelPart())
        nodal_area_process.Execute()

        ## Fix the pressure in one node (bottom left corner)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            if ((node.X<0.001) and (node.Y<0.001)):
                node.Fix(KratosMultiphysics.PRESSURE)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, 0.0)

        ## Initialize the buffer with the analytical solution
        for i_buff in range(self._GetSolver().GetComputingModelPart().GetBufferSize()):
            self._SetManufacturedSolutionValues(buffer_position = i_buff, fix = False)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        ## Apply manufactured solution BCs
        if (self.problem_type == "analytical_solution"):
            # Fix the manufactured solution values (only for visualization purposes)
            self._SetManufacturedSolutionValues(fix=True, set_only_boundaries=False)
        else:
            # Set the manufactured solution source terms
            self._SetManufacturedSolutionValues(fix=True, set_only_boundaries=True)
            self._SetManufacturedSolutionSourceValues()

    ## We enhance the class with these two methods to calculate the error norms
    def ComputeVelocityErrorNorm(self):
        err_v = 0

        for node in self._GetSolver().GetComputingModelPart().Nodes:
            weight = node.GetValue(KratosMultiphysics.NODAL_AREA)
            vel_x = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
            vel_y = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            end_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

            analytical_vel = self._ComputeNodalVelocityManufacturedSolution(node, end_time)

            err_x = analytical_vel[0] - vel_x
            err_y = analytical_vel[1] - vel_y
            err_node = err_x**2 + err_y**2
            err_v += weight*err_node

        return math.sqrt(err_v) # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    def ComputePressureErrorNorm(self):
        err_p = 0

        for node in self._GetSolver().GetComputingModelPart().Nodes:
            weight = node.GetValue(KratosMultiphysics.NODAL_AREA)
            pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            analytical_pres = self._ComputeNodalPressureManufacturedSolution(node)
            err_p += weight*(analytical_pres - pres)**2

        return math.sqrt(err_p) # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    ## Internal methods required for the manufactured solution calculation
    def _SetManufacturedSolutionValues(self, buffer_position = 0, fix=True, set_only_boundaries=False):
        ## Set the analytical solution for the manufactured solution computation
        time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

        if set_only_boundaries == False:
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                vel = self._ComputeNodalVelocityManufacturedSolution(node, time)
                pres = self._ComputeNodalPressureManufacturedSolution(node)

                if fix:
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.Fix(KratosMultiphysics.PRESSURE)

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, buffer_position, vel[0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, buffer_position, vel[1])
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, buffer_position, pres)
        else:
            root_model_part = self._GetSolver().GetComputingModelPart().GetRootModelPart()
            for node in root_model_part.GetSubModelPart("Inlet2D_Contour").Nodes:
                vel = self._ComputeNodalVelocityManufacturedSolution(node, time)

                if fix:
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, buffer_position, vel[0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, buffer_position, vel[1])

    def _SetManufacturedSolutionSourceValues(self):
        ## Set the body force as source term
        ## Note 1: it is assumed that the DENSITY is nodally stored
        ## Note 2: it is assumed that the DYNAMIC_VISCOSITY is constant and stored in the properties
        model_part = self._GetSolver().GetComputingModelPart()
        time = model_part.ProcessInfo[KratosMultiphysics.TIME]

        mu = model_part.GetElement(1).Properties[KratosMultiphysics.DYNAMIC_VISCOSITY]
        for node in model_part.Nodes:
            rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
            rho_f = self._ComputeNodalSourceTermManufacturedSolution(node, time, rho, mu)
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, 0, rho_f[0]/rho)  # Set the x-component body force field
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, 0, rho_f[1]/rho)  # Set the y-component body force field

    def _ComputeNodalSourceTermManufacturedSolution(self, node, time, rho, mu):
        if (self.analytical_solution_type == "sinusoidal_transient_field"):
            rhofx = -rho*math.pi*math.sin(math.pi*node.X)*math.cos(math.pi*node.Y)*math.sin(math.pi*time) + 2*mu*math.pi*math.pi*math.sin(math.pi*node.X)*math.cos(math.pi*node.Y)*math.cos(math.pi*time) + rho*math.pi*(math.cos(math.pi*time)**2)*math.sin(math.pi*node.X)*math.cos(math.pi*node.X)
            rhofy =  rho*math.pi*math.cos(math.pi*node.X)*math.sin(math.pi*node.Y)*math.sin(math.pi*time) - 2*mu*math.pi*math.pi*math.cos(math.pi*node.X)*math.sin(math.pi*node.Y)*math.cos(math.pi*time) + rho*math.pi*(math.cos(math.pi*time)**2)*math.sin(math.pi*node.Y)*math.cos(math.pi*node.Y)

        elif (self.analytical_solution_type == "nonlinear_transient_field"):
            rhofx =   rho*math.pi*(node.X**2)*node.Y*math.cos(math.pi*time) - 2*mu*node.Y*math.sin(math.pi*time) + rho*(node.X**3)*(node.Y**2)*((math.sin(math.pi*time))**2)
            rhofy =  -rho*math.pi*node.X*(node.Y**2)*math.cos(math.pi*time) + 2*mu*node.X*math.sin(math.pi*time) + rho*(node.X**2)*(node.Y**3)*((math.sin(math.pi*time))**2)

        elif (self.analytical_solution_type == "nonlinear_stationary_field"):
            rhofx = -2*mu*node.Y + rho*node.X**3*node.Y**2
            rhofy =  2*mu*node.X + rho*node.X**2*node.Y**3

        return [rhofx, rhofy]

    def _ComputeNodalVelocityManufacturedSolution(self, node, time):
        if (self.analytical_solution_type == "sinusoidal_transient_field"):
            vx =  math.sin(math.pi*node.X)*math.cos(math.pi*node.Y)*math.cos(math.pi*time)
            vy = -math.cos(math.pi*node.X)*math.sin(math.pi*node.Y)*math.cos(math.pi*time)

        elif (self.analytical_solution_type == "nonlinear_transient_field"):
            vx = node.X**2*node.Y*math.sin(math.pi*time)
            vy = -node.X*node.Y**2*math.sin(math.pi*time)

        elif (self.analytical_solution_type == "nonlinear_stationary_field"):
            vx = node.X**2*node.Y
            vy = -node.X*node.Y**2

        return [vx, vy]

    def _ComputeNodalPressureManufacturedSolution(self, node):
        # We consider solenoidal velocity fields in order to have a known zero pressure solution on the continuum
        return 0.0

    def _AddOutput(self):
        gid_output_settings = KratosMultiphysics.Parameters("""{
            "python_module" : "gid_output_process",
            "kratos_module" : "KratosMultiphysics",
            "process_name"  : "GiDOutputProcess",
            "help"          : "This process writes postprocessing files for GiD",
            "Parameters"    : {
                "model_part_name"        : "FluidModelPart.fluid_computational_model_part",
                "output_name"            : "TO_BE_SET_BELOW",
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
                        "nodal_results"               : ["VELOCITY","PRESSURE","BODY_FORCE"],
                        "gauss_point_results"         : [],
                        "nodal_nonhistorical_results" : []
                    },
                    "point_data_configuration"  : []
                }
            }
        }""")
        gid_output_settings["Parameters"]["output_name"].SetString(self.input_file_name)
        self.project_parameters["output_processes"]["gid_output"].Append(gid_output_settings)

if __name__ == '__main__':
    KratosUnittest.main()
