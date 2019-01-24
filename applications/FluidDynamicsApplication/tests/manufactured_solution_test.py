from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False

import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Python modules
import math
import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@KratosUnittest.skipUnless(have_external_solvers, "Missing required application: ExternalSolversApplication")
class ManufacturedSolutionTest(KratosUnittest.TestCase):
    def testManufacturedSolution(self):
        self.runTest()

    def setUp(self):
        self.print_output = False
        self.print_convergence_plot = False
        self.problem_type = "manufactured_solution" # Available problem types: "manufactured_solution" "analytical_solution"
        self.analytical_solution_type = "sinusoidal_transient_field" # Available fields: "nonlinear_transient_field" "sinusoidal_transient_field" "nonlinear_stationary_field"

        self.work_folder = "ManufacturedSolutionTest"
        self.settings = "ManufacturedSolutionTestParameters.json"

        self.meshes_list = ["manufactured_solution_ref0",
                            "manufactured_solution_ref1",
                            "manufactured_solution_ref2",
                            "manufactured_solution_ref3"]
                            #"manufactured_solution_ref4"]

    def tearDown(self):
        with WorkFolderScope(self.work_folder):
            for filename in self.meshes_list:
                try:
                    os.remove(filename + '.time')
                except FileNotFoundError as e:
                    pass

    def runTest(self):
        with WorkFolderScope(self.work_folder):
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
                FluidProblem.SetFluidProblem()
                FluidProblem.SolveFluidProblem()

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
            expected_velocity_errors = [0.01708951546622635, 0.005366727106714455, 0.0013142808355902074, 0.00032206907919625683, 8.037719698951708e-05]
            expected_pressure_errors = [44.03061907965929, 4.8775536490608316, 0.8950814197625788, 0.2200468445178847, 0.0666813658821848]

            for i in range(len(self.meshes_list)):
                self.assertAlmostEqual(err_v[i], expected_velocity_errors[i])
                self.assertAlmostEqual(err_p[i], expected_pressure_errors[i])

class ManufacturedSolutionProblem:

    def __init__(self, ProjectParameters, input_file_name, print_output, problem_type, analytical_solution_type):

        self.problem_type = problem_type
        self.print_output = print_output
        self.input_file_name = input_file_name
        self.ProjectParameters = ProjectParameters
        self.analytical_solution_type = analytical_solution_type
        self.model = KratosMultiphysics.Model()


    def SetFluidProblem(self):

        ## Set the current mesh case problem info
        if (self.problem_type == "analytical_solution"):
            self.ProjectParameters["problem_data"]["problem_name"].SetString(self.input_file_name+"_manufactured")
        else:
            self.ProjectParameters["problem_data"]["problem_name"].SetString(self.input_file_name)
        self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.input_file_name)

        ## Solver construction
        import python_solvers_wrapper_fluid
        self.solver = python_solvers_wrapper_fluid.CreateSolver(self.model, self.ProjectParameters)

        self.solver.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        self.solver.ImportModelPart()
        self.solver.PrepareModelPart()

        self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

        ## Add AddDofs
        self.solver.AddDofs()

        ## Initialize GiD  I/O
        if (self.print_output):
            from gid_output_process import GiDOutputProcess
            self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
                                               self.ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                               self.ProjectParameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

        ## Solver initialization
        self.solver.Initialize()

        ## Compute and set the nodal area
        self.SetNodalArea()

        ## Set the distance to 1 to have full fluid elements
        if (self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "Embedded"):
            for node in self.main_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)

        ## Fix the pressure in one node (bottom left corner)
        for node in self.main_model_part.Nodes:
            if ((node.X<0.001) and (node.Y<0.001)):
                node.Fix(KratosMultiphysics.PRESSURE)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, 0.0)

    def SolveFluidProblem(self):

        ## Stepping and time settings
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = 0.0

        if (self.print_output):
            self.gid_output.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            time = self.solver.AdvanceInTime(time)

            if (self.print_output):
                self.gid_output.ExecuteInitializeSolutionStep()

            if (self.problem_type == "analytical_solution"):
                # Fix the manufactured solution values (only for visualization purposes)
                self.SetManufacturedSolutionValues(fix=True, set_only_boundaries=False)
            else:
                # Set the manufactured solution source terms
                self.SetManufacturedSolutionValues(fix=True, set_only_boundaries=True)
                self.SetManufacturedSolutionSourceValues()

            if (self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] < 3):
                self.SetManufacturedSolutionValues(False) # Set the analytical solution in the two first steps
            else:
                if (self.problem_type != "analytical_solution"):
                    self.solver.InitializeSolutionStep()
                    self.solver.Predict()
                    self.solver.SolveSolutionStep()
                    self.solver.FinalizeSolutionStep()

            if (self.print_output):
                self.gid_output.ExecuteFinalizeSolutionStep()

                if self.gid_output.IsOutputStep():
                    self.gid_output.PrintOutput()

        if (self.print_output):
            self.gid_output.ExecuteFinalize()

    def SetManufacturedSolutionValues(self, fix=True, set_only_boundaries=False):
        ## Set the analytical solution for the manufactured solution computation
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        if (set_only_boundaries == False):
            for node in self.main_model_part.Nodes:
                vel = self.ComputeNodalVelocityManufacturedSolution(node, time)
                pres = self.ComputeNodalPressureManufacturedSolution(node)

                if (fix == True):
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)
                    node.Fix(KratosMultiphysics.PRESSURE)

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0, vel[0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0, vel[1])
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, pres)
        else:
            for node in self.main_model_part.GetSubModelPart("Inlet2D_Contour").Nodes:
                vel = self.ComputeNodalVelocityManufacturedSolution(node, time)

                if (fix == True):
                    node.Fix(KratosMultiphysics.VELOCITY_X)
                    node.Fix(KratosMultiphysics.VELOCITY_Y)

                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0, vel[0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0, vel[1])

    def SetNodalArea(self):
        # Compute nodal area
        for element in self.main_model_part.Elements:
            x = []
            y = []
            for node in element.GetNodes():
                x.append(node.X)
                y.append(node.Y)

            Area = 0.5*((x[1]*y[2]-x[2]*y[1])+(x[2]*y[0]-x[0]*y[2])+(x[0]*y[1]-x[1]*y[0])) # Element area (Jacobian/2)
            # print("Element "+str(element.Id)+" area: "+str(Area))

            for node in element.GetNodes():
                aux = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)  # Current nodal area (from other elements)
                aux += Area/3.0                              # Accumulate the current element nodal area
                node.SetSolutionStepValue(KratosMultiphysics.NODAL_AREA, 0, aux)
                node.SetValue(KratosMultiphysics.NODAL_AREA, aux)

        ## Check nodal area computation (squared shaped domain of 1x1 m)
        AreaTotal = 0.0
        for node in self.main_model_part.Nodes:
            # print("Node id "+str(node.Id)+" nodal area: "+str(node.GetValue(KratosMultiphysics.NODAL_AREA)))
            AreaTotal += node.GetValue(KratosMultiphysics.NODAL_AREA)

        if (abs(1.0-AreaTotal) > 1e-5):
            print("Obtained total area: "+str(AreaTotal))
            raise Exception("Error in NODAL_AREA computation.")

    def SetManufacturedSolutionSourceValues(self):
        ## Set the body force as source term
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

        for node in self.main_model_part.Nodes:
            rho = node.GetSolutionStepValue(KratosMultiphysics.DENSITY)
            # If VMS2D element is used, set mu as the Kinematic viscosity
            if (self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "Embedded"):
                mu = node.GetSolutionStepValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            elif (self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "Monolithic"):
                mu = rho*node.GetSolutionStepValue(KratosMultiphysics.VISCOSITY)

            rhof = self.ComputeNodalSourceTermManufacturedSolution(node, time, rho, mu)

            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_X, 0, rhof[0]/rho)  # Set the x-component body force field
            node.SetSolutionStepValue(KratosMultiphysics.BODY_FORCE_Y, 0, rhof[1]/rho)  # Set the y-component body force field

    def ComputeVelocityErrorNorm(self):
        err_v = 0

        for node in self.main_model_part.Nodes:
            weight = node.GetValue(KratosMultiphysics.NODAL_AREA)
            vel_x = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
            vel_y = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
            end_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

            analytical_vel = self.ComputeNodalVelocityManufacturedSolution(node, end_time)

            err_x = analytical_vel[0] - vel_x
            err_y = analytical_vel[1] - vel_y
            err_node = err_x**2 + err_y**2
            err_v += weight*err_node

        return math.sqrt(err_v) # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    def ComputePressureErrorNorm(self):
        err_p = 0

        for node in self.main_model_part.Nodes:
            weight = node.GetValue(KratosMultiphysics.NODAL_AREA)
            pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            analytical_pres = self.ComputeNodalPressureManufacturedSolution(node)
            err_p += weight*(analytical_pres - pres)**2

        return math.sqrt(err_p) # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    def ComputeNodalSourceTermManufacturedSolution(self, node, time, rho, mu):
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

    def ComputeNodalVelocityManufacturedSolution(self, node, time):
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

    def ComputeNodalPressureManufacturedSolution(self, node):
        # We consider solenoidal velocity fields in order to have a known zero pressure solution on the continuum
        return 0.0

if __name__ == '__main__':
    test = ManufacturedSolutionTest()
    test.setUp()
    test.runTest()
    test.tearDown()
