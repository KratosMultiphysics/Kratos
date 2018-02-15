from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import Python modules
import math


class KratosExecuteManufacturedSolutionTest(KratosUnittest.TestCase):

    def __init__(self, ProjectParameters):

        # self.meshes_list = ["manufactured_solution_ref0",
        #                     "manufactured_solution_ref1",
        #                     "manufactured_solution_ref2",
        #                     "manufactured_solution_ref3",
        #                     "manufactured_solution_ref4"]
        self.meshes_list = ["manufactured_solution_ref0",
                            "manufactured_solution_ref1",
                            "manufactured_solution_ref2",
                            "manufactured_solution_ref3"]

        self.OriginalProjectParameters = ProjectParameters
        self.print_convergence_plot = self.OriginalProjectParameters["manufactured_solution_settings"]["print_convergence_plot"].GetBool()


    def Solve(self):

        h = []
        err_p = []
        err_v = []

        den = 1
        mesh_0_characteristic_size = 0.2

        # Solve the manufactured solution problem for each one of the refinements
        for mesh_name in self.meshes_list:
            # Solve the problem imposing the previously obtained values
            CaseProjectParameters = self.OriginalProjectParameters.Clone()
            FluidProblem = ManufacturedSolutionProblem(CaseProjectParameters, mesh_name)
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
            h_1_p = []
            h_2_p = []

            den = 1.0
            for i in range(0,len(err_p)):
                h_1_v.append(err_v[0]/den)
                h_2_v.append(err_v[0]/den**2)
                h_1_p.append(err_p[0]/den)
                h_2_p.append(err_p[0]/den**2)
                den *= 2

            plt.loglog(h, err_p, '-+', color = 'r', label = 'Obtained pressure')
            plt.loglog(h, h_1_p, '--', color = 'r', label = 'Linear pressure')
            plt.loglog(h, h_2_p, ':' , color = 'r', label = 'Quadratic pressure')
            plt.loglog(h, err_v, '-x', color = 'k', label = 'Obtained velocity')
            plt.loglog(h, h_1_v, '--', color = 'k', label = 'Linear velocity')
            plt.loglog(h, h_2_v, ':' , color = 'k', label = 'Quadratic velocity')

            plt.title('L2 norm convergence')
            plt.ylabel('Absolute error')
            plt.xlabel('h')
            # plt.xlim([0.02, 0.25])
            # plt.ylim([1e-4, 1])
            plt.legend(loc=0)
            plt.savefig('l2_norm_convergence.png')

        # Check obtained solution
        expected_velocity_errors = [0.020910246816825257, 0.0062279017039999045, 0.0014846307453335115, 0.0003540805601027302, 8.621417044815537e-05]
        expected_pressure_errors = [46.48407227368183, 4.678777003089299, 0.8570316463968392, 0.2160365355817885, 0.06642008924417026]

        for i in range(len(self.meshes_list)):
            self.assertAlmostEqual(err_v[i], expected_velocity_errors[i])
            self.assertAlmostEqual(err_p[i], expected_pressure_errors[i])

class ManufacturedSolutionProblem:

    def __init__(self, ProjectParameters, input_file_name):

        self.ProjectParameters = ProjectParameters
        self.print_output = self.ProjectParameters["manufactured_solution_settings"]["print_output"].GetBool()
        self.problem_type = self.ProjectParameters["manufactured_solution_settings"]["problem_type"].GetString()
        self.analytical_solution_type = self.ProjectParameters["manufactured_solution_settings"]["analytical_solution_type"].GetString()
        self.input_file_name = input_file_name

    def SetFluidProblem(self):

        ## Set the current mesh case problem info
        if (self.problem_type == "analytical_solution"):
            self.ProjectParameters["problem_data"]["problem_name"].SetString(self.input_file_name+"_manufactured")
        else:
            self.ProjectParameters["problem_data"]["problem_name"].SetString(self.input_file_name)
        self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].SetString("ManufacturedSolutionTest/"+self.input_file_name)

        ## Fluid model part definition
        self.main_model_part = KratosMultiphysics.ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        ###TODO replace this "model" for real one once available
        Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString() : self.main_model_part}

        ## Solver construction
        import python_solvers_wrapper_fluid
        self.solver = python_solvers_wrapper_fluid.CreateSolver(self.main_model_part, self.ProjectParameters)

        self.solver.AddVariables()

        ## Read the model - note that SetBufferSize is done here
        self.solver.ImportModelPart()

        ## Add AddDofs
        self.solver.AddDofs()

        ## Initialize GiD  I/O
        if (self.print_output):
            from gid_output_process import GiDOutputProcess
            self.gid_output = GiDOutputProcess(self.solver.GetComputingModelPart(),
                                               self.ProjectParameters["problem_data"]["problem_name"].GetString() ,
                                               self.ProjectParameters["output_configuration"])

            self.gid_output.ExecuteInitialize()

        ## Get the list of the skin submodel parts in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["skin_parts"].size()):
            skin_part_name = self.ProjectParameters["solver_settings"]["skin_parts"][i].GetString()
            Model.update({skin_part_name: self.main_model_part.GetSubModelPart(skin_part_name)})
        for i in range(self.ProjectParameters["solver_settings"]["no_skin_parts"].size()):
            no_skin_part_name = self.ProjectParameters["solver_settings"]["no_skin_parts"][i].GetString()
            Model.update({no_skin_part_name: self.main_model_part.GetSubModelPart(no_skin_part_name)})

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
        step = 0
        out = 0.0

        if (self.print_output):
            self.gid_output.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            Dt = self.solver.ComputeDeltaTime()
            step += 1
            time += Dt
            self.main_model_part.CloneTimeStep(time)
            self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = step

            if (self.print_output):
                self.gid_output.ExecuteInitializeSolutionStep()

            if (self.problem_type == "analytical_solution"):
                # Fix the manufactured solution values (only for visualization purposes)
                self.SetManufacturedSolutionValues(fix=True, set_only_boundaries=False)
            else:
                # Set the manufactured solution source terms
                self.SetManufacturedSolutionValues(fix=True, set_only_boundaries=True)
                self.SetManufacturedSolutionSourceValues()

            if (step < 3):
                self.SetManufacturedSolutionValues(False) # Set the analytical solution in the two first steps
            else:
                if (self.problem_type != "analytical_solution"):
                    self.solver.Solve()

            if (self.print_output):
                self.gid_output.ExecuteFinalizeSolutionStep()

                if self.gid_output.IsOutputStep():
                    self.gid_output.PrintOutput()

            out = out + Dt

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
