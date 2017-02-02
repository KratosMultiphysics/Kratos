from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.MeshingApplication import *

from math import *

######################################################################################
######################################################################################
######################################################################################


class ManufacturedSolutionProblem:

    def __init__(self, input_file_name, settings_file_name, problem_type):

        self.problem_type = problem_type
        self.input_file_name = input_file_name
        self.settings_file_name = settings_file_name

    def SetFluidProblem(self):

        ## Parse the ProjectParameters
        parameter_file = open(self.settings_file_name,'r')
        self.ProjectParameters = Parameters( parameter_file.read())

        ## Set the current mesh case problem info
        if (self.problem_type == "Manufactured solution"):
            self.ProjectParameters["problem_data"]["problem_name"].SetString(self.input_file_name+"_manufactured")
        else:
            self.ProjectParameters["problem_data"]["problem_name"].SetString(self.input_file_name)
        self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].SetString(self.input_file_name)

        ## Fluid model part definition
        self.main_model_part = ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

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
                node.SetSolutionStepValue(DISTANCE, 0, 1.0)

        ## Fix the pressure in one node (bottom left corner)
        for node in self.main_model_part.Nodes:
            if ((node.X<0.001) and (node.Y<0.001)):
                node.Fix(PRESSURE)
                node.SetSolutionStepValue(PRESSURE, 0, 0.0)

    def SolveFluidProblem(self):

        ## Stepping and time settings
        Dt = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        time = 0.0
        step = 0
        out = 0.0

        self.gid_output.ExecuteBeforeSolutionLoop()

        while(time <= end_time):

            time = time + Dt
            step = step + 1
            self.main_model_part.CloneTimeStep(time)

            print("STEP = ", step)
            print("TIME = ", time)

            self.gid_output.ExecuteInitializeSolutionStep()

            if (self.problem_type == "Manufactured solution"):
                # Fix the manufactured solution values (only for visualization purposes)
                self.SetManufacturedSolutionValues(fix=True, set_on_boundary=False)
            else:
                # Set the manufactured solution source terms
                self.SetManufacturedSolutionValues(fix=True, set_on_boundary=True)
                self.SetManufacturedSolutionSourceValues()

            if (step < 3):
                self.SetManufacturedSolutionValues(False) # Set the analytical solution in the two first steps
            else:
                if (self.problem_type != "Manufactured solution"):
                    self.solver.Solve()

            self.gid_output.ExecuteFinalizeSolutionStep()

            if self.gid_output.IsOutputStep():
                self.gid_output.PrintOutput()

            out = out + Dt

        self.gid_output.ExecuteFinalize()

    def SetManufacturedSolutionValues(self, fix=True, set_on_boundary=False):
        ## Set the analytical solution for the manufactured solution computation
        time = self.main_model_part.ProcessInfo[TIME]
        if (set_on_boundary == False):
            for node in self.main_model_part.Nodes:
                x_vel = self.ComputeNodalVelocityXManufacturedSolution(node, time)
                y_vel = self.ComputeNodalVelocityYManufacturedSolution(node, time)
                press = self.ComputeNodalPressureManufacturedSolution(node)

                if (fix == True):
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)
                    node.Fix(PRESSURE)

                node.SetSolutionStepValue(VELOCITY_X, 0, x_vel)
                node.SetSolutionStepValue(VELOCITY_Y, 0, y_vel)
                node.SetSolutionStepValue(PRESSURE, 0, press)
        else:
            for node in self.main_model_part.GetSubModelPart("Inlet2D_Contour").Nodes:
                x_vel = self.ComputeNodalVelocityXManufacturedSolution(node, time)
                y_vel = self.ComputeNodalVelocityYManufacturedSolution(node, time)

                if (fix == True):
                    node.Fix(VELOCITY_X)
                    node.Fix(VELOCITY_Y)

                node.SetSolutionStepValue(VELOCITY_X, 0, x_vel)
                node.SetSolutionStepValue(VELOCITY_Y, 0, y_vel)

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
                aux = node.GetSolutionStepValue(NODAL_AREA)  # Current nodal area (from other elements)
                aux += Area/3.0                              # Accumulate the current element nodal area
                node.SetSolutionStepValue(NODAL_AREA, 0, aux)
                node.SetValue(NODAL_AREA, aux)

        ## Check nodal area computation (squared shaped domain of 1x1 m)
        AreaTotal = 0.0
        for node in self.main_model_part.Nodes:
            # print("Node id "+str(node.Id)+" nodal area: "+str(node.GetValue(NODAL_AREA)))
            AreaTotal += node.GetValue(NODAL_AREA)

        if (abs(1.0-AreaTotal) > 1e-5):
            print("Obtained total area: "+str(AreaTotal))
            raise Exception("Error in NODAL_AREA computation.")

    def SetManufacturedSolutionSourceValues(self):
        ## Set the body force as source term
        time = self.main_model_part.ProcessInfo[TIME]

        for node in self.main_model_part.Nodes:
            rho = node.GetSolutionStepValue(DENSITY)
            # nodal_area = node.GetValue(NODAL_AREA)
            # If VMS2D element is used, set mu as the Kinematic viscosity
            if (self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "Embedded"):
                mu = node.GetSolutionStepValue(DYNAMIC_VISCOSITY)
            elif (self.ProjectParameters["solver_settings"]["solver_type"].GetString() == "Monolithic"):
                mu = rho*node.GetSolutionStepValue(VISCOSITY)

            # Trigonometric transient field obtained with by hand derivation
            # rhofy = -rho*pi*cos(pi*node.X)*sin(pi*node.Y)*cos(pi*time) - 2*mu*pi*pi*cos(pi*node.X)*sin(pi*node.Y)*sin(pi*time) + rho*pi*(sin(pi*time)**2)*sin(pi*node.Y)*cos(pi*node.Y)
            # rhofx =  rho*pi*sin(pi*node.X)*cos(pi*node.Y)*cos(pi*time) + 2*mu*pi*pi*sin(pi*node.X)*cos(pi*node.Y)*sin(pi*time) + rho*pi*(sin(pi*time)**2)*sin(pi*node.X)*cos(pi*node.X)
            rhofx = -rho*pi*sin(pi*node.X)*cos(pi*node.Y)*sin(pi*time) + 2*mu*pi*pi*sin(pi*node.X)*cos(pi*node.Y)*cos(pi*time) + rho*pi*(cos(pi*time)**2)*sin(pi*node.X)*cos(pi*node.X)
            rhofy =  rho*pi*cos(pi*node.X)*sin(pi*node.Y)*sin(pi*time) - 2*mu*pi*pi*cos(pi*node.X)*sin(pi*node.Y)*cos(pi*time) + rho*pi*(cos(pi*time)**2)*sin(pi*node.Y)*cos(pi*node.Y)

            # Trigonometric transient field obtained with symbolic generation
            # rhofx =  2*pi**2*mu*sin(pi*node.X)*sin(pi*time)*cos(pi*node.Y) + rho*(pi*sin(pi*node.X)*sin(pi*node.Y)**2*sin(pi*time)**2*cos(pi*node.X) + pi*sin(pi*node.X)*sin(pi*time)**2*cos(pi*node.X)*cos(pi*node.Y)**2) + pi*rho*sin(pi*node.X)*cos(pi*node.Y)*cos(pi*time)
            # rhofy = -2*pi**2*mu*sin(pi*node.Y)*sin(pi*time)*cos(pi*node.X) + rho*(pi*sin(pi*node.X)**2*sin(pi*node.Y)*sin(pi*time)**2*cos(pi*node.Y) + pi*sin(pi*node.Y)*sin(pi*time)**2*cos(pi*node.X)**2*cos(pi*node.Y)) - pi*rho*sin(pi*node.Y)*cos(pi*node.X)*cos(pi*time)
            # print("Node id. ",node.Id," by hand: ",rhofx_an," ",rhofy_an," sym: ",rhofx," ",rhofy)

            # Non-linear transient field
            # rhofx =   rho*pi*(node.X**2)*node.Y*cos(pi*time) - 2*mu*node.Y*sin(pi*time) + rho*(node.X**3)*(node.Y**2)*((sin(pi*time))**2)
            # rhofy =  -rho*pi*node.X*(node.Y**2)*cos(pi*time) + 2*mu*node.X*sin(pi*time) + rho*(node.X**2)*(node.Y**3)*((sin(pi*time))**2)

            # Non-linear stationary field
            # rhofx = -2*mu*node.Y + rho*node.X**3*node.Y**2
            # rhofy =  2*mu*node.X + rho*node.X**2*node.Y**3

            node.SetSolutionStepValue(BODY_FORCE_X, 0, rhofx/rho)  # Set the x-component body force field
            node.SetSolutionStepValue(BODY_FORCE_Y, 0, rhofy/rho)  # Set the y-component body force field

    def ComputeVelocityErrorNorm(self):
        err_v = 0

        for node in self.main_model_part.Nodes:
            weight = node.GetValue(NODAL_AREA)
            vel_x = node.GetSolutionStepValue(VELOCITY_X)
            vel_y = node.GetSolutionStepValue(VELOCITY_Y)
            end_time = self.main_model_part.ProcessInfo[TIME]

            analytical_vel_x = self.ComputeNodalVelocityXManufacturedSolution(node, end_time)
            analytical_vel_y = self.ComputeNodalVelocityYManufacturedSolution(node, end_time)

            err_x = analytical_vel_x - vel_x
            err_y = analytical_vel_y - vel_y
            err_node = err_x**2 + err_y**2
            err_v += weight*err_node

        return sqrt(err_v) # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    def ComputePressureErrorNorm(self):
        err_p = 0

        for node in self.main_model_part.Nodes:
            weight = node.GetValue(NODAL_AREA)
            pres = node.GetSolutionStepValue(PRESSURE)
            analytical_pres = self.ComputeNodalPressureManufacturedSolution(node)
            err_p += weight*(analytical_pres - pres)**2

        return sqrt(err_p) # Note, there is no need of dividing by the total area (sum of weights) since it is 1

    def ComputeNodalDensityManufacturedSolution(self, node):
        # return 1.0 + 0.1*math.sin(0.75*math.pi*node.X) + 0.15*math.cos(1.0*math.pi*node.Y) + 0.08*math.cos(1.25*math.pi*node.X*node.Y)
        return 1.0

    # def ComputeNodalVelocityXManufacturedSolution(self, node):
    def ComputeNodalVelocityXManufacturedSolution(self, node, time):
        # return sin(pi*node.X)*cos(pi*node.Y)*sin(pi*time)   # Time dependent solution
        return sin(pi*node.X)*cos(pi*node.Y)*cos(pi*time)   # Time dependent solution
        # return node.X**2*node.Y*sin(pi*time)                # Non-linear transient field
        # return node.X**2*node.Y                             # Non-linear stationary field

    # def ComputeNodalVelocityYManufacturedSolution(self, node):
    def ComputeNodalVelocityYManufacturedSolution(self, node, time):
        # return -cos(pi*node.X)*sin(pi*node.Y)*sin(pi*time) # Time dependent solution
        return -cos(pi*node.X)*sin(pi*node.Y)*cos(pi*time) # Time dependent solution
        # return -node.X*node.Y**2*sin(pi*time)               # Non-linear transient field
        # return -node.X*node.Y**2                            # Non-linear stationary field

    def ComputeNodalPressureManufacturedSolution(self, node):
        # return 100000.0 - 30000.0*math.cos(1.0*math.pi*node.X) + 20000.0*math.sin(1.25*math.pi*node.Y) - 25000.0*math.sin(0.75*math.pi*node.X*node.Y)
        return 0.0


## List of cases to solve
meshes_list = ["manufactured_solution_ref0",
              "manufactured_solution_ref1",
              "manufactured_solution_ref2",
              "manufactured_solution_ref3",
              "manufactured_solution_ref4"]
# meshes_list = ["manufactured_solution_ref0",
#                "manufactured_solution_ref1",
#                "manufactured_solution_ref2",
#                "manufactured_solution_ref3"]

h = []
err_p = []
err_v = []

den = 1
for mesh in meshes_list:
    # Obtain the pressure and reaction for the given manufactured solution
    FluidProblemManufactured = ManufacturedSolutionProblem(mesh, "ProjectParametersManufactured.json", "Manufactured solution")
    FluidProblemManufactured.SetFluidProblem()
    FluidProblemManufactured.SolveFluidProblem()

    # Solve the problem imposing the previously obtained values
    FluidProblem = ManufacturedSolutionProblem(mesh, "ProjectParameters.json", "Source")
    FluidProblem.SetFluidProblem()
    FluidProblem.SolveFluidProblem()

    # Compute the obtained solution error
    h.append(0.2/den)
    err_p.append(FluidProblem.ComputePressureErrorNorm())
    err_v.append(FluidProblem.ComputeVelocityErrorNorm())
    den *= 2

print(err_p)
print(err_v)

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

# print(h_1_v)
# print(h_2_v)
# print(h_1_p)
# print(h_2_p)

# m_1_p = (log(h_1_p[0])-log(h_1_p[-1]))/(log(h[0])-log(h[-1]))
# m_2_v = (log(h_2_v[0])-log(h_2_v[-1]))/(log(h[0])-log(h[-1]))
#
# print(m_1_p)
# print(m_2_v)

m_avg_p = (log(err_p[0])-log(err_p[-1]))/(log(h[0])-log(h[-1]))
m_avg_v = (log(err_v[0])-log(err_v[-1]))/(log(h[0])-log(h[-1]))

print(m_avg_p)
print(m_avg_v)

plt.loglog(h, err_v,'-x', color ='k', label = 'Velocity')
plt.loglog(h, err_p,'-+', color ='k', label = 'Pressure')
plt.loglog(h,h_1_p, label = 'Pressure (exactly linear)')
plt.loglog(h,h_2_p, label = 'Pressure (exactly quadratic)')
plt.loglog(h,h_1_v, label = 'Velocity (exactly linear)')
plt.loglog(h,h_2_v, label = 'Velocity (exactly quadratic)')

plt.title('L2 norm convergence')
plt.ylabel('Absolute error')
plt.xlabel('h')
# plt.xlim([0.02, 0.25])
# plt.ylim([1e-4, 1])
plt.legend(loc=0)
plt.savefig('l2_norm_convergence.png')
