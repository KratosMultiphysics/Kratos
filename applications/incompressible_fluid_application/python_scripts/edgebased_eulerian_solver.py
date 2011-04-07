#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)

    print "variables for the edgebased incompressible fluid solver added correctly"

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);


class EdgeBasedLevelSetSolver:
    
    def __init__(self,model_part,domain_size,body_force,viscosity,density):

        print "entered in EdgeBasedEulerianSolver python constructor"
        #data of the problem
        self.model_part = model_part
        self.domain_size = domain_size
        self.body_force = body_force
        self.density = density
        self.viscosity = viscosity

        self.use_mass_correction = True

        self.stabdt_pressure_factor = 1.0;
        self.stabdt_convection_factor = 0.01;

        self.redistance_frequency = 5;
        self.step = 0;

        self.extrapolation_layers = 5
        self.tau2_factor = 1.0
        self.edge_detection_angle = 45.0
        self.assume_constant_pressure = False

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()

        pDiagPrecond = DiagonalPreconditioner()
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)




    def Initialize(self):
        print "entered in EdgeBasedEulerianSolver Initialize"
        #build the edge data structure
        if(self.domain_size == 2):
            self.matrix_container = MatrixContainer2D()
        else:
            self.matrix_container = MatrixContainer3D()

        self.matrix_container.ConstructCSRVector(self.model_part)
        self.matrix_container.BuildCSRData(self.model_part)

        ##for 3D problems we need to evaluate the condition's neighbours
        if(self.domain_size == 3):
            self.condition_neighbours_finder = FindConditionsNeighboursProcess(self.model_part,self.domain_size,10)
            self.condition_neighbours_finder.Execute()

        ##constructing the solver
        print "ln82"
        if(self.domain_size == 2):
            self.fluid_solver = FluidSolver2D(self.matrix_container,self.model_part,self.viscosity,self.density,self.body_force,self.use_mass_correction,self.edge_detection_angle,self.stabdt_pressure_factor,self.stabdt_convection_factor,self.edge_detection_angle,self.assume_constant_pressure)
        else:
            print "ln83"
            self.fluid_solver = FluidSolver3D(self.matrix_container,self.model_part,self.viscosity,self.density,self.body_force,self.use_mass_correction,self.edge_detection_angle,self.stabdt_pressure_factor,self.stabdt_convection_factor,self.edge_detection_angle,self.assume_constant_pressure)
            print "ln84"

        self.fluid_solver.Initialize()
        print "ln91"



        print "**********************************************"
        print "finished EdgeBasedLevelSetSolver initialize"

        

    ################################################################
    ################################################################
    def Solve(self):
##        (self.fluid_solver).UpdateFixedVelocityValues()
        (self.fluid_solver).SolveStep1();
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);
        (self.fluid_solver).SolveStep3();
   
    ################################################################
    ################################################################
    def EstimateTimeStep(self,safety_factor,max_Dt):
        dt = (self.fluid_solver).ComputeTimeStep(safety_factor,max_Dt);

        if(dt > max_Dt):
            dt = max_Dt

        print dt

        return dt




        
        

