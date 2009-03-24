#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS);
    print "variables for the incompressible fluid solver added correctly"

def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(PRESSURE);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
    print "dofs for the incompressible fluid solver added correctly"
    

class IncompressibleFluidSolver:
    
    def __init__(self,model_part,domain_size):

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.vel_toll = 0.001;
        self.press_toll = 0.001;
        self.max_vel_its = 4;
        self.max_press_its = 3;
        self.time_order = 2;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = False; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False;

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        pILUPrecond = ILU0Preconditioner()
##        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)


    def Initialize(self):
        (self.neighbour_search).Execute()
        
#        self.solver = ResidualBasedFluidStrategyCoupled(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        print "in python: okkio using Coupled Strategy"
        self.solver = ResidualBasedFluidStrategyCoupled(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   

        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
        
   
    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            (self.neighbour_search).Execute()

        print "just before solve"
        (self.solver).Solve()

