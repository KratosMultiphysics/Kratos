#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *



def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(RHS_VECTOR);
    model_part.AddNodalSolutionStepVariable(AUX_VECTOR);
    
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);

    model_part.AddNodalSolutionStepVariable(NORMAL);
    #this one is for Proj Dirichlet Conds
    model_part.AddNodalSolutionStepVariable(AUX_VEL);

    print "variables for the Runge Kutta Frac Step GLS solver added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);
        node.AddDof(AUX_VEL_X);
        node.AddDof(AUX_VEL_Y);
        node.AddDof(AUX_VEL_Z);

        node.AddDof(PRESSURE);
        
    print "dofs for the the Runge Kutta Frac Step GLS solver added correctly"

class RungeKuttaFracStepCompSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        #self.move_mesh_strategy = 2
        pDiagPrecond = DiagonalPreconditioner()

        #definition of the solvers
        self.linear_solver =  SkylineLUFactorizationSolver()
        #self.linear_solver = CGSolver(1e-3, 5000,pDiagPrecond)
        #self.conv_criteria = UPCriteria(1e-3,1e-9,1e-3,1e-6)

        self.max_iter = 10
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.domain_size = 2;
        #self.MoveMeshFlag = True

        self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)

        ##calculate normals
        self.normal_tools = NormalCalculationUtils()

        #self.Hfinder  = FindNodalHProcess(model_part);
        
    

        

        
    #######################################################################
    def Initialize(self):
        #calculate the normals to the overall domain
        self.normal_tools.CalculateOnSimplex(self.model_part.Conditions,self.domain_size);
        #for SLIP condition we need to save these Conditions in a list
        #by now SLIP conditions are identified by FLAG_VARIABLE=3.0. this is done in the constructir of the strategy
        
        #creating the solution strategy
        self.solver = RungeKuttaFracStepCompStrategy(self.model_part,self.linear_solver,self.CalculateReactionFlag,
                                                  self.ReformDofSetAtEachStep, self.CalculateNormDxFlag, self.domain_size )   
        (self.solver).SetEchoLevel(self.echo_level)

        (self.neigh_finder).Execute();
        #(self.Hfinder).Execute();
        
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()
        #(self.neigh_finder).Execute();
        #self.Remesh()

           

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    #######################################################################

    
