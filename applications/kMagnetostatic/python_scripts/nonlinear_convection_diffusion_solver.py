#importing the Kratos Library
from Kratos import *
from KratosR1ElectrostaticApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(ELECTROSTATIC_POTENTIAL);
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(HEAT_FLUX);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
    model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX);
    model_part.AddNodalSolutionStepVariable(ELECTRICAL_PERMITTIVITY);
    model_part.AddNodalSolutionStepVariable(ELECTROSTATIC_POINT_CHARGE);

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(ELECTROSTATIC_POTENTIAL);

    print "variables for the electromagnetic solver added correctly"


class ConvectionDiffusionSolver:
    
    def __init__(self,model_part,domain_size):

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.time_order = 1;
	self.prediction_order = 1;
        self.ReformDofAtEachIteration = False;
        self.max_iter = 100
        self.toll = 1e-6

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        pILUPrecond = ILU0Preconditioner()
##        self.linear_solver =  BICGSTABSolver(1e-8, 5000,pILUPrecond)


    def Initialize(self):
        (self.neighbour_search).Execute()
        
        self.solver = ResidualBasedConvectionDiffusionStrategyNonLinear(self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.time_order,self.max_iter,self.toll)   
        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
                 
   
    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            (self.neighbour_search).Execute()        
        
        (self.solver).Solve()

        if(self.ReformDofAtEachIteration == True):
            (self.solver).Clear()      
       

