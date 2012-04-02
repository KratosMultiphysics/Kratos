from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(MACH_NUMBER)
    model_part.AddNodalSolutionStepVariable(PRESSURE_COEFFICIENT)

    print "variables for the edgebased incompressible fluid solver added correctly"

def AddDofs(model_part):
      print "dofs for the edgebased incompressible fluid solver added correctly"
    
##def ReadRestartFile(FileName,nodes):
##   aaa = __import__(FileName)
##   aaa.Restart(nodes)

def ReadRestartFile(FileName,nodes):
   NODES = nodes
   aaa = open(FileName)
   for line in aaa:
       exec(line)
       
##   import start.pyinc
   
##   aaa = __import__(FileName)
##   aaa.Restart(nodes)

   

class IncompressibleFluidSolver:
    
    def __init__(self,model_part,domain_size,matrix_container):
        #data of the problem
        self.model_part = model_part
        self.domain_size = domain_size

        #edge data
        self.matrix_container = matrix_container

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()

        #default values for user-defined parameters 
        self.include_shock_capturing = False
        self.smooth_convection = False
        self.CFLnumber = 0.9

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.pressure_linear_solver =  CGSolver(1e-3, 5000,pDiagPrecond)


    def Initialize(self):
        
        
   
    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            (self.neighbour_search).Execute()

        print "just before solve"
        (self.solver).Solve()



        
        

