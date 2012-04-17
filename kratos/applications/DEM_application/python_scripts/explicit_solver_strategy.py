
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOURS)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_STIFFNESS)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COEF_RESTITUTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ZETA)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)

    print "variables for the explicit solver added correctly"

def AddDofs(model_part):
    print "Dofs added"
 
class ExplicitStrategy:
    
    def __init__(self,model_part,domain_size):

        self.model_part               = model_part  
        self.domain_size              = domain_size
        self.damping_ratio            = 0.00;
        self.penalty_factor           = 10.00  
        self.max_delta_time           = 0.05;
        self.fraction_delta_time      = 0.90;
        self.MoveMeshFlag             = True;
        self.time_scheme              = ResidualBasedIncrementalUpdateStaticScheme()
        
    #######################################################################
  

    def Initialize(self):      
        #creating the solution strategy
        self.solver = ExplicitSolverStrategy(self.model_part, self.domain_size,  self.damping_ratio, self.fraction_delta_time, self.max_delta_time,  self.MoveMeshFlag, self.time_scheme)

        
    #######################################################################   
    def Solve(self):    
        (self.solver).Solve()
    
    