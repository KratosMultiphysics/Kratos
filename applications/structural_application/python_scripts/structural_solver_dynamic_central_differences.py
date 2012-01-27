#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(RHS);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DAMAGE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(NODAL_VALUES);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_CONTACT_SLAVE);
    model_part.AddNodalSolutionStepVariable(IS_CONTACT_MASTER);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_VOLUME);
    model_part.AddNodalSolutionStepVariable(NODAL_DAMAGE);
    model_part.AddNodalSolutionStepVariable(IS_DUPLICATED);
    model_part.AddNodalSolutionStepVariable(NODAL_STRAIN);
    model_part.AddNodalSolutionStepVariable(REFINEMENT_LEVEL);
    model_part.AddNodalSolutionStepVariable(SPLIT_NODAL);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    #model_part.AddNodalSolutionStepVariable(KINETIC_ENERGY);
    #model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_STRESS);
    print "Variables for the dynamic structural solution added correctly"
    
    

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X,REACTION_X);
        node.AddDof(DISPLACEMENT_Y,REACTION_Y);
        node.AddDof(DISPLACEMENT_Z,REACTION_Z);
    print "Dofs for the dynamic structural solution added correctly"



class DynamicStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part               = model_part  
        self.domain_size              = domain_size
        self.damping_ratio            = 0.3;
        self.penalty_factor           = 50.00  
        self.max_delta_time           = 0.05;
        self.fraction_delta_time      = 0.90;
        self.CalculateReactionFlag    = True;
        self.MoveMeshFlag             = True;
        self.ComputeContactConditions = False;
        self.CE                       = Constraint_Enforcement.Penalty_Methods; #Constraint_Enforcement.Lagrange_Multiplie_Methods;

               
    #######################################################################
  
    def CriticalTime(self):
         (self.solver).Initialize();
         print "Calculating Time Step "
         (self.solver).ComputeCriticalTime()  


    def Initialize(self):
        
        #creating the solution strategy
        self.solver = ResidualBasedCentralDiferencesStrategy(self.model_part, self.CE, self.domain_size,  self.damping_ratio, self.fraction_delta_time, self.max_delta_time, self.penalty_factor, self.CalculateReactionFlag, self.ComputeContactConditions, self.MoveMeshFlag)

        
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
        
    def SetFractionDeltaTime(self, fraction):
      (self.solver).SetFractionDeltaTime(fraction)
      
    def SetConditionsFlag(self, ComputeContactConditions):
      (self.solver).SetConditionsFlag(ComputeContactConditions)
      
    def CalculateBoundaryContours(self, ComputeBoundary):
      (self.solver).CalculateBoundaryContours(ComputeBoundary)
    
    
      
      