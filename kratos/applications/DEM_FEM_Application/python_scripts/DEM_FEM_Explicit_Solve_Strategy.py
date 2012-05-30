#importing the Kratos Library

#######old ones
#from Kratos import *
#from KratosDEM_FEM_Application import *

from KratosMultiphysics import *
from KratosMultiphysics.DEM_FEM_Application import *


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
    model_part.AddNodalSolutionStepVariable(GRAVITY);

    model_part.AddNodalSolutionStepVariable(RADIUS);
    model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOURS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FAILURE_TYPE)

    model_part.AddNodalSolutionStepVariable(PARTICLE_INERTIA)
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MOMENT_OF_INERTIA)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ROTATION_ANGLE)

    print "Variables for the DEM&FEM Coupled solution added correctly"
    
    

def AddDofs(model_part):    
     for node in model_part.Nodes:
            node.AddDof(DISPLACEMENT_X,REACTION_X);
            node.AddDof(DISPLACEMENT_Y,REACTION_Y);
            node.AddDof(DISPLACEMENT_Z,REACTION_Z);
            node.AddDof(VELOCITY_X,REACTION_X);
            node.AddDof(VELOCITY_Y,REACTION_Y);
            node.AddDof(VELOCITY_Z,REACTION_Z);
     print "Dofs for the DEM&FEM Coupled solution added correctly"


class DynamicStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part               = model_part  
        self.domain_size              = domain_size
        self.damp_type                = 1
        self.damping_ratio            = 0.8
        self.virtual_mass             = True
        self.contact_stiffness_ratio  = 1.00
        self.max_delta_time           = 0.9
        self.CalculateReactionFlag    = True
        self.ComputeFemFemContact     = False
        self.MoveMeshFlag             = True       
        self.structure_linear_solver  = SkylineLUFactorizationSolver()
        self.time_scheme              = ResidualBasedIncrementalUpdateStaticScheme()
        self.builder                  = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
        self.ConvUnbalForceRatio      = 1.0e-5
        self.Particle_If_Cal_Rotate   = 1
        self.Particle_If_Cal_Rotate_Spring = 1

          
    #######################################################################
  
    def CriticalTime(self):
         (self.solver).Initialize();
         print "Calculating Time Step "
         (self.solver).ComputeCriticalTime()  


    def Initialize(self):
        
        #store processInfo variables:
        self.model_part.ProcessInfo.SetValue(DEM_FEM_DAMP_TYPE, self.damp_type)
        self.model_part.ProcessInfo.SetValue(DEM_FEM_DAMP_RATIO,self.damping_ratio)


        self.model_part.ProcessInfo.SetValue(DEM_FEM_DIMENSION,self.domain_size)
        self.model_part.ProcessInfo.SetValue(DEM_FEM_TIME_STEP,self.max_delta_time)
        self.model_part.ProcessInfo.SetValue(DEM_FEM_CONVERGENCE_RATIO,self.ConvUnbalForceRatio)
        self.model_part.ProcessInfo.SetValue(PARTICLE_IF_CAL_ROTATE,self.Particle_If_Cal_Rotate)
        self.model_part.ProcessInfo.SetValue(PARTICLE_IF_CAL_ROTATE_SPRING,self.Particle_If_Cal_Rotate_Spring)
        
        #creating the solution strategy
        self.solver = DEM_FEM_Explicit_Forward_Differential_Strategy(self.model_part, self.domain_size,self.damp_type, self.damping_ratio, self.virtual_mass, self.contact_stiffness_ratio, self.max_delta_time, self.CalculateReactionFlag, self.ComputeFemFemContact, self.MoveMeshFlag, self.structure_linear_solver, self.time_scheme, self.builder)

        
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

  
    
      
      