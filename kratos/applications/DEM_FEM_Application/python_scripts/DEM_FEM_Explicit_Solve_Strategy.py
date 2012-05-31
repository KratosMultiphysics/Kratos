#importing the Kratos Library

#######old ones
#from Kratos import *
#from KratosDEM_FEM_Application import *

from KratosMultiphysics import *
from KratosMultiphysics.DEM_FEM_Application import *


def AddVariables(model_part):

    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(APPLIED_FORCE)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);

    model_part.AddNodalSolutionStepVariable(PARTICLE_NUMBER_OF_NEIGHBOURS)

    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(RADIUS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_DENSITY)
    model_part.AddNodalSolutionStepVariable(PARTICLE_STIFFNESS)
    model_part.AddNodalSolutionStepVariable(YOUNG_MODULUS)
    model_part.AddNodalSolutionStepVariable(POISSON_RATIO)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COEF_RESTITUTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_ZETA)
    model_part.AddNodalSolutionStepVariable(PARTICLE_MATERIAL)
    model_part.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COHESION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_TENSION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_LOCAL_DAMP_RATIO)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FAILURE_ID)

    
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
        self.gravity                  = (0.0,-9.81,0.0)

         #type of problem:

        self.delta_OPTION                   = True
        self.continuum_simulating_OPTION    = True
        self.case_OPTION                    = 0  #aixo es una xapuza fins que pooyan permeti bools a pyton o tinguis flags.

        self.rotation_OPTION                = 1 #CANVIAR AIXO, ENTRADA GID
        self.rotation_spring_OPTION         = 1 #CANVIAR AIXO, ENTRADA GID

        #problem specific parameters

        self.force_calculation_type_id      =1
        self.damp_id                        =1
        self.search_radius_extension        = 0.0

        self.dummy_switch                   =0

          
    #######################################################################
  
    def CriticalTime(self):
         (self.solver).Initialize();
         print "Calculating Time Step "
         (self.solver).ComputeCriticalTime()  


    def Initialize(self):


        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)

        if(self.delta_OPTION==True):
            if(self.continuum_simulating_OPTION==True): self.case_OPTION = 2
            else: self.case_OPTION = 1
        elif(self.delta_OPTION==False):
            if(self.continuum_simulating_OPTION==False): self.case_OPTION = 0
            else: self.case_OPTION = 3

        #store processInfo variables:
        # Miquel
        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)    #M: = a type_id
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)                                   #M: = a damp_type
        self.model_part.ProcessInfo.SetValue(SEARCH_RADIUS_EXTENSION, self.search_radius_extension)
        self.model_part.ProcessInfo.SetValue(DUMMY_SWITCH, self.dummy_switch)

        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_OPTION)
        self.model_part.ProcessInfo.SetValue(ROTATION_OPTION, self.rotation_OPTION)
        self.model_part.ProcessInfo.SetValue(ROTATION_SPRING_OPTION, self.rotation_spring_OPTION)

        # Cfeng      
        self.model_part.ProcessInfo.SetValue(DEM_DELTA_TIME,self.max_delta_time)
        self.model_part.ProcessInfo.SetValue(DEM_FEM_CONVERGENCE_RATIO,self.ConvUnbalForceRatio)
     
        #creating the solution strategy
        self.solver = DEM_FEM_Explicit_Forward_Differential_Strategy(self.model_part, self.domain_size,self.damp_type, self.damping_ratio, self.virtual_mass, self.contact_stiffness_ratio, self.max_delta_time, self.CalculateReactionFlag, self.ComputeFemFemContact, self.MoveMeshFlag, self.structure_linear_solver, self.time_scheme, self.builder)
     
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

  
    
      
      