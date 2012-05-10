
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
# Check that KratosMultiphysics was imported in the main script
#CheckForPreviousImport(

def AddVariables(model_part):
    #model_part.AddNodalSolutionStepVariable(NUMBER_OF_NEIGHBOURS)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(APPLIED_FORCE)
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
    model_part.AddNodalSolutionStepVariable(PARTICLE_CONTINUUM)
    model_part.AddNodalSolutionStepVariable(PARTICLE_COHESION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FRICTION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_TENSION)
    model_part.AddNodalSolutionStepVariable(PARTICLE_LOCAL_DAMP_RATIO)
    model_part.AddNodalSolutionStepVariable(PARTICLE_FAILURE_ID)

    print "variables for the explicit solver added correctly"

def AddDofs(model_part):
    
    for node in model_part.Nodes:
    #adding dofs
      node.AddDof(DISPLACEMENT_X,REACTION_X);
      node.AddDof(DISPLACEMENT_Y,REACTION_Y);
      node.AddDof(DISPLACEMENT_Z,REACTION_Z);
    print "dofs for the dynamic structural solution added correctly"

 
class ExplicitStrategy:
    
    def __init__(self,model_part,domain_size):

        self.model_part                     = model_part
        self.domain_size                    = domain_size
        self.damping_ratio                  = 0.00;
        self.penalty_factor                 = 10.00
        self.max_delta_time                 = 0.05;
        self.fraction_delta_time            = 0.90;
        self.MoveMeshFlag                   = True;
        self.time_scheme                    = FowardEulerScheme();
        self.gravity                        = (0.0,-9.81,0.0)
        self.delta_time                     = 0.0001;

        #type of problem:

        self.delta_OPTION                   = False
        self.continuum_simulating_OPTION    = False
        self.case_OPTION                    = 0  #aixo es una xapuza fins que pooyan permeti bools a pyton o tinguis flags.

        #problem specific parameters

        self.force_calculation_type_id      =1
        self.damp_id                        =1
        self.solver_id                      =1
        self.search_radius_extension        = 0.0

      
    #######################################################################
  

    def Initialize(self):

        #definir les variables del ProcessInfo:
        self.model_part.ProcessInfo.SetValue(GRAVITY, self.gravity)
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, self.delta_time)

        #POOOYAAAAN NO EM VAN AKESTS NO ELS CONEIX PYTON DE BOOL A BOOL

        ##################XAPUZA##XAPUZA##XAPUZA#############################


        #LA CONTINUUM OPTION SHA DE LLEGIR MIRANT SI TOTS ELS CONTINUUM GROUPS SON 0 LLAVORS LA OPTION ES FALSE
        #SI HI HA ALGUN KE ES DIFERENT DE 0 LLAVORS LA OPTION ES TRUE.
        
        if(self.delta_OPTION==True):
            if(self.continuum_simulating_OPTION==True): self.case_OPTION = 2
            else: self.case_OPTION = 1
        elif(self.delta_OPTION==False):
            if(self.continuum_simulating_OPTION==False): self.case_OPTION = 0
            else: self.case_OPTION = 3

        self.model_part.ProcessInfo.SetValue(CASE_OPTION, self.case_OPTION)

        
        #####

        self.model_part.ProcessInfo.SetValue(FORCE_CALCULATION_TYPE, self.force_calculation_type_id)    #M: = a type_id
        self.model_part.ProcessInfo.SetValue(DAMP_TYPE, self.damp_id)                                   #M: = a damp_type
        self.model_part.ProcessInfo.SetValue(SEARCH_RADIUS_EXTENSION, self.search_radius_extension)
        #self.model_part.ProcessInfo.SetValue(SOLVER_ID, self.solver_id) AIXO NO ES TA FET PERO NOSE SI SA DE FER AIXI O NO.


        #creating the solution strategy
        self.solver = ExplicitSolverStrategy(self.model_part, self.domain_size,  self.damping_ratio, self.fraction_delta_time, self.delta_time,
                                            self.MoveMeshFlag, self.delta_OPTION, self.continuum_simulating_OPTION, self.time_scheme)
        #self.solver.Check() #es sa fer sempre un check despres de montar una estrategia.
        self.solver.Initialize() #aqui definirem el initialize dels elements pero tamb funcions que vulguem fer en el primer pras.
       
           
              
            
        #self.solver.SetCohesiveContacts()



    #######################################################################   
    def Solve(self):
        (self.solver).Solve()
    
    