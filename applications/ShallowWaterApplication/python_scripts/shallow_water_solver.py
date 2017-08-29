from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShallowWaterApplication import *
CheckForPreviousImport()


def AddVariables(model_part):  #this way er only need one command to add all the variables to our problem 
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY);
    model_part.AddNodalSolutionStepVariable(DELTA_VELOCITY);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(PROJECTED_HEIGHT);
    model_part.AddNodalSolutionStepVariable(DELTA_HEIGHT);
    # Physic problem parameters
    model_part.AddNodalSolutionStepVariable(BATHYMETRY);
    model_part.AddNodalSolutionStepVariable(GRAVITY);
    model_part.AddNodalSolutionStepVariable(RAIN);
    # Specific variables to convect particles
    model_part.AddNodalSolutionStepVariable(YP);
    model_part.AddNodalSolutionStepVariable(MEAN_SIZE);

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(HEIGHT);
    print ("variables for the SWE solver added correctly")

class ShallowWaterSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):  # Constructor of the class 
        self.model_part = model_part
        self.domain_size = domain_size
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # Definition of the mesh stage solver
        self.swe_linear_solver =  SkylineLUFactorizationSolver()  # We set the type of solver that we want 
        
        # definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)  # Tolerance for the solver 
        
        # For the pfem2
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= FindElementalNeighboursProcess(model_part,domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()
    
    #######################################################################
    def Initialize(self):
        # Creating the solution strategy for the mesh stage
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.mesh_solver = strategy_python.SolvingStrategyPython(self.model_part,
                                            self.time_scheme, self.swe_linear_solver,
                                            self.conv_criteria, CalculateReactionFlag,
                                            ReformDofSetAtEachStep, MoveMeshFlag)

        # Creating the solution strategy for the particle stage
        self.VariableUtils = VariableUtils()
        maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = MoveShallowWaterParticleUtility(self.model_part,maximum_number_of_particles)  
        self.moveparticles.MountBin()
        
        # Initialize dry/wet state utility
        #~ self.drybedutility = DryBedUtility(self.model_part)

    #######################################################################   
    def Solve(self):
        # Check dry bed
        #~ (self.drybedutility).CheckConservedVariables(self.model_part.Nodes)
        
        # Move particles
        (self.moveparticles).CalculateVelOverElemSize();
        (self.moveparticles).MoveParticles();
        pre_minimum_number_of_particles=self.domain_size;
        (self.moveparticles).PreReseed(pre_minimum_number_of_particles);    
        (self.moveparticles).TransferLagrangianToEulerian();
        (self.VariableUtils).CopyScalarVar(PROJECTED_HEIGHT,HEIGHT,self.model_part.Nodes)
        (self.VariableUtils).CopyVectorVar(PROJECTED_VELOCITY,VELOCITY,self.model_part.Nodes)
        (self.moveparticles).ResetBoundaryConditions()
        #(self.moveparticles).CopyScalarVarToPreviousTimeStep(PROJECTED_HEIGHT,self.model_part.Nodes)
        #(self.moveparticles).CopyVectorVarToPreviousTimeStep(PROJECTED_VELOCITY,self.model_part.Nodes)
        
        # Solve equations on mesh
        (self.mesh_solver).Solve()
        
        # Update particles
        (self.moveparticles).CalculateDeltaVariables();        
        (self.moveparticles).CorrectParticlesWithoutMovingUsingDeltaVariables();
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles); 

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.mesh_solver).SetEchoLevel(level)




