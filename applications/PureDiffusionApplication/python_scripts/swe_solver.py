#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PureDiffusionApplication import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
CheckForPreviousImport()

def AddVariables(model_part):  #this way er only need one command to add all the variables to our problem 
    model_part.AddNodalSolutionStepVariable(BATHYMETRY);
    model_part.AddNodalSolutionStepVariable(PROJECTED_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PROJECTED_HEIGHT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    model_part.AddNodalSolutionStepVariable(SCALAR_VELOCITY_X);
    model_part.AddNodalSolutionStepVariable(SCALAR_VELOCITY_Y);
    model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED_VELOCITY_X);
    model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED_VELOCITY_Y);
    #~ model_part.AddNodalSolutionStepVariable(POINT_HEAT_SOURCE);
    #~ model_part.AddNodalSolutionStepVariable(TEMPERATURE);

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(HEIGHT);

    print ("variables for the SWE solver added correctly")

class StaticSWESolver:
    #######################################################################
    def __init__(self,model_part,domain_size):  # Constructor of the class 

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        self.swe_linear_solver =  SkylineLUFactorizationSolver()  # We set the type of solver that we want 
        
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)  # Tolerance for the solver 
        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.model_part,
                                                            self.time_scheme,
                                                            self.swe_linear_solver,
                                                            self.conv_criteria,
                                                            CalculateReactionFlag,
                                                            ReformDofSetAtEachStep,
                                                            MoveMeshFlag)
      
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

class BFECCConvectionSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):  # Constructor of the class
        
        self.model_part = model_part
        self.domain_size = domain_size
        
        # Definition of parameters
        #~ self.scalar_var_convected = 1
        
        # For the bfecc
        # Mount the search structure
        self.locator = BinBasedFastPointLocator2D(self.model_part)
        self.locator.UpdateSearchDatabase()
        
    #######################################################################
    def Initialize(self):
        # Convection only
        #~ self.thermal_settings  = (self.model_part.ProcessInfo).GetValue(CONVECTION_DIFFUSION_SETTINGS)
        #~ self.unknown_var = (self.thermal_settings).GetUnknownVariable()
        #~ self.projection_var = (self.thermal_settings).GetProjectionVariable()
        #~ self.velocity_var = (self.thermal_settings).GetVelocityVariable()

        
        # Now tools for pfem2:
        self.VariableUtils = VariableUtils()
        
        # Construct the utility to move the points
        if self.domain_size ==2: 
            self.bfecc_utility = BFECCConvection2D(self.locator)
        
        print("Finished Initialize")
        
    #######################################################################
    def Solve(self):
        substepping  = 10.0
        
        # Copy components to scalar auxiliary variables
        for node in self.model_part.Nodes:
            # x-component
            value = node.GetSolutionStepValue(VELOCITY_X)
            node.SetSolutionStepValue(SCALAR_VELOCITY_X,value)
            value = node.GetSolutionStepValue(PROJECTED_VELOCITY_X)
            node.SetSolutionStepValue(SCALAR_PROJECTED_VELOCITY_X,value)
            # y-component
            value = node.GetSolutionStepValue(VELOCITY_Y)
            node.SetSolutionStepValue(SCALAR_VELOCITY_Y,value)
            value = node.GetSolutionStepValue(PROJECTED_VELOCITY_Y)
            node.SetSolutionStepValue(SCALAR_PROJECTED_VELOCITY_Y,value)
        
        # Copy unknown to projection_unknown
        (self.VariableUtils).CopyScalarVar(HEIGHT,PROJECTED_HEIGHT,self.model_part.Nodes)
        (self.VariableUtils).CopyScalarVar(SCALAR_VELOCITY_X,SCALAR_PROJECTED_VELOCITY_X,self.model_part.Nodes)
        #~ (self.VariableUtils).CopyScalarVar(VELOCITY_Y,PROJECTED_VELOCITY_Y,self.model_part.Nodes)
        # Copy projection_unknown to previous time step
        self.bfecc_utility.CopyScalarVarToPreviousTimeStep(self.model_part,PROJECTED_HEIGHT)
        self.bfecc_utility.CopyScalarVarToPreviousTimeStep(self.model_part,SCALAR_PROJECTED_VELOCITY_X)
        self.bfecc_utility.CopyScalarVarToPreviousTimeStep(self.model_part,SCALAR_PROJECTED_VELOCITY_Y)
        # Convect projection_unknown from previous time step to current time step
        self.bfecc_utility.BFECCconvect(self.model_part,PROJECTED_HEIGHT,VELOCITY,substepping)
        self.bfecc_utility.BFECCconvect(self.model_part,SCALAR_PROJECTED_VELOCITY_X,VELOCITY,substepping)
        #~ self.bfecc_utility.BFECCconvect(self.model_part,PROJECTED_VELOCITY_Y,VELOCITY,substepping)
        
        # Return auxiliar variables value to vector components
        for node in self.model_part.Nodes:
            # x-component
            value = node.GetSolutionStepValue(SCALAR_PROJECTED_VELOCITY_X)
            node.SetSolutionStepValue(PROJECTED_VELOCITY_X,value)
            # y-component
            value = node.GetSolutionStepValue(SCALAR_PROJECTED_VELOCITY_Y)
            node.SetSolutionStepValue(PROJECTED_VELOCITY_Y,value)
        
        
        
