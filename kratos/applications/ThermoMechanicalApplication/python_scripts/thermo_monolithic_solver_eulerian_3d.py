#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *

# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def AddVariables(model_part,settings):
    model_part.AddNodalSolutionStepVariable(settings.GetConvectionVariable());
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);    
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());   
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());    
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());
    #model_part.AddNodalSolutionStepVariable(CONVECTION_COEFFICIENT);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    print "variables for the MONOLITHIC_SOLVER_EULERIAN added correctly"
        
def AddDofs(model_part,settings):
    for node in model_part.Nodes:
	node.AddDof(settings.GetUnknownVariable());
        
    print "dofs for the THERMO monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size,my_settings):

        self.model_part = model_part

        #self.alpha = -0.3
        #self.move_mesh_strategy = 0
        #self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme( self.alpha,self.move_mesh_strategy )
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        
        self.settings = my_settings
        #definition of the solvers
        #self.linear_solver =  SkylineLUFactorizationSolver()
##        self.linear_solver =SuperLUSolver()
##        self.linear_solver = MKLPardisoSolver()

        pPrecond = DiagonalPreconditioner()
##        pPrecond = ILU0Preconditioner()
        self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)


       # self.conv_criteria = UPCriteria(1e-12,1e-14,1e-15,1e-17)
        #self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.1);

        self.dynamic_tau = 1.0
        #self.oss_switch  = 0

        
        #self.max_iter = 30
                            
        #default settings
        self.echo_level = 0
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
    
        self.neigh_finder = FindNodalNeighboursProcess(self.model_part,9,18)
        self.elem_neighbor_finder = FindElementalNeighboursProcess(self.model_part, 2, 20)
        self.Nmax = len(model_part.Properties)
        self.contact_matrix = Matrix()
        ##calculate normals
        self.normal_tools = BodyNormalCalculationUtils()
        
    #######################################################################
    def Initialize(self):

        self.duplicate_and_create_conditions = DuplicateInterfaceNodesCreateConditionsProcess(self.model_part,"HeatContact3D", self.Nmax, self.contact_matrix )

        (self.neigh_finder).ClearNeighbours();
        (self.neigh_finder).Execute();
        
        (self.elem_neighbor_finder).ClearNeighbours()
        (self.elem_neighbor_finder).Execute() 
        print "INSIDE INITIALIZE"           
 	(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)    
 	
        self.solver = ResidualBasedLinearStrategy(self.model_part,self.time_scheme,self.linear_solver,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.CalculateNormDxFlag,self.MoveMeshFlag)   
        (self.solver).SetEchoLevel(self.echo_level)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);

        #self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        #self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch );
        #self.model_part.ProcessInfo.SetValue(M, self.regularization_coef );
        #(self.duplicate_and_create_conditions).Execute()
        
        self.normal_tools.CalculateBodyNormals(self.model_part,3);        

##        print "Initialization monolithic solver finished"
	                     
    #######################################################################   
    def Solve(self):        
        print "*****************entering solve?????????????"
        (self.solver).Solve()
##        print "solving step monolithic solver finished"
       

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################

        




