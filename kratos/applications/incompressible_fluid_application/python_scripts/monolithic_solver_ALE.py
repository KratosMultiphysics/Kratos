#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE); 
    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
        
    print "dofs for the monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
##        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.move_mesh_strategy = 0
        self.time_scheme = ResidualBasedLagrangianMonolithicScheme( self.move_mesh_strategy)

        #definition of the solvers
        self.linear_solver =  SkylineLUFactorizationSolver()
        
        #definition of the convergence criteria
##        self.conv_criteria = UPCriteria(1e-7,1e-9,1e-7,1e-9)
        self.conv_criteria = UPCriteria(1e-3,1e-9,1e-3,1e-6)

        self.max_iter = 10
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        ####MESH CHANGES
##        self.UlfUtils = UlfUtils()
##        self.ulf_time_step_dec_process = UlfTimeStepDecProcess(model_part);
##        self.mark_close_nodes_process = MarkCloseNodesProcess(model_part);
##        self.node_erase_process = NodeEraseProcess(model_part);
##        self.h_multiplier = 0.8
        
##        self.Mesher = TriGenPFEMModeler()
##        self.Mesher = TriGenModeler()
##        self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)

        #detect initial size distribution - note that initially the fluid model part contains
        #all the elements of both structure and fluid ... this is only true after reading the input
##        (self.neigh_finder).Execute();

##        self.remeshing_flag = True
##        self.alpha_shape = 1.5
##        for node in self.model_part.Nodes:
##            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
##                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)
        #U NEED IT FOR ALPHA-shape
##        Hfinder  = FindNodalHProcess(model_part);
##        Hfinder.Execute();
    

        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        
        self.solver = NewtonRaphsonOssStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
        (self.solver).SetEchoLevel(self.echo_level)

##        (self.neigh_finder).Execute();
##        self.Remesh()        
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()
##        (self.neigh_finder).Execute();
##        self.Remesh()

    def EstimateDeltaTime(self,max_dt,domain_size):
        print "Estimating delta time"
        return (self.ulf_time_step_dec_process).EstimateDeltaTime(max_dt,domain_size)


        

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################
##    def Remesh(self):
        #reset the flag
##        for node in self.model_part.Nodes:
##            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)
##
##        for node in self.model_part.Nodes:            
##            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
##                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)
        #and erase bad nodes
##        (self.mark_close_nodes_process).MarkCloseNodes(self.h_multiplier);
        #(self.UlfUtils).MarkNodesCloseToWall(self.model_part, 2, self.alpha_shape)
             
        ##erase all conditions and elements prior to remeshing
##        if (self.remeshing_flag==True):
##            print "DI"
##            ((self.model_part).Elements).clear();
##            ((self.model_part).Conditions).clear();            
## 
##            (self.node_erase_process).Execute();
##
##        if (self.remeshing_flag==True):
####            (self.Mesher).ReGeneratePFEMASGS(self.model_part,self.alpha_shape)
##            (self.Mesher).ReGenerateASGS(self.model_part,self.alpha_shape)
##
####            (self.node_erase_process).Execute();
##            ##calculating fluid neighbours before applying boundary conditions
##            (self.neigh_finder).Execute();
##        for node in self.model_part.Nodes:
##            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)
##
##        for node in self.model_part.Nodes:            
##            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
##                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)
            
        ##calculating fluid neighbours before applying boundary conditions
        #(self.neigh_finder).Execute();
        ##################################################################
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
        




