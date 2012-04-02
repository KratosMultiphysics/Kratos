from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.PFEMApplication import *
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_POROUS);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(THAWONE);
    model_part.AddNodalSolutionStepVariable(THAWTWO); 
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE); 
    model_part.AddNodalSolutionStepVariable(ARRHENIUS); 
    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
	node.AddDof(AIR_PRESSURE,REACTION_WATER_PRESSURE);
        
    print "dofs for the monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size,box_corner1,box_corner2):

        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2

        self.model_part = model_part
        self.move_mesh_strategy = 2
        self.time_scheme = ResidualBasedLagrangianMonolithicScheme( self.move_mesh_strategy)

        #definition of the solvers
        self.linear_solver =  SkylineLUFactorizationSolver()
        
        #definition of the convergence criteria
        self.conv_criteria = UPCriteria(1e-7,1e-9,1e-7,1e-9)
      #  self.conv_criteria = UPCriteria(1e-12,1e-14,1e-9,1e-10)

        self.max_iter = 100
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = True
        self.remove_nodes = True
        self.add_nodes = True

        self.node_erase_process = NodeEraseProcess(model_part);
        self.h_multiplier = 0.4
        
##        self.Mesher = MSuitePFEMModeler()
        self.Mesher = TriGenPFEMModeler()
#	self.ChooseElement = ChooseElementProcess(model_part, 2)

        self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)

        #detect initial size distribution - note that initially the fluid model part contains
        #all the elements of both structure and fluid ... this is only true after reading the input
        (self.neigh_finder).Execute();
        print "nana"
        self.remeshing_flag = True
        
        self.alpha_shape = 1.2
        self.length_factor = .3
        
        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)
        #U NEED IT FOR ALPHA-shape
        Hfinder  = FindNodalHProcess(model_part);
        Hfinder.Execute();
    
        #helper functions - by Riccardo
        self.PfemUtils = PfemUtils()
        self.MeshMover= MoveMeshProcess(self.model_part);

        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        
        self.solver = NewtonRaphsonOssStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
        (self.solver).SetEchoLevel(self.echo_level)

        (self.neigh_finder).Execute();
        self.Remesh() 
	    
                 
    #######################################################################   
    def Solve(self):
        self.Remesh()
        (self.solver).Solve()
        (self.solver).Clear()

    def EstimateDeltaTime(self,min_dt,max_dt):
        print "Estimating delta time"
        return (self.PfemUtils).EstimateDeltaTime(min_dt,max_dt,self.model_part)

        

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################
    def Remesh(self):

        #and erase bad nodes
        print "BEFOR TOUCHING"
        if (self.remeshing_flag==True):
            print "DI"
            (self.neigh_finder).ClearNeighbours();

            ((self.model_part).Elements).clear();
            ((self.model_part).Conditions).clear();

            (self.PfemUtils).MoveLonelyNodes(self.model_part)
            (self.MeshMover).Execute();
            (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );

 
            (self.node_erase_process).Execute();
            print "AFTER MkkkRK"

        if (self.remeshing_flag==True):
            (self.Mesher).ReGenerateMesh("Fluid2DASGS","Condition2D",self.model_part,self.node_erase_process,self.remove_nodes, self.add_nodes, self.alpha_shape,0.4)
#	    (self.ChooseElement).Execute();   


        (self.neigh_finder).Execute();
        
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)

        for node in self.model_part.Nodes:            
            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)

        ##################################################################
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
        




