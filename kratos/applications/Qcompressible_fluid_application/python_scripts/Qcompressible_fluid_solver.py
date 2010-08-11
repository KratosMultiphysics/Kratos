#importing the Kratos Library
from Kratos import *
from KratosQcompressibleFluidApplication import *
from KratosConvectionDiffusionSpeciesApplication import *
from KratosMeshingApplication import *
from KratosULFApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(DESP);
    model_part.AddNodalSolutionStepVariable(VELOCITY);  
    model_part.AddNodalSolutionStepVariable(ERASE_FLAG);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(VEL_i);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESSUREAUX_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(MASSQ);
    model_part.AddNodalSolutionStepVariable(NODAL_PRESS);
    model_part.AddNodalSolutionStepVariable(NODAL_MASSAUX);
    model_part.AddNodalSolutionStepVariable(NODAL_PRESSAUX);
    model_part.AddNodalSolutionStepVariable(NODAL_MAUX);
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(NODAL_DENSITYAUX);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(BULK_MODULUS);
    model_part.AddNodalSolutionStepVariable(TC);
    model_part.AddNodalSolutionStepVariable(IS_BURN);
    model_part.AddNodalSolutionStepVariable(IS_DRIPPING);
    model_part.AddNodalSolutionStepVariable(MATERIAL_VARIABLE);
    model_part.AddNodalSolutionStepVariable(IS_WALL);
##flag variables
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_PERMANENT);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS);
    model_part.AddNodalSolutionStepVariable(ARRHENIUSAUX);

    model_part.AddNodalSolutionStepVariable(PRESSUREAUX);

    print "variables for the incompressible fluid solver added correctly"

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);

    print "dofs for the incompressible fluid solver added correctly"
    

def ReadRestartFile(FileName,nodes):
   NODES = nodes
   aaa = open(FileName)
   for line in aaa:
       exec(line)
       
 

class QcompressibleFluidSolver:
    
    def __init__(self,model_part,domain_size, gid_io):

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
	self.gid_io = gid_io
	self.counter=0
        ##self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
	print "LALALA"
	print domain_size
        self.domain_size = domain_size
	
        #assignation of parameters to be used
        self.vel_toll = 0.001;
        self.press_toll = 0.001;
        self.max_vel_its = 4;
        self.max_press_its = 3;
        self.time_order = 2;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = True; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False;

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        pILUPrecond = ILU0Preconditioner()
##        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)
        self.velocity_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)##pDiagPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)


        print"MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM"
        print (model_part.Elements).Size()
        print"KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
        
###############MESH CHANGES
        self.UlfUtils = UlfUtils()
        self.qUtils = qUtils()
        self.ulf_time_step_dec_process = UlfTimeStepDecProcess(model_part);
        self.mark_close_nodes_process = MarkCloseNodesProcess(model_part);
        self.node_erase_process = NodeEraseProcess(model_part);
        self.h_multiplier = 0.8
###############MESH CHANGES

        if(domain_size == 2):
            self.Mesher = TriGenPFEMModeler() 
            self.node_erase_process = NodeEraseProcess(model_part);
            self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)
	elif (domain_size == 3):
	    self.Mesher =TetGenPfemModeler() 
            self.neigh_finder = FindNodalNeighboursProcess(model_part,20,30)

        
        (self.neigh_finder).Execute();
#########################
	
        self.remeshing_flag = True
        self.alpha_shape = 1.4



#############MESH CHANGES

        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)
        #U NEED IT FOR ALPHA-shape
        Hfinder  = FindNodalHProcess(model_part);
	 
       	Hfinder.Execute();

        for node in self.model_part.Nodes:
            temp=node.GetSolutionStepValue(NODAL_H)
            node.SetSolutionStepValue(NODAL_H, 1,temp )
            node.SetSolutionStepValue(NODAL_H, 2,temp )
    	

#########################

    def Initialize(self):
       
        
        print "in python: okkio using Coupled Strategy"
        self.solver = ResidualBasedFluidStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)   
        
        (self.solver).SetEchoLevel(self.echo_level)
 	(self.neigh_finder).Execute()
        (self.qUtils).IdentifyFluidNodes(self.model_part);
        (self.qUtils).IdentifyInterfaceNodes(self.model_part);
#        (self.qUtils).IdentifyWallNodes(self.model_part);
        self.Remesh() #comentadoooooooooooooo
        print "finished initialization of the fluid strategy"
        
   
    def Solve(self):

        if(self.ReformDofAtEachIteration == True):
            (self.neigh_finder).Execute()

        #(self.qUtils).IdentifyFluidNodes(self.model_part);

        (self.solver).Solve()
  
    def Pressure(self):

	(self.solver).SolveStep4();

    def Move(self):
	(self.qUtils).QuasiLagrangianMove(self.model_part)#//lo


    def WriteRestartFile(self,FileName):
        backupfile = open(FileName+".py",'w')
        
        import restart_utilities
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(PRESSURE,"PRESSURE",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(DENSITY,"DENSITY",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(VISCOSITY,"VISCOSITY",self.model_part.Nodes,backupfile)

        restart_utilities.PrintRestartFixity_PyFormat(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(PRESSURE,"PRESSURE",self.model_part.Nodes,backupfile)
        
        backupfile.close()



    def Remesh(self):
	print "AAlll"

        (self.solver).Clear();

	print "AAAA"

	(self.neigh_finder).ClearNeighbours();
        ((self.model_part).Elements).clear();
        ((self.model_part).Conditions).clear();            
	print "PORAQUI??????"
#     	for node in self.model_part.Nodes:
            
            #node.SetSolutionStepValue(NODAL_H, 1,0.0048 )
            #node.SetSolutionStepValue(NODAL_H, 2,0.0048 )
        if(self.domain_size == 2):
            (self.Mesher).ReGenerateMesh("QFluid2D","Condition2D", self.model_part, self.node_erase_process,False,False,self.alpha_shape, 0.5)################0.5
	elif (self.domain_size == 3):
            (self.Mesher).ReGenerateMeshElementsQcomp(self.model_part, self.alpha_shape)#self.alpha_shape)

        (self.solver).Clear();

	print "AAAA"

        (self.neigh_finder).Execute();

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)

        for node in self.model_part.Nodes:            
            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)



#############################CHANGE
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
#############################CHANGE


        
            
        print "end of remesh fucntion"
        print self.model_part
        print self.model_part



      ######################################################################
    def CheckForInvertedElements(self):
        volume = (self.qUtils).CalculateVolume(self.model_part,self.domain_size)
        inverted_elements = False
        if(volume < 0.0):
            volume = - volume
            inverted_elements = True
        return [inverted_elements,volume]
        
    ######################################################################   
        
    def Remesh2(self):

        #if(self.ReformDofAtEachIteration == True):
        #    (self.neigh_finder).Execute()

        #print "just before solve"
        #print "probelama"

        #(self.qUtils).IdentifyFluidNodes(self.model_part);
	
	#print "AAAA"

	(self.neigh_finder).ClearNeighbours();
        ((self.model_part).Elements).clear();
        ((self.model_part).Conditions).clear();            
	#print "PORAQUI??????"
     
        if(self.domain_size == 2):
            (self.Mesher).ReGenerateMesh("QFluid2D","Condition2D", self.model_part, self.node_erase_process,True,False,1.6, 0.5)################0.5
	elif (self.domain_size == 3):
            (self.Mesher).ReGenerateMeshElementsQcomp(self.model_part, self.alpha_shape)#self.alpha_shape)

        #(self.qUtils).IdentifyFluidNodes(self.model_part);

#        (self.neigh_finder).Execute();

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)

        for node in self.model_part.Nodes:            
            if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
                node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)

