#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.PFEMApplication import PfemUtils
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MeshingApplication import *
# Check that KratosMultiphysics was imported in the main script
#CheckForPreviousImport()

import time

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);    
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_LAGRANGIAN_INLET);    
    model_part.AddNodalSolutionStepVariable(SYMMETRY_CUT);
    model_part.AddNodalSolutionStepVariable(CENTER_LINE);
    model_part.AddNodalSolutionStepVariable(FIXED_WALL);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(BULK_MODULUS);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(IS_WATER);
    model_part.AddNodalSolutionStepVariable(DENSITY_WATER);    
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(CONDUCTIVITY);
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(EMISSIVITY);

def AddDofs(model_part, compute_reactions):
    for node in model_part.Nodes:
      node.AddDof(DISPLACEMENT_X);
      node.AddDof(DISPLACEMENT_Y);
      node.AddDof(DISPLACEMENT_Z);    
      node.AddDof(PRESSURE); 
      node.AddDof(TEMPERATURE); 
     
class ULF_FSISolver:

    def __init__(self, model_part, box_corner1,box_corner2, domain_size, add_nodes, bulk_modulus, density):
        
	#THIS SOLVER IS MEANT FOR FLUID-ONLY PROBLEMS
	self.FSI=0	
        self.domain_size=domain_size;
        self.echo_level = 0
        self.counter = int(0)

        # TO REMESH OR NOT: 0 - no remeshing, 1 - remeshing
        self.add_nodes = bool(add_nodes)                
        self.bulk_modulus = bulk_modulus
        self.density = density
        
        #saving the different model parts     
        self.model_part     = model_part


        #time integration scheme (dam_factor= alpha_bossak)
        self.damp_factor = -0.3
        self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.damp_factor)        
        
        #self.conv_criteria = DisplacementCriteria(1e-6,1e-9)
        self.conv_criteria = IncrementalDisplacementCriteria(1e-3,1e-6)

        #self.pressure_calculate_process = PressureCalculateProcess(fluid_model_part,domain_size);
        self.ulf_apply_bc_process = UlfApplyBCProcess(model_part);
        self.ulf_time_step_dec_process = UlfTimeStepDecProcess(model_part);
        self.mark_fluid_process = MarkFluidProcess(model_part);
        self.mark_close_nodes_process = MarkCloseNodesProcess(model_part);
        self.mark_outer_nodes_process = MarkOuterNodesProcess(model_part);
        self.node_erase_process = NodeEraseProcess(model_part);
        
        
        ###temporary ... i need it to calculate the nodal area
        self.UlfUtils = UlfUtils()
        self.PfemUtils = PfemUtils()

        #self.save_structural_elements
        self.alpha_shape = 1.5;
        self.h_multiplier = 0.1

        ##saving the limits of the box (all the nodes external to this will be erased)
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2
        
        if(domain_size == 2):
            self.Mesher = TriGenModeler()
            #self.Mesher = TriGenPFEMModeler()            
            self.neigh_finder = FindNodalNeighboursProcess(model_part,9,18)            
             #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,2,10)
        elif (domain_size == 3):
            #improved mesher
            #self.Mesher =TetGenCDT()
            self.Mesher = TetGenPfemModeler()
            #self.Mesher = TetGenModeler()
            self.neigh_finder = FindNodalNeighboursProcess(model_part,20,30)           
            #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,3, 20)
     

        print "after reading all the model contains:"
        print self.model_part

        #detect initial size distribution - note that initially the fluid model part contains
        #all the elements of both structure and fluid ... this is only true after reading the input
        (self.neigh_finder).Execute();
        self.Hfinder  = FindNodalHProcess(model_part);
        self.Hfinder.Execute();   

       
    #######################################################################
    #delta time estimation based on the non-negativity of the jacobian
    def EstimateDeltaTime(self,max_dt,domain_size):
    	#return (self.UlfUtils).EstimateDeltaTime(min_dt,max_dt,self.combined_model_part)
        return (self.ulf_time_step_dec_process).EstimateDeltaTime(max_dt,domain_size)    
    
    #######################################################################
    def Initialize(self):	
        #creating the solution strategy        
        ReformDofSetAtEachStep = True
        MoveMeshFlag = True
        
        import ulf_CDT_strategy
        
        self.solver = ulf_CDT_strategy.ULFFracStrategyPython(self.model_part, self.time_scheme, self.conv_criteria, ReformDofSetAtEachStep,MoveMeshFlag,self.domain_size, self.bulk_modulus, self.density)
        
        print "self.echo_level = " , self.echo_level
        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
        
        #saving the structural elements
        (self.mark_fluid_process).Execute(); #we need this before saving the structrural elements       

        #marking the fluid
        (self.neigh_finder).Execute();
                
        (self.ulf_apply_bc_process).Execute();  
        (self.mark_fluid_process).Execute();
        #caluclating nodal area in order to calculate pressures
        (self.UlfUtils).CalculateNodalArea(self.model_part,self.domain_size);
        #remeshing before the first solution
        self.Remesh();    
        

    ######################################################################
    def CheckForInvertedElements(self):
        volume = (self.UlfUtils).CalculateVolume(self.model_part,self.domain_size)
        inverted_elements = False
        if(volume < 0.0):
            volume = - volume
            inverted_elements = True
        return [inverted_elements,volume]
                         
    #######################################################################
    def Solve(self):
        print "solving the fluid problem"
        inverted_elements = (self.solver).Solve(self.domain_size,self.UlfUtils, self.FSI)
        print "succesful solution of the fluid "

        reduction_factor = 0.5
        max_reduction_steps = 5
        time_reduction_step = 0
        while(inverted_elements == True and time_reduction_step <= max_reduction_steps):
            print " *************************************************** "
            print "inverted element found ... reducing the time step"            
            (self.UlfUtils).ReduceTimeStep(self.model_part,reduction_factor);            
            print "reduction_step = ", time_reduction_step
            time_reduction_step = time_reduction_step + 1


            self.solver.MoveMesh()
            print "time step reduction completed"
            print " *************************************************** "

            (self.solver).Solve(self.domain_size,self.UlfUtils, self.FSI)
            [inverted_elements,vol] = self.CheckForInvertedElements()            

        if(inverted_elements == True):
            
            print "***********************************************************************"
            print "***********************************************************************"
            print "CRITICAL: ... element is still inverted after reducing the time step"
            print "***********************************************************************"
            print "***********************************************************************"
            factor = 2.0**5 #this is the original time step
            (self.UlfUtils).ReduceTimeStep(self.model_part,factor);
            
            self.solver.MoveMesh()

            print "advancing in time without doing anything..."
            (self.solver).PredictionStep(self.domain_size,self.UlfUtils)         
          
        (self.neigh_finder).Execute();       
        self.Remesh();


   ######################################################################
   #  the parameter in this function is set to 1 IN CASE WE DO NOT WANT TO CREATE THE NEW MESH, BUT JUST
   #  the operations on model part
   #  This is done to make switching off/on of remeshing easier
    def Remesh(self):    	
	timeRemesh=time.time()
        ##preventing the nodes from coming tooo close to wall
        self.PfemUtils.MarkNodesTouchingWall(self.model_part, self.domain_size, 0.05)
        #self.PfemUtils.MarkNodesTouchingWall(self.model_part, self.domain_size, 0.05)
	#self.PfemUtils.MarkNodesTouchingWall(self.model_part, self.domain_size, 0.15)
	self.PfemUtils.MarkExcessivelyCloseNodes(self.model_part.Nodes, 0.0001)
	
	#((self.model_part).Elements).clear();
        #((self.model_part).Conditions).clear();
        
        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);
        #adaptivity=True
        
        h_factor=0.25
        if (self.domain_size == 2):	  
	  (self.Mesher).ReGenerateMesh("UlfFrac2D","Condition2D", self.model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
	elif (self.domain_size == 3):
	  #(self.Mesher).GenerateCDT(self.model_part, "UlfFrac3D", False, 0, self.node_erase_process, h_factor)
	  (self.Mesher).ReGenerateMesh("UlfFrac3D","Condition3D", self.model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
       
	print ("REMESHING COMPLETED!")
	##calculating fluid neighbours before applying boundary conditions
        (self.neigh_finder).Execute();
        #(self.condition_neigh_finder).Execute();

        #print "marking fluid" and applying fluid boundary conditions
        (self.ulf_apply_bc_process).Execute();	    
        (self.mark_fluid_process).Execute();
        
        print self.model_part
        #calculating the neighbours for the overall model        
        #(self.fluid_only_neigh_finder).Execute();
        (self.UlfUtils).CalculateNodalArea(self.model_part,self.domain_size);
        
        self.UlfUtils.MarkLonelyNodesForErasing(self.model_part)
        (self.mark_fluid_process).Execute();
        self.node_erase_process.Execute()        
        timeRemeshEnd=time.time()
        print "Remeshing took ", str(timeRemeshEnd-timeRemesh)
        print "end of remesh function"
	#BodyNormalCalculationUtils().CalculateBodyNormals(self.combined_model_part, self.domain_size)
	#NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.domain_size, FLAG_VARIABLE) 
	
	NormalCalculationUtils().CalculateOnSimplex(self.model_part.Conditions, self.domain_size) 
	#print "aaaa"
	AssignPointNeumannConditions().AssignPointNeumannConditionsDisp(self.model_part) 
	#print "bbbbb"
	#print "COMPUTED NORMALS"
	##to treat the sides in the correct way
	

	#AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.model_part) 
	##########################################################################################################3
##	for node in self.model_part.Nodes:	  
##	  #if (node.Z==0.0 and node.X>0.099 and node.GetSolutionStepValue(FIXED_WALL)==0):
##	  if (node.Z>=-0.0000000000001 and node.X>=-0.198 and node.GetSolutionStepValue(FIXED_WALL)==0):
##	    node.Fix(DISPLACEMENT_X)
##	    node.Fix(DISPLACEMENT_Z)
##	    #node.Fix(DISPLACEMENT_Y)
##	    node.Free(DISPLACEMENT_Y)
##	    node.SetSolutionStepValue(SYMMETRY_CUT, 0, 10.0)
##	    print ("Symmetry CUT 100000000000000000000000000000000000000000000")
##	    print ("Symmetry CUT 100000000000000000000000000000000000000000000")
##	    print ("Symmetry CUT 100000000000000000000000000000000000000000000")
##	    print ("Symmetry CUT 100000000000000000000000000000000000000000000")
##	    print ("Symmetry CUT 100000000000000000000000000000000000000000000")
##	    print ("Symmetry CUT 100000000000000000000000000000000000000000000")
##	    #a=node.GetSolutionStepValue(NORMAL,0)
##	    #b=2.0*a;
##	    #node.SetSolutionStepValue(NORMAL,0, b)
##	    #print "YYYYYYYYYYYYYYYAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"	    
##       
    ######################################################################
    
    ######################################################################
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
        
