#importing the Kratos Library
from Kratos import *
from KratosULFApplication import *
from KratosPFEMApplication import PfemUtils
from KratosStructuralApplication import *
from KratosMeshingApplication import *
#import time

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);

    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    #model_part.AddNodalSolutionStepVariable(DISP_FRAC);
    model_part.AddNodalSolutionStepVariable(VAUX);
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
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(BULK_MODULUS);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(NORMAL);

def AddDofs(model_part, compute_reactions):
    if (compute_reactions==0):
	for node in model_part.Nodes:
	      #adding dofs
	      node.AddDof(DISPLACEMENT_X);
	      node.AddDof(DISPLACEMENT_Y);
	      node.AddDof(DISPLACEMENT_Z);
	      node.AddDof(PRESSURE); 
	
    elif (compute_reactions==1):
	for node in model_part.Nodes:
	      #adding dofs
	      node.AddDof(DISPLACEMENT_X, REACTION_X);
	      node.AddDof(DISPLACEMENT_Y, REACTION_Y);
	      node.AddDof(DISPLACEMENT_Z, REACTION_Z);
	      node.AddDof(PRESSURE); 
        

class ULF_FSISolver:

    def __init__(self, out_file, fluid_only_model_part, fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1,box_corner2, domain_size, add_nodes, bulk_modulus, density):
        self.out_file=out_file

	self.compute_reactions=compute_reactions
        self.domain_size=domain_size;
        self.echo_level = 0
        self.counter = int(0)

        # TO REMESH OR NOT: 0 - no remeshing, 1 - remeshing
        self.add_nodes = bool(add_nodes)
        # K - the bulk modulus
        self.bulk_modulus = bulk_modulus
        self.density = density
        
        #saving the different model parts
        self.combined_model_part  = combined_model_part; #contains both structure and fluid
        self.fluid_model_part     = fluid_model_part; #contains only fluid elements, but all nodes!
        self.structure_model_part = structure_model_part; #contains only structural elements

        self.fluid_only_model_part     = fluid_only_model_part; #contains only fluid elements and nodes

        #time integration scheme (dam_factor= alpha_bossak)
        self.damp_factor = -0.3
        self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.damp_factor)

        #definition of the solvers
        self.pres_time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.pres_linear_solver =  SkylineLUFactorizationSolver()
                
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)

        #self.pressure_calculate_process = PressureCalculateProcess(fluid_model_part,domain_size);
        self.ulf_apply_bc_process = UlfApplyBCProcess(fluid_model_part);
        self.ulf_time_step_dec_process = UlfTimeStepDecProcess(fluid_model_part);
        self.mark_fluid_process = MarkFluidProcess(fluid_model_part);
        self.mark_close_nodes_process = MarkCloseNodesProcess(fluid_model_part);
        self.mark_outer_nodes_process = MarkOuterNodesProcess(fluid_model_part);
        self.node_erase_process = NodeEraseProcess(fluid_model_part);

        #tools to save and merge the structural contributions
        self.save_structure_model_part_process = SaveStructureModelPartProcess();
        self.save_structure_conditions_process = SaveStructureConditionsProcess();
        self.merge_model_parts_process = MergeModelPartsProcess();
        #this will save fluid only model part
        self.save_fluid_only_process=SaveFluidOnlyProcess();

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
            #self.Mesher = TriGenModeler()
            self.Mesher = TriGenPFEMModeler()
            self.combined_neigh_finder = FindNodalNeighboursProcess(combined_model_part,9,18)
            self.fluid_neigh_finder = FindNodalNeighboursProcess(fluid_model_part,9,18)
            self.fluid_only_neigh_finder = FindNodalNeighboursProcess(fluid_only_model_part,9,18)
             #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(fluid_model_part,2,10)
        elif (domain_size == 3):
            #improved mesher
            self.Mesher = TetGenPfemModeler()
            #self.Mesher = TetGenModeler()
            self.combined_neigh_finder = FindNodalNeighboursProcess(combined_model_part,20,30)
            self.fluid_neigh_finder = FindNodalNeighboursProcess(fluid_model_part,20,30)
            self.fluid_only_neigh_finder = FindNodalNeighboursProcess(fluid_only_model_part,20,30)
            #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(fluid_model_part,3, 20)
     

        print "after reading all the model contains:"
        print self.fluid_model_part

        #detect initial size distribution - note that initially the fluid model part contains
        #all the elements of both structure and fluid ... this is only true after reading the input
        (self.fluid_neigh_finder).Execute();
        Hfinder  = FindNodalHProcess(fluid_model_part);
        Hfinder.Execute();
    

       
    #######################################################################
    #delta time estimation based on the non-negativity of the jacobian
    def EstimateDeltaTime(self,max_dt,domain_size):
    	#return (self.UlfUtils).EstimateDeltaTime(min_dt,max_dt,self.combined_model_part)
        return (self.ulf_time_step_dec_process).EstimateDeltaTime(max_dt,domain_size)
    
    #######################################################################
    #this function is needed only in case there is a lagrangian nodes-inlet in the problem
    #three numbers - are veolicties in x y and z directions
    #def MoveInletNodes(self, model_part):
    #    (self.UlfUtils).MoveInletNodes(self.fluid_model_part, 0.1, 0.0, 0.0)
        
    
    #######################################################################
    def Initialize(self):

        #creating the solution strategy
        CalculateReactionFlag = bool(self.compute_reactions)
        ReformDofSetAtEachStep = True
        MoveMeshFlag = True
        
        import ulf_frac_strategy
        
        self.solver = ulf_frac_strategy.ULFFracStrategyPython(self.fluid_only_model_part, self.combined_model_part, self.fluid_model_part, self.time_scheme, self.pres_time_scheme, self.pres_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag,self.domain_size, self.bulk_modulus, self.density)
        
        print "self.echo_level = " , self.echo_level
        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
        
        #saving the structural elements
        (self.mark_fluid_process).Execute(); #we need this before saving the structrural elements

        #we specify domain size, to deal with problems involving membarnes in 3D in a specific way (see save_structure_model_part_process.h
        (self.save_structure_model_part_process).SaveStructure(self.fluid_model_part, self.structure_model_part);
        print "STRUCTURE PART!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!", self.structure_model_part
        #(self.save_structure_conditions_process).SaveStructureConditions(self.fluid_model_part, self.structure_model_part, self.domain_size);

        #creating initially empty container for lagrangian inlet-nodes

        #marking the fluid
        (self.fluid_neigh_finder).Execute();
        (self.combined_neigh_finder).Execute();
        
        (self.ulf_apply_bc_process).Execute();  
        (self.mark_fluid_process).Execute();
        #caluclating nodal area in order to calculate pressures
        (self.UlfUtils).CalculateNodalArea(self.fluid_model_part,self.domain_size);
        #remeshing before the first solution

        self.Remesh();
               

    ######################################################################
    def CheckForInvertedElements(self):
        volume = (self.UlfUtils).CalculateVolume(self.combined_model_part,self.domain_size)
        inverted_elements = False
        if(volume < 0.0):
            volume = - volume
            inverted_elements = True
        return [inverted_elements,volume]
                         
    #######################################################################
    def Solve(self, lagrangian_inlet_process):
      
	self.lagrangian_inlet_process=lagrangian_inlet_process

        print "solving the fluid problem"
        inverted_elements = (self.solver).Solve(self.domain_size,self.UlfUtils)
        print "succesful solution of the fluid "

        reduction_factor = 0.5
        max_reduction_steps = 5
        time_reduction_step = 0
        while(inverted_elements == True and time_reduction_step <= max_reduction_steps):
            print " *************************************************** "
            print "inverted element found ... reducing the time step"
            (self.UlfUtils).ReduceTimeStep(self.combined_model_part,reduction_factor);
            (self.UlfUtils).ReduceTimeStep(self.fluid_model_part,reduction_factor);
            (self.UlfUtils).ReduceTimeStep(self.structure_model_part,reduction_factor);
            print "reduction_step = ", time_reduction_step
            time_reduction_step = time_reduction_step + 1


            self.solver.MoveMesh()
            print "time step reduction completed"
            print " *************************************************** "

            (self.solver).Solve(self.domain_size,self.UlfUtils)
            [inverted_elements,vol] = self.CheckForInvertedElements()            

        if(inverted_elements == True):
            
            print "***********************************************************************"
            print "***********************************************************************"
            print "CRITICAL: ... element is still inverted after reducing the time step"
            print "***********************************************************************"
            print "***********************************************************************"
            factor = 2.0**5 #this is the original time step
            (self.UlfUtils).ReduceTimeStep(self.combined_model_part,factor);
            (self.UlfUtils).ReduceTimeStep(self.fluid_model_part,factor);
            (self.UlfUtils).ReduceTimeStep(self.structure_model_part,factor);
            
            self.solver.MoveMesh()

            print "advancing in time without doing anything..."
            (self.solver).PredictionStep(self.domain_size,self.UlfUtils)     
                
       
        self.lagrangian_inlet_process.Execute()       
          
        (self.fluid_neigh_finder).Execute();       
        self.Remesh();


   ######################################################################
   #  the parameter in this function is set to 1 IN CASE WE DO NOT WANT TO CREATE THE NEW MESH, BUT JUST
   #  the operations on model part
   #  This is done to make switching off/on of remeshing easier
    def Remesh(self):                 
        ##erase all conditions and elements prior to remeshing
        #self.PfemUtils.MarkNodesTouchingWall(self.fluid_model_part, self.domain_size, 0.05)
        ((self.combined_model_part).Elements).clear();
        ((self.combined_model_part).Conditions).clear();
        ((self.combined_model_part).Nodes).clear();
        ((self.fluid_model_part).Elements).clear();
        ((self.fluid_model_part).Conditions).clear();
            
	#self.UlfUtils.MarkNodesCloseToWall(self.fluid_model_part, self.domain_size, 2.5000)
        #mark outer nodes for erasing
        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);
        #adaptivity=True
        #time=self.combined_model_part.ProcessInfo.GetValue(TIME);
        #print "TIMEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE", time
        
        h_factor=0.2
        if (self.domain_size == 2):
	  (self.Mesher).ReGenerateMesh("UlfFrac2D","Condition2D", self.fluid_model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
	elif (self.domain_size == 3):
	  (self.Mesher).ReGenerateMesh("UlfFrac3D","Condition3D", self.fluid_model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
       

        ##calculating fluid neighbours before applying boundary conditions
        (self.fluid_neigh_finder).Execute();
        (self.condition_neigh_finder).Execute();

        #print "marking fluid" and applying fluid boundary conditions
        (self.ulf_apply_bc_process).Execute();
        (self.mark_fluid_process).Execute();

        #merging the structural elements back (they are saved in the Initialize)
        (self.merge_model_parts_process).MergeParts(self.fluid_model_part, self.structure_model_part, self.combined_model_part);
        #this one is used to perfrom some operations exclusively on the fluid elements/nodes
        (self.save_fluid_only_process).SaveFluidOnly(self.fluid_model_part, self.fluid_only_model_part);

        print self.fluid_model_part
        print self.fluid_only_model_part
        print self.combined_model_part
        #calculating the neighbours for the overall model
        (self.combined_neigh_finder).Execute();
        #(self.fluid_only_neigh_finder).Execute();
        (self.UlfUtils).CalculateNodalArea(self.fluid_model_part,self.domain_size);
        
        self.UlfUtils.MarkLonelyNodesForErasing(self.fluid_model_part)
        self.node_erase_process.Execute()
        (self.mark_fluid_process).Execute();
        print "end of remesh function"

       
    ######################################################################
        
               
    ######################################################################
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
        
