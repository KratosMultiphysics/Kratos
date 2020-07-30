#importing the Kratos Library
####THIS IS A SOLVER THAT WAS DERIVED FROM ULF_FSI.. BUT THIS ONE IS ESPECIALLY DEVELOPED FOR THE GLASS SIMULATION
####SPECIAL FEATURES: APPLICATION OF EXTERNAL PRESSURE AT THE INBLOW
#### USE OF THE AXISYMMETRIC FORMULATION (SEE ELEMENT TYPE AND THE PRESSURE CALCULATE PROCESS)
#THE FORMULATION IS THE ELEMENTWISE CONSTANT PRESSURE
from KratosMultiphysics import *
from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.MeshingApplication import *

import time

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_LAGRANGIAN_INLET);
    model_part.AddNodalSolutionStepVariable(BULK_MODULUS);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);

def AddDofs(model_part, compute_reactions):
    if (compute_reactions==0):
        for node in model_part.Nodes:
            #adding dofs
            node.AddDof(DISPLACEMENT_X);
            node.AddDof(DISPLACEMENT_Y);
            node.AddDof(DISPLACEMENT_Z);
    elif (compute_reactions==1):
        for node in model_part.Nodes:
  #adding dofs
              node.AddDof(DISPLACEMENT_X, REACTION_X);
              node.AddDof(DISPLACEMENT_Y, REACTION_Y);
              node.AddDof(DISPLACEMENT_Z, REACTION_Z);

##THIS FUNCTION STORES THE NODES OF THE SECOND MOULD IN A SEPARATE SUBPART (THAT WILL BE DISACTIVATED IN THE STEPS OF THE FIRST BLOW)
##THE NODES OF THE SECOND MOULD ARE DISTINGUISHED BY THE FLAG
#def CreateSubModelPartsFirstSecondBlow(total_model_part, second_mould_flag_value):
  #init_domain_model_part=total_model_part.CreateSubModelPart("InitialDomain");
  #second_mould_model_part=total_model_part.CreateSubModelPart("SecondMould");
  #for node in total_model_part.Nodes:
    #if (node.GetSolutionStepValue(FLAG_VARIABLE)==second_mould_flag_value):
      #second_mould_model_part.AddNode(node, 0 );
    #else:
      #init_domain_model_part.AddNode(node, 0 );

  #for element in total_model_part.Elements:
    #init_domain_model_part.AddElement(element, 0 );
  #for condition in total_model_part.Conditions:
    #init_domain_model_part.AddCondition(condition, 0 );

class ULF_FSISolver:

    def __init__(self, fluid_model_part, structure_model_part, combined_model_part, compute_reactions, box_corner1,box_corner2, domain_size, add_nodes, blow_pressure):
        self.domain_size=domain_size;
        self.compute_reactions=compute_reactions
        self.echo_level = 0
        self.blow_pressure=blow_pressure

        #saving the different model parts
        self.combined_model_part  = combined_model_part; #contains both structure and fluid
        self.fluid_model_part     = fluid_model_part; #contains only fluid elements
        self.structure_model_part = structure_model_part; #contains only structural elements


        #adaptivity options
        self.add_nodes=bool(add_nodes)
        print(("Add nodes? ", self.add_nodes))
        #time integration scheme
        damp_factor = -0.3
        self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(damp_factor)

        #definition of the solvers
#        self.model_linear_solver =  SkylineLUFactorizationSolver()
        pDiagPrecond = DiagonalPreconditioner()
        #pILUPrecond = ILU0Preconditioner()
        #self.model_linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)
        #self.model_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
        self.model_linear_solver =  CGSolver(1e-6, 5000, pDiagPrecond)

        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-4,1e-6)

        #self.pressure_calculate_process = PressureCalculateProcess(fluid_model_part,domain_size);
        self.pressure_calculate_process_axisym = PressureCalculateProcessAxisym(fluid_model_part,domain_size);
        self.mark_free_surface_process = MarkFreeSurfaceProcess(fluid_model_part);
        self.ulf_time_step_dec_process = UlfTimeStepDecProcess(fluid_model_part);
        self.mark_fluid_process = MarkFluidProcess(fluid_model_part);
        self.mark_close_nodes_process = MarkCloseNodesProcess(fluid_model_part);

        self.mark_outer_nodes_process = MarkOuterNodesProcess(fluid_model_part);
        self.node_erase_process = NodeEraseProcess(fluid_model_part);

        #tools to save and merge the structural contributions
        self.save_structure_model_part_process = SaveStructureModelPartProcess();
        self.save_structure_conditions_process = SaveStructureConditionsProcess();

        self.merge_model_parts_process = MergeModelPartsProcess();

        ###temporary ... i need it to calculate the nodal area
        self.UlfUtils = UlfUtils()
        #self.PfemUtils = PfemUtils()

        #self.save_structural_elements
        self.alpha_shape = 1.5;
        self.h_multiplier = 0.1

        ##saving the limits of the box (all the nodes external to this will be erased)
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2

        if(domain_size == 2):
            #self.Mesher = TriGenPFEMModeler()
            self.Mesher = TriGenGLASSModeler()
            self.combined_neigh_finder = FindNodalNeighboursProcess(combined_model_part)
            self.fluid_neigh_finder = FindNodalNeighboursProcess(fluid_model_part)
            #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(fluid_model_part,2, 10)

        elif (domain_size == 3):
            #self.Mesher = TetGenModeler()
            #improved mesher
            self.Mesher = TetGenPfemModeler()
            self.combined_neigh_finder = FindNodalNeighboursProcess(combined_model_part)
            self.fluid_neigh_finder = FindNodalNeighboursProcess(fluid_model_part)
            #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(fluid_model_part,3, 20)


        print("after reading all the model contains:")
        print((self.fluid_model_part))

        #detect initial size distribution - note that initially the fluid model part contains
        #all the elements of both structure and fluid ... this is only true after reading the input
        #Finding the nodal_h for the nodes that are not lonely
        (self.fluid_neigh_finder).Execute();
        Hfinder  = FindNodalHProcess(self.fluid_model_part);
        Hfinder.Execute();
        self.ResetNodalHAtLonelyNodes()
        #assigning average nodal h to lonely nodes
        self.AssignHtoLonelyStructureNodes()
    #######################################################################
    #delta time estimation based on the non-negativity of the jacobian
    def EstimateDeltaTime(self,max_dt,domain_size):
    	#return (self.UlfUtils).EstimateDeltaTime(min_dt,max_dt,self.combined_model_part)
        return (self.ulf_time_step_dec_process).EstimateDeltaTime(max_dt,domain_size)


    #######################################################################
    def Initialize(self):

        #creating the solution strategy
        CalculateReactionFlag = bool(self.compute_reactions)
        ReformDofSetAtEachStep = True
        MoveMeshFlag = True
        import ulf_strategy_PGLASS
        self.solver = ulf_strategy_PGLASS.ULFStrategyPython(self.combined_model_part,self.time_scheme,self.model_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag,self.domain_size)

        print(("self.echo_level = " , self.echo_level))
        (self.solver).SetEchoLevel(self.echo_level)
        print("finished initialization of the fluid strategy")

        #saving the structural elements
        (self.mark_fluid_process).Execute(); #we need this before saving the structrural elements
        print("Saving STRUCTURE")
        (self.save_structure_model_part_process).SaveStructure(self.fluid_model_part, self.structure_model_part);
        (self.save_structure_conditions_process).SaveStructureConditions(self.fluid_model_part, self.structure_model_part);

        #marking the fluid
        (self.fluid_neigh_finder).Execute();
        (self.mark_free_surface_process).Execute();
        (self.mark_fluid_process).Execute();

        #remeshing before the first solution
        #for node in self.fluid_model_part.Nodes:
        #  node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)
        self.Remesh();
        #for node in self.fluid_model_part.Nodes:
        #  if (node.GetSolutionStepValue(IS_FREE_SURFACE)==1.0):
        #    node.SetSolutionStepValue(FLAG_VARIABLE,0,1.0)
        #AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.combined_model_part)


    ######################################################################
    def CheckForInvertedElements(self):
        volume = (self.UlfUtils).CalculateVolume(self.combined_model_part,self.domain_size)
        #volume = (self.UlfUtils).CalculateVolume(self.fluid_model_part,self.domain_size)
        inverted_elements = False
        if(volume < 0.0):
            volume = - volume
            inverted_elements = True

        return [inverted_elements, volume]

    #
    def Solve(self, lagrangian_inlet_process):

        print("solving the fluid problem")
        inverted_elements = (self.solver).Solve(self.domain_size, self.UlfUtils)
        print("successful solution of the fluid ")
        self.lagrangian_inlet_process = lagrangian_inlet_process
        reduction_factor = 0.5
        max_reduction_steps = 5
        time_reduction_step = 0
        while(inverted_elements == True and time_reduction_step <= max_reduction_steps):
            print(" *************************************************** ")
            print("inverted element found ... reducing the time step")
            (self.UlfUtils).ReduceTimeStep(self.combined_model_part,reduction_factor);
            (self.UlfUtils).ReduceTimeStep(self.fluid_model_part,reduction_factor);
            (self.UlfUtils).ReduceTimeStep(self.structure_model_part,reduction_factor);
            print(("reduction_step = ", time_reduction_step))
            time_reduction_step = time_reduction_step + 1

            #copying vars from the old step
##            for node in (self.combined_model_part).Nodes:
##                pold = node.GetSolutionStepValue(PRESSURE,1);
##                dispold = node.GetSolutionStepValue(DISPLACEMENT,1);
##                velold = node.GetSolutionStepValue(VELOCITY,1);
##                accold = node.GetSolutionStepValue(ACCELERATION,1);
##
##                node.SetSolutionStepValue(PRESSURE,0,pold);
##                node.SetSolutionStepValue(DISPLACEMENT,0,dispold);
##                node.SetSolutionStepValue(VELOCITY,0,velold);
##                node.SetSolutionStepValue(ACCELERATION,0,accold);

            self.solver.MoveMesh()
            print("time step reduction completed")
            print(" *************************************************** ")

            (self.solver).Solve(self.domain_size,self.UlfUtils)
            [inverted_elements,vol] = self.CheckForInvertedElements()

        if(inverted_elements == True):

            print("***********************************************************************")
            print("***********************************************************************")
            print("CRITICAL: ... element is still inverted after reducing the time step")
            print("***********************************************************************")
            print("***********************************************************************")
            factor = 2.0**5 #this is the original time step
            (self.UlfUtils).ReduceTimeStep(self.combined_model_part,factor);
            (self.UlfUtils).ReduceTimeStep(self.fluid_model_part,factor);
            (self.UlfUtils).ReduceTimeStep(self.structure_model_part,factor);

##            for node in (self.combined_model_part).Nodes:
##                pold = node.GetSolutionStepValue(PRESSURE,1);
##                dispold = node.GetSolutionStepValue(DISPLACEMENT,1);
##                velold = node.GetSolutionStepValue(VELOCITY,1);
##                accold = node.GetSolutionStepValue(ACCELERATION,1);
##
##                node.SetSolutionStepValue(PRESSURE,0,pold);
##                node.SetSolutionStepValue(DISPLACEMENT,0,dispold);
##                node.SetSolutionStepValue(VELOCITY,0,velold);
##                node.SetSolutionStepValue(ACCELERATION,0,accold);

            self.solver.MoveMesh()

            print("advancing in time without doing anything...")
            (self.solver).PredictionStep(self.domain_size,self.UlfUtils)


        #print "pressure contribution process" - to be executed using exclusively fluid elements
        #and neighbouring relationships
        self.lagrangian_inlet_process.Execute()

        (self.fluid_neigh_finder).Execute();
        (self.UlfUtils).CalculateNodalArea(self.fluid_model_part,self.domain_size);

        #(self.pressure_calculate_process).Execute();
        self.pressure_calculate_process_axisym.Execute();
        #self.lagrangian_inlet_process.Execute()
        #print "remeshing"
        #print "Solution toook time =========================================================",time.time()-step_initialization_start
        self.Remesh();


   ######################################################################
######################################################################
    def Remesh(self):
        self.UlfUtils.MarkNodesTouchingWall(self.fluid_model_part, self.domain_size, 0.1)
        self.UlfUtils.MarkExcessivelyCloseNodes(self.fluid_model_part.Nodes, 0.00005)

        ##erase all conditions and elements prior to remeshing
        ((self.combined_model_part).Elements).clear();
        ((self.combined_model_part).Conditions).clear();
        ((self.combined_model_part).Nodes).clear();
        ((self.fluid_model_part).Elements).clear();
        #((self.fluid_model_part).Conditions).clear();


        #self.UlfUtils.DeleteFreeSurfaceNodesBladder(self.fluid_model_part)
        #marking nodes outside of the bounding box
        #self.UlfUtils.MarkLonelyNodesForErasing(self.fluid_model_part)
        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);
        self.node_erase_process.Execute()


        h_factor=0.25;
        ##remesh CHECK for 3D or 2D
        if (self.domain_size == 2):
            #(self.Mesher).ReGenerateMesh("UpdatedLagrangianFluid2D","LineCondition2D2N", self.fluid_model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
            #(self.Mesher).ReGenerateMesh("UlfAxisym","LineCondition2D2N", self.fluid_model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
            (self.Mesher).ReGenerateMeshGlass("UlfAxisym","LineCondition2D2N", self.fluid_model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)
        elif (self.domain_size == 3):
            (self.Mesher).ReGenerateMesh("UpdatedLagrangianFluid3D","SurfaceCondition3D3N", self.fluid_model_part, self.node_erase_process, True, self.add_nodes, self.alpha_shape, h_factor)

        ##calculating fluid neighbours before applying boundary conditions
        (self.fluid_neigh_finder).Execute();
        (self.condition_neigh_finder).Execute();

        #print "marking fluid" and applying fluid boundary conditions
        (self.mark_free_surface_process).Execute();
        (self.mark_fluid_process).Execute();

        #merging the structural elements back (they are saved in the Initialize)
        (self.merge_model_parts_process).MergeParts(self.fluid_model_part, self.structure_model_part, self.combined_model_part);

        #calculating the neighbours for the overall model
        (self.combined_neigh_finder).Execute();
        (self.UlfUtils).CalculateNodalArea(self.fluid_model_part,self.domain_size);
        #self.UlfUtils.MarkLonelyNodesForErasing(self.fluid_model_part)

        self.node_erase_process.Execute()
        (self.mark_fluid_process).Execute();

        #self.UlfUtils.MarkLonelyNodesForErasing(self.fluid_model_part)#, self.domain_size)
        self.node_erase_process.Execute()
        ##############THIS IS FOR EMBEDDED"""""""""""""""""""""""""

######################################################################################################
        for node in self.fluid_model_part.Nodes:
            #node.SetSolutionStepValue(FLAG_VARIABLE,0,0.0)
            node.SetSolutionStepValue(EXTERNAL_PRESSURE,0,0.0)
#
        #for node in self.fluid_model_part.Nodes:
        #    if (node.GetSolutionStepValue(IS_FREE_SURFACE)==1.0 and node.GetSolutionStepValue(IS_INTERFACE)==0.0 and node.X>0.0001):
        #        node.SetSolutionStepValue(FLAG_VARIABLE, 0, 1.0)

        #time_curr=self.fluid_model_part.ProcessInfo[TIME]
        #print (time_curr)
        #h= 0.1+0.1*time_curr
        #for node in self.fluid_model_part.Nodes:
        #    if (node.Y>h):
        #      node.SetSolutionStepValue(FLAG_VARIABLE,0,0)
        #
        for node in self.fluid_model_part.Nodes:
            if (node.GetSolutionStepValue(FLAG_VARIABLE)==1.0 and node.GetSolutionStepValue(IS_FREE_SURFACE)==1.0):
                node.SetSolutionStepValue(EXTERNAL_PRESSURE,0, self.blow_pressure)
        #FOR PFEM
        NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part.Conditions, self.domain_size)

        for node in self.fluid_model_part.Nodes:
          if (node.GetSolutionStepValue(FLAG_VARIABLE)!=1.0):
            node.SetSolutionStepValue(NORMAL_X,0,0.0)
            node.SetSolutionStepValue(NORMAL_Y,0,0.0)
            node.SetSolutionStepValue(NORMAL_Z,0,0.0)
        AssignPointNeumannConditions().AssignPointNeumannConditionsDisp(self.combined_model_part)
######################################################################################################
        #FOR MEMBRANE
        #print "00000"
        #NormalCalculationUtils().CalculateOnSimplex(self.combined_model_part, self.domain_size)
        #print "aaaa"
        #print self.fluid_model_part
        #AssignPointNeumannConditions().AssignPointNeumannConditions3D(self.combined_model_part)
        #print self.fluid_model_part
        #print "bbbb"


        ##############THIS IS FOR EMBEDDED"""""""""""""""""""""""""
        print("end of remesh fucntion")
    ######################################################################
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
    ######################################################################
    #this function is included since at the findNodalHProcess of Kratos core assigns max value (1.7*e^309) to all nodes
    #and then computes the correct value only for the nodes of elements. Thus lonely nodes remain with this enormous values
    def ResetNodalHAtLonelyNodes(self):
       for node in self.fluid_model_part.Nodes:
          if (node.GetSolutionStepValue(NODAL_H)>100000000.0):
            node.SetSolutionStepValue(NODAL_H,0,0.0)
    ######################################################################
    def AssignHtoLonelyStructureNodes(self):
        nnodes=0
        nodal_h=0.0
        av_nodal_h=0.0
        for node in self.fluid_model_part.Nodes:
            if (node.GetSolutionStepValue(NODAL_H)!=0.0):
               nnodes=nnodes+1;
               nodal_h=nodal_h+node.GetSolutionStepValue(NODAL_H);

        av_nodal_h=nodal_h/nnodes
        for node in self.fluid_model_part.Nodes:
           if (node.GetSolutionStepValue(IS_STRUCTURE)==1.0 and node.GetSolutionStepValue(IS_FLUID)==0.0):
              node.SetSolutionStepValue(NODAL_H, av_nodal_h)
    ######################################################################

