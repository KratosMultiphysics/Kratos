from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *

#from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEM2Application import *
from KratosMultiphysics.MeshingApplication import *
#from KratosMultiphysics.PFEMApplication import *
#from KratosMultiphysics.ConvectionDiffusionApplication import *
#from KratosMultiphysics.ULFApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

#from KratosMultiphysics.StructuralApplication import *



def AddVariables(model_part,p_model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);

    model_part.AddNodalSolutionStepVariable(PRESSURE);

    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);

    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)

    model_part.AddNodalSolutionStepVariable(BODY_FORCE);

    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);




    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);



    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(DISTANCE);
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(FORCE_CM)

    model_part.AddNodalSolutionStepVariable(MOMENTUM_CM)

   
    
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    #model_part.AddNodalSolutionStepVariable(TO_ERASE)
    model_part.AddNodalSolutionStepVariable(MATERIAL_VARIABLE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_WATER_ELEMENT)
    model_part.AddNodalSolutionStepVariable(IS_WATER)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    

    
    model_part.AddNodalSolutionStepVariable(IS_VISITED)
    model_part.AddNodalSolutionStepVariable(NODAL_H)

    model_part.AddNodalSolutionStepVariable(IS_POROUS)
    model_part.AddNodalSolutionStepVariable(IS_EROSIONABLE)
    

    model_part.AddNodalSolutionStepVariable(PRESSUREAUX)

    model_part.AddNodalSolutionStepVariable(BULK_MODULUS)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(VELOCITY)






    p_model_part.AddNodalSolutionStepVariable(VELOCITY);
    #p_model_part.AddNodalSolutionStepVariable(FRACT_VEL);

    p_model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);


    p_model_part.AddNodalSolutionStepVariable(PRESSURE);
    p_model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);

    p_model_part.AddNodalSolutionStepVariable(PRESS_PROJ);

    p_model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    p_model_part.AddNodalSolutionStepVariable(NODAL_AREA)

    p_model_part.AddNodalSolutionStepVariable(BODY_FORCE);

    p_model_part.AddNodalSolutionStepVariable(DENSITY);
    p_model_part.AddNodalSolutionStepVariable(VISCOSITY);




    p_model_part.AddNodalSolutionStepVariable(ACCELERATION);
    p_model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);



    p_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    p_model_part.AddNodalSolutionStepVariable(DISTANCE);
    p_model_part.AddNodalSolutionStepVariable(FORCE)
    p_model_part.AddNodalSolutionStepVariable(FORCE_CM)

    p_model_part.AddNodalSolutionStepVariable(MOMENTUM_CM)

   
    
    p_model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    #p_model_part.AddNodalSolutionStepVariable(ERASE_FLAG)
    p_model_part.AddNodalSolutionStepVariable(MATERIAL_VARIABLE);
    p_model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    p_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    p_model_part.AddNodalSolutionStepVariable(IS_FLUID)
    p_model_part.AddNodalSolutionStepVariable(IS_WATER_ELEMENT)
    p_model_part.AddNodalSolutionStepVariable(IS_WATER)
    p_model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    p_model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    

    
    p_model_part.AddNodalSolutionStepVariable(IS_VISITED)
    p_model_part.AddNodalSolutionStepVariable(NODAL_H)

    p_model_part.AddNodalSolutionStepVariable(IS_POROUS)
    p_model_part.AddNodalSolutionStepVariable(IS_EROSIONABLE)
    

    p_model_part.AddNodalSolutionStepVariable(PRESSUREAUX)

    p_model_part.AddNodalSolutionStepVariable(BULK_MODULUS)
    p_model_part.AddNodalSolutionStepVariable(NORMAL)
    p_model_part.AddNodalSolutionStepVariable(VELOCITY)


def AddDofs(model_part,p_model_part):
  
    for node in model_part.Nodes:
        node.AddDof(PRESSURE); 
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);

   

def ReadRestartFile(FileName,nodes):
   NODES = nodes
   aaa = open(FileName)
   for line in aaa:
       exec(line)
       
  

class FracStepSolver:
    
    def __init__(self,model_part,p_model_part,domain_size):

        
        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part)

        self.model_part = model_part
        self.p_model_part = p_model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.vel_toll = 0.0001;
        self.press_toll = 0.001;
        self.max_vel_its = 6;
        self.max_press_its = 5;
        self.time_order = 2;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = False; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False;

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        #self.velocity_linear_solver = BICGSTABSolver(1e-6, 5000, pDiagPrecond)
        
        self.velocity_linear_solver =  SuperLUSolver()

        self.pressure_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)

        self.dynamic_tau = 0.001
        self.activate_tau2 = False


        ##handling slip condition
        self.slip_conditions_initialized = False
        self.neigh_finder = FindNodalNeighboursProcess(model_part)
        self.compute_reactions=False
        self.timer=Timer()    
        #self.particle_utils = ParticleUtils2D()

        #self.normal=NormalCalculationUtils()
        #self.actnormal=NormalToWallCalculationUtils()
        self.Mesher1 = TriGenPFEMModeler()

        #self.PfemUtils = PfemUtils()

        self.node_erase_process = NodeEraseProcess(model_part);

        self.Pfem2_apply_bc_process = Pfem2ApplyBCProcess(model_part);
        self.Pfem2Utils = Pfem2Utils()
        self.mark_outer_nodes_process = MarkOuterNodesProcess(model_part);
        
        if(domain_size == 2):
            self.particle_utils = ParticleUtils2D()	
            self.Mesher = TriGenPFEMModeler()
            
            self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part)
            #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,2, 10)
            self.elem_neighbor_finder = FindElementalNeighboursProcess(model_part, 2, 10)	
            
        elif (domain_size == 3):
            #self.Mesher = TetGenModeler()
            #improved mesher
            self.particle_utils = ParticleUtils3D()
            
            self.Mesher = TetGenPfemModeler()
            self.fluid_neigh_finder = FindNodalNeighboursProcess(model_part)
            #this is needed if we want to also store the conditions a node belongs to
            self.condition_neigh_finder = FindConditionsNeighboursProcess(model_part,3, 20)
            self.elem_neighbor_finder = FindElementalNeighboursProcess(model_part, 20, 30)
     
        self.mark_fluid_process = MarkFluidProcess(model_part);

    def Initialize(self):
        (self.neighbour_search).Execute()

        self.domain_size = int(self.domain_size)
        self.laplacian_form = int(self.laplacian_form)

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.vel_toll = float(self.vel_toll)
        self.press_toll = float(self.press_toll)
        self.max_vel_its = int(self.max_vel_its)
        self.max_press_its = int(self.max_press_its)
        self.time_order = int(self.time_order)
        self.domain_size = int(self.domain_size)


        self.predictor_corrector = bool(self.predictor_corrector)
        self.solver = FracStepStrategy( self.model_part, self.velocity_linear_solver,self.pressure_linear_solver, self.ReformDofAtEachIteration, self.vel_toll, self.press_toll, self.max_vel_its, self.max_press_its, self.time_order, self.domain_size,self.predictor_corrector)


        (self.solver).SetEchoLevel(self.echo_level)

        (self.fluid_neigh_finder).Execute();
        (self.Pfem2_apply_bc_process).Execute();  

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)

		#self.Remesh2();
        self.RemeshAux(); #solo multifluido

        
    def Solve(self):
	

        (self.solver).Solve()

    def solve3(self):
	

        (self.solver).SolveStep3(1.0)

    def solve5(self):
	

        (self.solver).SolveStep5(1.0)

    def solve4(self):
	

        (self.solver).SolveStep4(1.0)

    def solve5(self):
	

        (self.solver).SolveStep5(1.0)


    def Projections(self):
	
        (self.solver).SolveStep5()

    def Reactions(self):


        (self.solver).Compute()


    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True


    def Remesh2(self):
	
        h_factor=1.0;
        alpha_shape=1.2;


        self.node_erase_process = NodeEraseProcess(self.model_part);

        (self.Mesher).ReGenerateMesh("Fluid3DGLS","SurfaceCondition3D3N", self.model_part, self.node_erase_process,False, False, alpha_shape, h_factor)


        (self.fluid_neigh_finder).Execute();
        (self.condition_neigh_finder).Execute();


    def RemeshAux(self):
	
        
        h_factor=0.1 
        alpha_shape=1.4;


        for node in (self.model_part).Nodes:
            node.Set(TO_ERASE, False)
        #dambreak
        for node in (self.model_part).Nodes: 
            node.SetSolutionStepValue(NODAL_H,0,0.0011) 

        self.node_erase_process = NodeEraseProcess(self.model_part);
        box_corner1 = Vector(3); 
        box_corner2 = Vector(3); 

        

        box_corner1[0]=-0.014653
        box_corner1[1]=-0.16#0.31315#-0.1689
        box_corner1[1]=-0.15#0.31315#-0.1689
        box_corner1[2]=-10.0
        box_corner2[0]=0.015702
        box_corner2[1]=0.01189
        box_corner2[2]=10.0




        box_corner1[0]=-10.0
        box_corner1[1]=-10.0#0.31315#-0.1689
        box_corner1[2]=-10.0
        box_corner2[0]=10.0
        box_corner2[1]=10.0
        box_corner2[2]=10.0

        for node in (self.model_part).Nodes: 
            node.SetSolutionStepValue(NODAL_H,0,0.15) 


        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2

        self.Pfem2Utils.MarkLonelyNodesForErasing(self.model_part)

        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);

        self.node_erase_process.Execute()


        if (self.domain_size == 2):
            (self.Mesher).ReGenerateMesh("QFluid2D","LineCondition2D2N", self.model_part, self.node_erase_process, True, False, alpha_shape, h_factor)
        elif (self.domain_size == 3):
            
            (self.Mesher).ReGenerateMesh("QFluid3D","SurfaceCondition3D3N", self.model_part, self.node_erase_process, True, False, alpha_shape, h_factor)            


        (self.fluid_neigh_finder).Execute();
        (self.elem_neighbor_finder).Execute()
        (self.condition_neigh_finder).Execute();

        (self.Pfem2_apply_bc_process).Execute();
        (self.mark_fluid_process).Execute();
        
    def FindNeighbours(self):
        (self.neigh_finder).Execute();
   

    def WriteRestartFile(self,FileName):
        restart_file = open(FileName + ".mdpa",'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintElements("Fluid2DGLS",self.model_part.Elements,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(PRESSURE,"PRESSURE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VISCOSITY,"VISCOSITY",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(DENSITY,"DENSITY",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(TEMPERATURE,"TEMPERATURE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(IS_STRUCTURE,"IS_STRUCTURE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(IS_FREE_SURFACE,"IS_FREE_SURFACE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(IS_FREE_SURFACE,"IS_LAGRANGIAN_INLET",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(IS_BOUNDARY,"IS_BOUNDARY",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(YCH4,"YCH4",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(TEMPERATURE,"TEMPERATURE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(FACE_HEAT_FLUX,"FACE_HEAT_FLUX",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(ENTHALPY,"ENTHALPY",self.model_part.Nodes,restart_file)
        restart_file.close() 



    def CalculateDistanceAndDiviedSet(self,domain_size):
        
        (self.neigh_finder).Execute();
        distance_tools = ElemBasedDistanceUtilities(self.model_part)
        distance_calculator = BodyDistanceCalculationUtils()

        #assign IS_VISITED1 to elem with DISTANCE>=0 and change DSITANCE to posetive for external ones
        ##Assign Zero distance to interface nodes
        for node in (self.model_part).Nodes:
           
            node.SetSolutionStepValue(DISTANCE,0,0.0)
            node.SetValue(IS_VISITED,0.0)
            
        for node in (self.model_part).Nodes:
            if(node.GetSolutionStepValue(IS_INTERFACE)== 1.0):
                
                node.SetSolutionStepValue(DISTANCE,0,0.0)
                node.SetValue(IS_VISITED,1.0)

            else:
                node.SetValue(IS_VISITED,0.0)



        ##        distance_tools.MarkExternalAndMixedNodes()
        ##        distance_tools.ChangeSignToDistance()

        #calculate distances towards the interior of the domain
        if(domain_size == 2):
            distance_calculator.CalculateDistances2D((self.model_part).Elements,DISTANCE, True);
        else:
            distance_calculator.CalculateDistances3D((self.model_part).Elements,DISTANCE, True);

        #change sign 
        distance_tools.ChangeSignToDistance()

        #mark as visited all of the nodes inside the fluid domain
        ####        distance_tools.MarkInternalAndMixedNodes()
        print ((self.model_part).Elements).Size()

        
    def DistToH(self):
         
         for node in (self.model_part).Nodes: 
             node.SetSolutionStepValue(NODAL_H,0,0.05)

    def DistToH1(self):
         
         for node in (self.model_part).Nodes: 
             current_dist = node.GetSolutionStepValue(DISTANCE,0)
             if(abs(current_dist) <= 0.03):
                 node.SetSolutionStepValue(NODAL_H,0,0.02)
             if(0.03<abs(current_dist) and abs(current_dist)<0.09):
                 node_H = 0.03 #+ slope*abs(current_dist)0.005
                 node.SetSolutionStepValue(NODAL_H,0,0.012)
             else:
                 node.SetSolutionStepValue(NODAL_H,0,0.012)

         


         for node in (self.model_part).Nodes:
             if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                 node.SetSolutionStepValue(NODAL_H,0,0.020)

    def DistToH2(self):
         #possible_h = self.CalculateRadius()
         #print possible_h

         
         for node in (self.model_part).Nodes: 
             current_dist = node.GetSolutionStepValue(DISTANCE,0)
             if(abs(current_dist) <= 0.06): #0.15
                 node.SetSolutionStepValue(NODAL_H,0,0.20)
             #if(0.06<abs(current_dist) and abs(current_dist)<0.3):
             #    node_H = 0.03 #+ slope*abs(current_dist)0.005
             #    node.SetSolutionStepValue(NODAL_H,0,0.20)
             else:
                 node.SetSolutionStepValue(NODAL_H,0,0.012)

         


         for node in (self.model_part).Nodes:
             if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                 node.SetSolutionStepValue(NODAL_H,0,0.20)


    def DistToH4(self):
         for node in (self.model_part).Nodes: 
             current_dist = node.GetSolutionStepValue(DISTANCE,0)
             if(abs(current_dist) <= 0.06): #0.15
                 node.SetSolutionStepValue(NODAL_H,0,0.010)
             #if(0.06<abs(current_dist) and abs(current_dist)<0.3):
             #    node_H = 0.03 #+ slope*abs(current_dist)0.005
             #    node.SetSolutionStepVale(NODAL_H,0,0.20)
             else:
                 node.SetSolutionStepValue(NODAL_H,0,0.010)

         


         for node in (self.model_part).Nodes:
             if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                 node.SetSolutionStepValue(NODAL_H,0,0.012)

    def CalculateRadius(self):             
        max_radi = 0.0
        for node in (self.model_part).Nodes:
            if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                X_ref = node.X
                Y_ref = node.Y

        for node in (self.model_part).Nodes:
            if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                radi = pow(node.X-X_ref , 2)+pow(node.Y-Y_ref , 2)
                if(radi>max_radi):
                    max_radi = radi

        max_radi = pow(max_radi,0.5)

        if (max_radi == 0.0):
            max_radi = 0.076
        return max_radi








