# -*- coding: utf-8 -*-
#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *
from KratosMeshingApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(FRACT_VEL);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ);
    model_part.AddNodalSolutionStepVariable(CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_MASS);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS);

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
       
  

class FracStepSolver:
    
    def __init__(self,model_part,domain_size):

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.vel_toll = 0.001;
        self.press_toll = 0.001;
        self.max_vel_its = 6;
        self.max_press_its = 3;
        self.time_order = 2;
        self.CalculateReactions = False;
        self.ReformDofAtEachIteration = False; 
        self.CalculateNormDxFlag = True;
        self.laplacian_form = 2; #1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False;

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        #self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)

        self.dynamic_tau = 0.001
        self.activate_tau2 = False


        ##handling slip condition
        self.slip_conditions_initialized = False

        self.compute_reactions=False
        self.timer=Timer()    
        
        


    def Initialize(self):
        (self.neighbour_search).Execute()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        self.model_part.ProcessInfo.SetValue(ACTIVATE_TAU2, self.activate_tau2);

        self.domain_size = int(self.domain_size)
        self.laplacian_form = int(self.laplacian_form)
        solver_configuration = FractionalStepConfiguration(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.domain_size,self.laplacian_form )

        self.ReformDofAtEachIteration = bool(self.ReformDofAtEachIteration)
        self.vel_toll = float(self.vel_toll)
        self.press_toll = float(self.press_toll)
        self.max_vel_its = int(self.max_vel_its)
        self.max_press_its = int(self.max_press_its)
        self.time_order = int(self.time_order)
        self.domain_size = int(self.domain_size)
        self.predictor_corrector = bool(self.predictor_corrector)
        self.solver = FracStepStrategy( self.model_part, solver_configuration, self.ReformDofAtEachIteration, self.vel_toll, self.press_toll, self.max_vel_its, self.max_press_its, self.time_order, self.domain_size,self.predictor_corrector)


        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
        
   
    def Solve(self):
	

	self.timer.Start("Update Fixed Velocity Values")
        (self.solver).Solve()
	#self.timer.Stop("Update Fixed Velocity Values")
	print "####################"
	print "####################"
	
	print self.timer.Stop("Update Fixed Velocity Values")
	print self.timer
	print "####################"
	print "####################"

	

    def solve1(self):
	
 	self.timer.Start("primero")
        print self.model_part
        (self.solver).SolveStep3()
	self.timer.Stop("primero")

   
    def solve2(self):
	
 	self.timer.Start("segundo")
        print self.model_part
        (self.solver).SolveStep4()
	self.timer.Stop("segundo")

    def Reactions(self):
	
 	self.timer.Start("reactions")
	#sssssssssss
	print "ssssssssssssssssssssssssssss"
        print self.model_part
        (self.solver).Compute()
	self.timer.Stop("reactions")


    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True

       

    def WriteRestartFile(self,FileName):
        restart_file = open(FileName+".mdpa",'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintElements("Fluid2D",self.model_part.Elements,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,restart_file)
	new_restart_utilities.PrintRestart_ScalarVariable(ACCELERATION_X,"ACCELERATION_X",self.model_part.Nodes,restart_file)
	new_restart_utilities.PrintRestart_ScalarVariable(ACCELERATION_Y,"ACCELERATION_Y",self.model_part.Nodes,restart_file)
	new_restart_utilities.PrintRestart_ScalarVariable(ACCELERATION_Z,"ACCELERATION_Z",self.model_part.Nodes,restart_file)

	new_restart_utilities.PrintRestart_ScalarVariable(FORCE_X,"FORCE_X",self.model_part.Nodes,restart_file)
	new_restart_utilities.PrintRestart_ScalarVariable(FORCE_Y,"FORCE_Y",self.model_part.Nodes,restart_file)
	new_restart_utilities.PrintRestart_ScalarVariable(FORCE_Z,"FORCE_Z",self.model_part.Nodes,restart_file)
	#new_restart_utilities.PrintRestart_ScalarVariable(ANGULAR_ACCELERATION_X,"ANGULAR_ACCELERATION_X",self.model_part.Nodes,restart_file)
	#new_restart_utilities.PrintRestart_ScalarVariable(ANGULAR_ACCELERATION_Y,"ANGULAR_ACCELERATION_Y",self.model_part.Nodes,restart_file)
	#new_restart_utilities.PrintRestart_ScalarVariable(ANGULAR_ACCELERATION_Z,"ANGULAR_ACCELERATION_Z",self.model_part.Nodes,restart_file)
	#new_restart_utilities.PrintRestart_ScalarVariable(ANGULAR_VELOCITY_X,"ANGULAR_VELOCITY_X",self.model_part.Nodes,restart_file)
	#new_restart_utilities.PrintRestart_ScalarVariable(ANGULAR_VELOCITY_Y,"ANGULAR_VELOCITY_Y",self.model_part.Nodes,restart_file)
	#new_restart_utilities.PrintRestart_ScalarVariable(ANGULAR_VELOCITY_Z,"ANGULAR_VELOCITY_Z",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(PRESSURE,"PRESSURE",self.model_part.Nodes,restart_file)

        new_restart_utilities.PrintRestart_ScalarVariable(VISCOSITY,"VISCOSITY",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(DENSITY,"DENSITY",self.model_part.Nodes,restart_file)
        restart_file.close() 




        
        











