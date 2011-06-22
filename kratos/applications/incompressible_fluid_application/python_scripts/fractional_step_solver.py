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

    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS);

    print "variables for the incompressible fluid solver added correctly"

def AddDofs(model_part):
  
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(PRESSURE);
        node.AddDof(FRACT_VEL_X);
        node.AddDof(FRACT_VEL_Y);
        node.AddDof(FRACT_VEL_Z);
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(VELOCITY_Z);

    print "dofs for the incompressible fluid solver added correctly"
    
##def ReadRestartFile(FileName,nodes):
##   aaa = __import__(FileName)
##   aaa.Restart(nodes)

def ReadRestartFile(FileName,nodes):
   NODES = nodes
   aaa = open(FileName)
   for line in aaa:
       exec(line)
       
##   import start.pyinc
   
##   aaa = __import__(FileName)
##   aaa.Restart(nodes)

   

class IncompressibleFluidSolver:
    
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
#        pILUPrecond = ILU0Preconditioner()
##        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-9, 5000,pILUPrecond)
        self.velocity_linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pILUPrecond)
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)

        self.dynamic_tau = 0.001


        ##handling slip condition
        self.slip_conditions_initialized = False
        self.create_slip_conditions = GenerateSlipConditionProcess(self.model_part,domain_size)

        self.compute_reactions=False
        


    def Initialize(self):
        (self.neighbour_search).Execute()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);

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
        self.solver = FractionalStepStrategy( self.model_part, solver_configuration, self.ReformDofAtEachIteration, self.vel_toll, self.press_toll, self.max_vel_its, self.max_press_its, self.time_order, self.domain_size,self.predictor_corrector)

        self.solver.ApplyFractionalVelocityFixity()

        ##generating the slip conditions
        self.create_slip_conditions.Execute()
        (self.solver).SetSlipProcess(self.create_slip_conditions);
        self.slip_conditions_initialized = True

        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
        
   
    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            self.solver.ApplyFractionalVelocityFixity()
            (self.neighbour_search).Execute()
            self.slip_conditions_initialized = False

        if(self.slip_conditions_initialized == False):
            self.create_slip_conditions.Execute()
            (self.solver).SetSlipProcess(self.create_slip_conditions);
            self.slip_conditions_initialized = True            

        print "just before solve"
        print self.model_part
        (self.solver).Solve()

        if(self.compute_reactions == True):
            self.solver.ComputeReactions(REACTION)

##        (self.create_slip_conditions).SetNormalVelocityToZero()
##        (self.create_slip_conditions).ApplyEdgeConstraints()

    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True

    def AdaptMesh(self):
        admissible_ratio = 0.05
        max_levels = 2
        refinement_utils = RefinementUtilities()
        if(self.domain_size == 2):
            raise "error refine in 2d not yet implemented"
        else:
            Refine = LocalRefineTetrahedraMesh(self.model_part)
        (self.model_part).ProcessInfo[FRACTIONAL_STEP]=10; ##just to be sure nothign is done
        refinement_utils.MarkForRefinement(ERROR_RATIO,self.model_part,admissible_ratio,max_levels)
        self.Clear()
        refine_on_reference = False;
        interpolate_internal_variables = False;
        Refine.LocalRefineMesh(refine_on_reference,interpolate_internal_variables)
        
        (self.neighbour_search).Execute()
        self.slip_conditions_initialized = False
        print "Refining finished"
        

    def WriteRestartFile(self,FileName):
        restart_file = open(FileName+".mdpa",'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintElements("Fluid3D",self.model_part.Elements,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(PRESSURE,"PRESSURE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VISCOSITY,"VISCOSITY",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(DENSITY,"DENSITY",self.model_part.Nodes,restart_file)
        restart_file.close() 

    def ActivateSmagorinsky(self,C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY,C)

    def ActivateSpalartAllmaras(self,wall_nodes,DES):
        from KratosFluidDynamicsApplication import *
        for node in wall_nodes:
            node.SetValue(IS_VISITED,1.0)

        distance_calculator = BodyDistanceCalculationUtils()
        distance_calculator.CalculateDistances2D(self.model_part.Elements,DISTANCE,100.0)

        non_linear_tol = 0.001
        max_it = 5
        reform_dofset = self.ReformDofAtEachIteration
        time_order = self.time_order
        pPrecond = DiagonalPreconditioner()
        turbulence_linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)
        turbulence_model = SpalartAllmarasTurbulenceModel(self.model_part,turbulence_linear_solver,self.domain_size,non_linear_tol,max_it,reform_dofset,time_order);
        turbulence_model.AdaptForFractionalStep()
        if(DES==True):
            turbulence_model.ActivateDES(1.0);

        self.solver.AddInitializeIterationProcess(turbulence_model);



        
        

