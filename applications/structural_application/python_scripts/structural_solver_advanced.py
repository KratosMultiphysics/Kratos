#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import sys

import structural_solver_static

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(ELASTIC_LEFT_CAUCHY_GREEN_OLD);
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS);
    model_part.AddNodalSolutionStepVariable(PRESTRESS);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL_DT);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS_DT);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_NULL_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_EINS_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION_AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_DT);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_DT);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_DT);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(EXCESS_PORE_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    #auxiliary variables misused for mesh rezoning ;-)
    model_part.AddNodalSolutionStepVariable(IS_VISITED);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(LAGRANGE_AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(LAGRANGE_WATER_PRESSURE);
    #model_part.AddNodalSolutionStepVariable(INTERNAL_VARIABLES);
    model_part.AddNodalSolutionStepVariable(MOMENTUM);
    model_part.AddNodalSolutionStepVariable(ROTATION);
    model_part.AddNodalSolutionStepVariable(PRESSURE);        
    model_part.AddNodalSolutionStepVariable(ERROR_RATIO);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    print "variables for the dynamic structural solution added correctly"
    
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs 
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(WATER_PRESSURE);
        node.AddDof(AIR_PRESSURE);
        node.AddDof(LAGRANGE_DISPLACEMENT_X);
        node.AddDof(LAGRANGE_DISPLACEMENT_Y);
        node.AddDof(LAGRANGE_DISPLACEMENT_Z);
        #node.AddDof(ROTATION_X);
        #node.AddDof(ROTATION_Y);
        #node.AddDof(ROTATION_Z);
        #node.AddDof(LAGRANGE_AIR_PRESSURE);
        #node.AddDof(LAGRANGE_WATER_PRESSURE);
    print "dofs for the dynamic structural solution added correctly"
        
#######################################################################
class SolverAdvanced(structural_solver_static.StaticStructuralSolver):
    def __init__( self, model_part, domain_size, time_steps, analysis_parameters, abs_tol, rel_tol ):
        structural_solver_static.StaticStructuralSolver.__init__( self, model_part, domain_size )
        self.time_steps = time_steps
        self.analysis_parameters = analysis_parameters
        self.echo_level = 0
        self.damp_factor = 1.0
        self.toll = rel_tol
        self.absolute_tol = abs_tol
        #definition of the solvers
        self.structure_linear_solver =  SkylineLUFactorizationSolver()
        #pDiagPrecond = ParallelDiagonalPreconditioner()
        #self.structure_linear_solver =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(0.000001,1e-9)
        #self.conv_criteria = ParallelDisplacementCriteria(0.000001,1e-9)
        self.CalculateReactionFlag = False
        #######################################################################
    def Initialize(self):
        #definition of time integration scheme
        if( self.time_steps == 1 ):
            self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
            #self.time_scheme = ParallelResidualBasedIncrementalUpdateStaticScheme()
            self.MoveMeshFlag = True
        else:
            print "using newmark scheme"
            self.time_scheme = ResidualBasedNewmarkScheme(self.damp_factor)
            #self.time_scheme = ParallelResidualBasedNewmarkScheme(self.damp_factor)
            self.MoveMeshFlag = True
        #definition of the convergence criteria
        self.conv_criteria = MultiPhaseFlowCriteria(self.toll,self.absolute_tol)
        #self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
        builder_and_solver = ResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)
        #builder_and_solver = MultiPhaseBuilderAndSolver(self.structure_linear_solver)
        #builder_and_solver = ParallelResidualBasedEliminationBuilderAndSolverDeactivation(self.structure_linear_solver)
        #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
        #creating the solution strategy
        self.ReformDofSetAtEachStep = True
        #KLUDGE: this has to be True!
        self.MoveMeshFlag = True
        self.space_utils = UblasSparseSpace()
        #self.space_utils = ParallelUblasSparseSpace()
        #importing strategy
        #import ekate_strategy
        import uzawa_contact_strategy
        self.solver = uzawa_contact_strategy.SolvingStrategyPython( self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, builder_and_solver )

