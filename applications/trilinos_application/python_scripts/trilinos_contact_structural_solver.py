#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
from KratosTrilinosApplication import *
import mpi
import sys

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
    #auxiliary variables misused for mesh rezoning ;-)
    model_part.AddNodalSolutionStepVariable(IS_VISITED);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(LAGRANGE_DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(LAGRANGE_AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(LAGRANGE_WATER_PRESSURE);
    #model_part.AddNodalSolutionStepVariable(INTERNAL_VARIABLES);
    model_part.AddNodalSolutionStepVariable(MOMENTUM);
    model_part.AddNodalSolutionStepVariable(PRESSURE);        
    model_part.AddNodalSolutionStepVariable(ERROR_RATIO);
    model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);
    print "variables for the dynamic structural solution added correctly"
    
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs 
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(LAGRANGE_DISPLACEMENT_X);
        node.AddDof(LAGRANGE_DISPLACEMENT_Y);
        node.AddDof(LAGRANGE_DISPLACEMENT_Z);
        #node.AddDof(LAGRANGE_AIR_PRESSURE);
        #node.AddDof(LAGRANGE_WATER_PRESSURE);
    print "dofs for the dynamic structural solution added correctly"
        
#######################################################################
class SolverAdvanced:
    def __init__( self, model_part, domain_size, time_steps, analysis_parameters, abs_tol, rel_tol, application_path ):
        self.model_part = model_part
        self.time_steps = time_steps
        self.analysis_parameters = analysis_parameters

        self.echo_level = 0
        self.damp_factor = 1.0
        self.toll = rel_tol
        self.absolute_tol = abs_tol


        self.time_scheme = TrilinosResidualBasedIncrementalUpdateStaticScheme()

        self.Comm = CreateCommunicator()

        #self.buildertype="ML3D"
        #self.buildertype="superludist"
        #self.buildertype="MLdeactivation"
        self.buildertype="superludist_deactivation"

        #definition of the solvers
        #self.structure_linear_solver =  TrilinosLinearSolver()
        self.solver_parameters = ParameterList()
        self.preconditioner_parameters = ParameterList()
        #self.structure_linear_solver =  AmesosSolver("Superludist",self.solver_parameters);
        self.solver_parameters.set("AZ_solver","AZ_gmres")
        self.solver_parameters.set("AZ_kspace",100)
        self.solver_parameters.set("AZ_output",32)
        self.solver_parameters.set("AZ_precond","AZ_none")
        self.structure_linear_solver = AztecSolver(self.solver_parameters,"Amesos",self.preconditioner_parameters,1.0e-9,300,1);


        
        #definition of the convergence criteria
        self.conv_criteria = TrilinosDisplacementCriteria(1e-6,1e-9,self.Comm)

        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True
        self.calculate_norm_dx_flag = False
        self.max_iterations = 10

        if(domain_size == 2):
            self.guess_row_size = 15
        else:
            self.guess_row_size = 45
            

        self.guess_row_size = 18


        #######################################################################
    def Initialize(self):
        #definition of time integration scheme
        if( self.time_steps == 1 ):
            self.time_scheme = TrilinosResidualBasedIncrementalUpdateStaticScheme()
            self.MoveMeshFlag = False
        else:
            print "using newmark scheme"
            self.time_scheme = TrilinosResidualBasedNewmarkScheme(self.damp_factor)
            self.MoveMeshFlag = True
        #definition of the convergence criteria
        self.conv_criteria = TrilinosDisplacementCriteria(1e-6,1e-9,self.Comm)
        #self.conv_criteria = MultiPhaseFlowCriteria(self.toll,self.absolute_tol)
        #self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
        #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
        #creating the solution strategy
        self.ReformDofSetAtEachStep = True
        #KLUDGE: this has to be True!
        self.MoveMeshFlag = True
        self.space_utils = TrilinosSparseSpace()
        
        import trilinos_contact_strategy
        self.solver = trilinos_contact_strategy.SolvingStrategyPython( self.buildertype, self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag, self.analysis_parameters, self.space_utils, self.Comm, self.guess_row_size )

        
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()
        
    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
