from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
CheckForPreviousImport()

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(NODAL_MASS) #MSI, i included the variable becouse i calculate energy
    # add displacements
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    # add dynamic variables
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    # add reactions for the displacements
    model_part.AddNodalSolutionStepVariable(REACTION)
    # add nodal force variables
    model_part.AddNodalSolutionStepVariable(INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(CONTACT_FORCE)
    # add specific variables for the problem conditions
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)
    model_part.AddNodalSolutionStepVariable(LINE_LOAD)
    model_part.AddNodalSolutionStepVariable(SURFACE_LOAD)
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION)

    print("*********************************************************************** ")
    print("Variables for the Static Structural Arc Length Solution added correctly")
    print("*********************************************************************** ")

def IncreasePointLoad(forcing_nodes_list, Load):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(POINT_LOAD_X, 0, Load[0])
        node.SetSolutionStepValue(POINT_LOAD_Y, 0, Load[1])
        node.SetSolutionStepValue(POINT_LOAD_Z, 0, Load[2])

def IncreaseDisplacement(forcing_nodes_list, disp):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(DISPLACEMENT_X, 0, disp[0])
        node.SetSolutionStepValue(DISPLACEMENT_Y, 0, disp[1])
        node.SetSolutionStepValue(DISPLACEMENT_Z, 0, disp[2])
        
def ChangeCondition(model_part):
    for node in model_part.Nodes:
        new_load = node.GetSolutionStepValue(POINT_LOAD) * model_part.ProcessInfo[LAMBDA];
        node.SetSolutionStepValue(POINT_LOAD, 0, new_load)

def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
    # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);

    print("*********************************************************************** ")
    print("Dofs for the Static Structural Arc Length Solution added correctly")
    print("*********************************************************************** ")

class StaticArcLengthStructuralSolver:
    #
    def __init__(self, model_part, domain_size):

        # Default settings
        self.echo_level = 0
        self.model_part = model_part
        self.domain_size = domain_size
        self.buffer_size = 3 #default buffer_size

        # Default varibles de Control de Arc Lenght Method
        self.Ide = 5
        self.factor_delta_lmax = 1.00
        self.max_iteration = 20
        self.toler = 1.0E-10
        self.norm = 1.0E-7
        self.MaxLineSearchIterations = 20
        self.tolls = 0.000001
        self.amp = 1.618
        self.etmxa = 5
        self.etmna = 0.1
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True

        # Default linear solver
        self.linear_solver = SkylineLUFactorizationSolver()

        # definition of the convergence criteria
        self.rel_disp_tol = 1e-4
        self.abs_disp_tol = 1e-9
        self.rel_res_tol = 1e-4
        self.abs_res_tol = 1e-9
        
        Displacement = DisplacementCriteria(self.rel_disp_tol, self.abs_disp_tol)
        Residual = ResidualCriteria(self.rel_res_tol, self.abs_res_tol)

        self.convergence_criterion_type = "And_criteria"
        self.mechanical_convergence_criterion = AndCriteria(Residual, Displacement)
    
        # Definition of the default builder_and_solver:
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.block_builder = False

        # Definition of the component wise calculation "computation is slower"
        #(it affects to strategy, builder_and_solver, scheme and convergence_criterion)
        self.component_wise = False
    
    #
    def Initialize(self):

        # Creating the convergence criterion:
        self.SetConvergenceCriterion()

        # Creating the solver
        self.solver = ResidualBasedArcLengthStrategy(self.model_part,  self.time_scheme, self.linear_solver, 
            self.mechanical_convergence_criterion,  self.Ide,  self.max_iteration,  self.factor_delta_lmax,  self.CalculateReactionFlag, self.ReformDofSetAtEachStep,  self.MoveMeshFlag )
    #
    def Solve(self):
        
        (self.solver).Solve()
        
        print("LAMBDA: ", self.model_part.ProcessInfo[LAMBDA])
        ChangeCondition(self.model_part)
    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
    #
    def SetConvergenceCriterion(self):

        # Mechanical convergence criteria
        D_RT = self.rel_disp_tol;
        D_AT = self.abs_disp_tol;
        R_RT = self.rel_res_tol;
        R_AT = self.abs_res_tol;

        if(self.echo_level > 1):
            print("::[Arc length Solver]:: CONVERGENCE CRITERION : ", self.convergence_criterion_type)

        if(self.convergence_criterion_type == "Displacement_criteria"):
            self.mechanical_convergence_criterion = DisplacementConvergenceCriterion(D_RT, D_AT)
        elif(self.convergence_criterion_type == "Residual_criteria"):
            if(self.component_wise):
                self.mechanical_convergence_criterion = ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
            else:
                self.mechanical_convergence_criterion = ResidualCriteria(R_RT, R_AT)
        elif(self.convergence_criterion_type == "And_criteria"):
            Displacement = DisplacementConvergenceCriterion(D_RT, D_AT)
            if(self.component_wise):
                Residual = ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
            else:
                Residual = ResidualCriteria(R_RT, R_AT)
            self.mechanical_convergence_criterion = AndCriteria(Residual, Displacement)
        elif(self.convergence_criterion_type == "Or_criteria"):
            Displacement = DisplacementConvergenceCriterion(D_RT, D_AT)
            if(self.component_wise):
                Residual = ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
            else:
                Residual = ResidualCriteria(R_RT, R_AT)
            self.mechanical_convergence_criterion = OrCriteria(Residual, Displacement)
        elif(self.convergence_criterion_type == "Mixed_criteria"):
            Displacement = MixedElementConvergeCriteria(D_RT, D_AT)
            Residual = ResidualCriteria(R_RT, R_AT)
            self.mechanical_convergence_criterion = AndCriteria(Residual, Displacement)
        
    #
    def ChangeCondition(self, model_part, lamda):
        for cond in model_part.Conditions:
            print(cond)
            
def CreateSolver(model_part, config):
    structural_solver = StaticArcLengthStructuralSolver(model_part, config.domain_size)
    model_part.ProcessInfo[LAMBDA] = 0.00;
    
    # Definition of the control variables of the Arc Lenght Method
    if(hasattr(config, "Ide")): 
        structural_solver.Ide = config.Ide
    if(hasattr(config, "factor_delta_lmax")):
        structural_solver.factor_delta_lmax = config.factor_delta_lmax
    if(hasattr(config, "max_iteration")):
        structural_solver.max_iteration = config.max_iteration
    if(hasattr(config, "toler")):
        structural_solver.toler = config.toler
    if(hasattr(config, "norm")):
        structural_solver.norm = config.norm
    if(hasattr(config, "MaxLineSearchIterations")):
        structural_solver.MaxLineSearchIterations = config.MaxLineSearchIterations
    if(hasattr(config, "tolls")):
        structural_solver.tolls = config.tolls
    if(hasattr(config, "amp")):
        structural_solver.amp = config.amp
    if(hasattr(config, "etmxa")):
        structural_solver.etmxa = config.etmxa
    if(hasattr(config, "etmna")): 
        structural_solver.etmna = config.etmna
    if(hasattr(config, "CalculateReactionFlag")): # bool
        structural_solver.CalculateReactionFlag = config.CalculateReactionFlag
    if(hasattr(config, "ReformDofSetAtEachStep")): # bool
        structural_solver.ReformDofSetAtEachStep = config.ReformDofSetAtEachStep
    if(hasattr(config, "MoveMeshFlag")): # bool
        structural_solver.MoveMeshFlag = config.MoveMeshFlag

    # Definition of the convergence criteria
    if(hasattr(config, "convergence_criterion")):
        structural_solver.convergence_criterion_type = config.convergence_criterion
    if(hasattr(config, "displacement_relative_tolerance")):
        structural_solver.rel_disp_tol = config.displacement_relative_tolerance
    if(hasattr(config, "displacement_absolute_tolerance")):
        structural_solver.abs_disp_tol = config.displacement_absolute_tolerance
    if(hasattr(config, "residual_relative_tolerance")):
        structural_solver.rel_res_tol = config.residual_relative_tolerance
    if(hasattr(config, "residual_absolute_tolerance")):
        structural_solver.abs_res_tol = config.residual_absolute_tolerance

    # definition of the echo level
    if(hasattr(config, "echo_level")):
        structural_solver.echo_level = config.echo_level

    # definition of the linear solver
    import linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        if(config.echo_level > 1):
            print("::[Mechanical Solver]:: LINEAR SOLVER : ", config.linear_solver_config.solver_type)
        structural_solver.linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)

        if(config.linear_solver_config.solver_type == "AMGCL"):
            structural_solver.block_builder = True
        else:
            structural_solver.block_builder = False
    
    return structural_solver
