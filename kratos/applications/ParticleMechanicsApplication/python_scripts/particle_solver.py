from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
#from KratosMultiphysics.PfemSolidMechanicsApplication import *
from KratosMultiphysics.ParticleMechanicsApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, config=None):  
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
    #model_part.AddNodalSolutionStepVariable(CONTACT_FORCE)  
     
    # add specific variables for the problem conditions
    #model_part.AddNodalSolutionStepVariable(IMPOSED_DISPLACEMENT)
    #model_part.AddNodalSolutionStepVariable(IMPOSED_ROTATION)
    #model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    #model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)
    model_part.AddNodalSolutionStepVariable(LINE_LOAD)
    model_part.AddNodalSolutionStepVariable(SURFACE_LOAD)
    #model_part.AddNodalSolutionStepVariable(POINT_TORQUE)
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(NODAL_MOMENTUM)
    model_part.AddNodalSolutionStepVariable(NODAL_INERTIA)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_AUX)
    model_part.AddNodalSolutionStepVariable(AUX_VELOCITY)
    model_part.AddNodalSolutionStepVariable(AUX_ACCELERATION)
    if config is not None:
        if hasattr(config, "RotationDofs"):
            if config.RotationDofs:
                # add specific variables for the problem (rotation dofs)
                model_part.AddNodalSolutionStepVariable(ROTATION);
                model_part.AddNodalSolutionStepVariable(TORQUE);
                model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY)
                model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION)
        if hasattr(config, "PressureDofs"):
            if config.PressureDofs:
                # add specific variables for the problem (pressure dofs)
                model_part.AddNodalSolutionStepVariable(PRESSURE)
                model_part.AddNodalSolutionStepVariable(PRESSURE_REACTION);
        if hasattr(config, "time_integration_method"):
            if config.time_integration_method == "Explicit" :
                model_part.AddNodalSolutionStepVariable(NODAL_MASS)
                model_part.AddNodalSolutionStepVariable(FORCE_RESIDUAL)
                model_part.AddNodalSolutionStepVariable(MIDDLE_VELOCITY)
                        
    print("::[Particle Solver]:: Variables ADDED")


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);

    if config is not None:
        if hasattr(config, "RotationDofs"):
            if config.RotationDofs:
                for node in model_part.Nodes:
                    node.AddDof(ROTATION_X, TORQUE_X);
                    node.AddDof(ROTATION_Y, TORQUE_Y);
                    node.AddDof(ROTATION_Z, TORQUE_Z);
        if hasattr(config, "PressureDofs"):
            if config.PressureDofs:
                for node in model_part.Nodes:
                    node.AddDof(PRESSURE, PRESSURE_REACTION);

    print("::[Particle Solver]:: DOF's ADDED")


class ParticleSolver:
    #

    def __init__(self, model_part1, model_part2, new_element, domain_size, geometry_element):

        # default settings
        
        self.echo_level = 2
        self.model_part1 = model_part1  #grid_model_part
        self.model_part2 = model_part2  #mpm_model_part
        
        self.new_element = new_element
        
        self.domain_size = domain_size
        self.geometry_element = geometry_element
        self.buffer_size = 3 #default buffer_size

        self.pressure_dofs = False
        self.rotation_dofs = False
        self.stabilization_factor = 1.0

        # definition of the solvers
        self.scheme_type = "DynamicSolver"
        self.time_integration_method = "Implicit" 
        self.explicit_integration_scheme = "CentralDifferences" #for explicit strategies
        
        self.arlequin = 0
        
        self.linear_solver = SkylineLUFactorizationSolver()

        # definition of the convergence criteria
        self.rel_disp_tol = 1e-4
        self.abs_disp_tol = 1e-9
        self.rel_res_tol = 1e-4
        self.abs_res_tol = 1e-9
        self.max_iters = 30

        self.convergence_criterion_type = "Residual_criteria"
        self.particle_convergence_criterion = ResidualCriteria(self.rel_res_tol, self.abs_res_tol)
        # self.mechanical_convergence_criterion = DisplacementCriteria(self.rel_disp_tol,self.abs_disp_tol)
        # self.mechanical_convergence_criterion = ComponentWiseResidualConvergenceCriterion(self.rel_res_tol,self.abs_res_tol)
        # self.mechanical_convergence_criterion = DisplacementConvergenceCriterion(self.rel_disp_tol,self.abs_disp_tol)

        # definition of the default builder_and_solver:
        #self.block_builder = False
        #self.builder_and_solver = ResidualBasedBuilderAndSolver(self.linear_solver)

        # definition of the component wise calculation "computation is slower"
        #(it affects to strategy, builder_and_solver, scheme and convergence_criterion)
        self.component_wise = False

        # definition of computing flags
        self.compute_reactions = True
        self.compute_contact_forces = False
        self.line_search = False
        self.implex = False
        self.reform_step_dofs = True

        # definition of the (x_n+1 = x_n+dx) in each iteration -> strategy base class includes de MoveMesh method
        self.move_mesh_flag = False
        
        self.max_delta_time = 1e-5
        self.fraction_delta_time = 0.9
        self.time_step_prediction_level = 0
        self.rayleigh_damping = False

        print("::[Particle Solver]:: -START-")

    #
    
    def DefineBufferSize(self):
      
      if( self.time_integration_method == "Explicit"):
        if(self.explicit_integration_scheme == "CentralDifferences"):
            self.buffer_size = 3 
            #could be 1 becouse it uses a extra variable - MIDDLE_VELOCITY for previous time step velocity and avoids extra buffer_size. However
            #the prediction/application of B.C. need last step values 
        if(self.arlequin == 1):
            self.buffer_size = 3 
                  
    def Initialize(self):

       
        if(self.time_integration_method == "Implicit"):

          # creating the convergence criterion:
          self.SetConvergenceCriterion()

          # creating the solution strategy (application strategy or python strategy)

          # self.reform_step_dofs = False;
          print('set solving strategy')
          #if(self.component_wise):
              #self.particle_solver = ComponentWiseNewtonRaphsonStrategy(self.model_part1, self.particle_scheme, self.linear_solver, self.particle_convergence_criterion, self.builder_and_solver, self.max_iters, self.compute_reactions, self.reform_step_dofs, self.move_mesh_flag)
          
         # else:
              #if(self.line_search):
                  #self.particle_solver = ResidualBasedNewtonRaphsonLineSearchStrategy(self.model_part1, self.particle_scheme, self.linear_solver, self.particle_convergence_criterion, self.builder_and_solver, self.max_iters, self.compute_reactions, self.reform_step_dofs, self.move_mesh_flag)
              #else:
        if(self.domain_size==2):
            self.particle_solver = MPM2D(self.model_part1, self.model_part2, self.linear_solver, self.new_element, self.move_mesh_flag, self.scheme_type, self.geometry_element)
        else:
            self.particle_solver = MPM3D(self.model_part1, self.model_part2, self.linear_solver, self.new_element, self.move_mesh_flag, self.scheme_type, self.geometry_element)
        #if(self.time_integration_method == "Explicit"):
            #self.particle_solver = ExplicitStrategy(self.model_part1, self.particle_scheme, self.linear_solver, self.compute_reactions, self.reform_step_dofs, self.move_mesh_flag)
          
            #self.particle_solver.SetRebuildLevel(0) # 1 to recompute the mass matrix in each explicit step 
        
        
        # creating the builder and solver:
        print('set builder and solver')
        # creating the solution scheme:
        #self.SetBuilderAndSolver()
        print('set solution scheme')
        #self.SetSolutionScheme()
       
        self.SetStabilizationFactor(self.stabilization_factor)
        
        (self.particle_solver).SetEchoLevel(self.echo_level)
        
        # check if everything is assigned correctly
        print('check that everything is correctly assigned')
        self.Check()
        
        print("::[Particle Solver INITIALIZATION]:: -END- ")

    #
    def Solve(self):
        print("In the solver ")
        (self.particle_solver).Solve()

    #
    def SetEchoLevel(self, level):
        self.echo_level = level
        (self.particle_solver).SetEchoLevel(level)

    #
    def SetStabilizationFactor(self, factor):
        self.model_part1.ProcessInfo[STABILIZATION_FACTOR] = factor

    #
    def SetRestart(self, load_restart):
        # check if is a restart file is loaded
        if(load_restart):
            # set solver as initialized if is a run which is restarted
            self.particle_solver.SetInitializePerformedFlag(True)
        else:
            # initialize strategy solver
            print('initializing the solver, load_restart = false')
            self.particle_solver.Initialize()
        
            
    #
    def Clear(self):
        (self.solver).Clear()

    #
    def Check(self):
        self.particle_solver.Check();

    #
    def SetBuilderAndSolver(self):

        # type of builder and solver
        if(self.component_wise):
            self.builder_and_solver = ComponentWiseBuilderAndSolver(self.linear_solver)
        else:
            if(self.block_builder):
                # to keep matrix blocks in builder
                self.builder_and_solver = BlockResidualBasedBuilderAndSolver(self.linear_solver)
            else:
                self.builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)

    #
    def SetSolutionScheme(self):

        if(self.time_integration_method == "Explicit"):
            
            if(self.explicit_integration_scheme == "CentralDifferences"):

                  self.particle_scheme = ExplicitCentralDifferencesScheme(self.max_delta_time, self.fraction_delta_time, self.time_step_prediction_level, self.rayleigh_damping)
          
        else:
          
          # type of solver for implicit solution method (static,dynamic,quasi-static,pseudo-dynamic)
               
          if(self.scheme_type == "StaticSolver"):
              self.particle_scheme = ResidualBasedIncrementalUpdateStaticScheme()   
          elif(self.scheme_type == "DynamicSolver"):
              # definition of time scheme
              self.damp_factor_f = 0.00;
              self.damp_factor_m = -0.01;
              self.dynamic_factor = 1;
              
              self.model_part1.ProcessInfo[RAYLEIGH_ALPHA] = 0.0
              self.model_part1.ProcessInfo[RAYLEIGH_BETA ] = 0.0

              if(self.component_wise):
                  self.particle_scheme = ComponentWiseBossakScheme(self.damp_factor_m, self.dynamic_factor)
              else:
                  if(self.compute_contact_forces):
                      self.particle_scheme = ResidualBasedContactBossakScheme(self.damp_factor_m, self.dynamic_factor)
                  else:
                      self.particle_scheme = ResidualBasedBossakScheme(self.damp_factor_m, self.dynamic_factor)
          elif(self.scheme_type == "QuasiStaticSolver"):
              # definition of time scheme
              self.damp_factor_f = 0.00;
              self.damp_factor_m = 0.00;
              self.dynamic_factor = 0;  # quasi-static
              if(self.component_wise):
                  self.particle_scheme = ComponentWiseBossakScheme(self.damp_factor_m, self.dynamic_factor)
              else:
                  if(self.compute_contact_forces):
                      self.particle_scheme = ResidualBasedContactBossakScheme(self.damp_factor_m, self.dynamic_factor)
                  else:
                      self.particle_scheme = ResidualBasedBossakScheme(self.damp_factor_m, self.dynamic_factor)
                  # self.particle_scheme = ResidualBasedNewmarkScheme(self.dynamic_factor)
          elif(self.scheme_type == "PseudoDynamicSolver"):
              # definition of time scheme
              self.damp_factor_f = -0.3;
              self.damp_factor_m = 10.0;
              self.particle_scheme = ResidualBasedRelaxationScheme(self.damp_factor_f, self.damp_factor_m)

    #
    def SetConvergenceCriterion(self):

        # particle convergence criteria
        D_RT = self.rel_disp_tol;
        D_AT = self.abs_disp_tol;
        R_RT = self.rel_res_tol;
        R_AT = self.abs_res_tol;

        if(self.rotation_dofs):

            if(self.convergence_criterion_type == "Displacement_criteria"):
                self.particle_convergence_criterion = DisplacementCriteria(D_RT, D_AT)
            elif(self.convergence_criterion_type == "Residual_criteria"):
                self.particle_convergence_criterion = ResidualCriteria(R_RT, R_AT)
            elif(self.convergence_criterion_type == "And_criteria"):
                Displacement = DisplacementCriteria(D_RT, D_AT)
                Residual = ResidualCriteria(R_RT, R_AT)
                self.particle_convergence_criterion = AndCriteria(Residual, Displacement)
            elif(self.convergence_criterion_type == "Or_criteria"):
                Displacement = DisplacementCriteria(D_RT, D_AT)
                Residual = ResidualCriteria(R_RT, R_AT)
                self.particle_convergence_criterion = OrCriteria(Residual, Displacement)
            elif(self.convergence_criterion_type == "Mixed_criteria"):
                Displacement = MixedElementCriteria(D_RT, D_AT)
                Residual = ResidualCriteria(R_RT, R_AT)
                self.particle_convergence_criterion = AndCriteria(Residual, Displacement)

        else:

            if(self.echo_level > 1):
                print("::[Particle Solver]:: CONVERGENCE CRITERION : ", self.convergence_criterion_type)

            if(self.convergence_criterion_type == "Displacement_criteria"):
                self.particle_convergence_criterion = DisplacementConvergenceCriterion(D_RT, D_AT)
            elif(self.convergence_criterion_type == "Residual_criteria"):
                if(self.component_wise):
                    self.particle_convergence_criterion = ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    self.particle_convergence_criterion = ResidualCriteria(R_RT, R_AT)
            elif(self.convergence_criterion_type == "And_criteria"):
                Displacement = DisplacementConvergenceCriterion(D_RT, D_AT)
                if(self.component_wise):
                    Residual = ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = ResidualCriteria(R_RT, R_AT)
                self.particle_convergence_criterion = AndCriteria(Residual, Displacement)
            elif(self.convergence_criterion_type == "Or_criteria"):
                Displacement = DisplacementConvergenceCriterion(D_RT, D_AT)
                if(self.component_wise):
                    Residual = ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = ResidualCriteria(R_RT, R_AT)
                self.particle_convergence_criterion = OrCriteria(Residual, Displacement)
            elif(self.convergence_criterion_type == "Mixed_criteria"):
                Displacement = MixedElementConvergeCriteria(D_RT, D_AT)
                Residual = ResidualCriteria(R_RT, R_AT)
                self.particle_convergence_criterion = AndCriteria(Residual, Displacement)


#
#
def CreateSolver(model_part1, model_part2, new_element, config, geometry_element):

    structural_solver = ParticleSolver(model_part1, model_part2,new_element, config.domain_size,geometry_element)

    #Explicit scheme parameters
    if(hasattr(config, "max_delta_time")):
        structural_solver.max_delta_time = config.max_delta_time
    if(hasattr(config, "time_step_prediction_level")):
        value = 0
        if(str(config.time_step_prediction_level) == "Automatic"):
          value = 1
        elif(str(config.time_step_prediction_level) == "RefreshEveryTimeStep"):
          value = 2
        structural_solver.time_step_prediction_level = value

    # definition of the convergence criteria
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
    if(hasattr(config, "max_iteration")):
        structural_solver.max_iters = config.max_iteration

    # definition of the global solver type
    if(hasattr(config, "scheme_type")):
        structural_solver.scheme_type = config.scheme_type
    if(hasattr(config, "time_integration_method")):
        structural_solver.time_integration_method = config.time_integration_method
    if(hasattr(config, "explicit_integration_scheme")):
        structural_solver.explicit_integration_scheme = config.explicit_integration_scheme
        if(hasattr(config, "rayleigh_damping")):
            structural_solver.rayleigh_damping = config.rayleigh_damping
    if(hasattr(config, "arlequin")):
        structural_solver.arlequin = config.arlequin

    structural_solver.DefineBufferSize()
    
    # definition of the solver parameters
    if(hasattr(config, "ComputeReactions")):
        structural_solver.compute_reactions = config.ComputeReactions  # bool
    if(hasattr(config, "ComputeContactForces")):
        structural_solver.compute_contact_forces = config.ComputeContactForces  # bool
    if(hasattr(config, "ReformDofSetAtEachStep")):
        structural_solver.reform_step_dofs = config.ReformDofSetAtEachStep  # bool
    if(hasattr(config, "RotationDofs")):
        structural_solver.rotation_dofs = config.RotationDofs  # bool
    if(hasattr(config, "PressureDofs")):
        structural_solver.pressure_dofs = config.PressureDofs  # bool
    if(hasattr(config, "LineSearch")):
        structural_solver.line_search = config.LineSearch  # bool
    if(hasattr(config, "Implex")):
        structural_solver.implex = config.Implex  # bool
    if(hasattr(config, "ComponentWise")):
        structural_solver.component_wise = config.ComponentWise  # bool
    if(hasattr(config, "StabilizationFactor")):
        structural_solver.stabilization_factor = config.StabilizationFactor #double

    # definition of the echo level
    if(hasattr(config, "echo_level")):
        structural_solver.echo_level = config.echo_level

    # definition of the linear solver
    import linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        if(config.echo_level > 1):
            print("::[Particle Solver]:: LINEAR SOLVER : ", config.linear_solver_config.solver_type)
        structural_solver.linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)

        if(config.linear_solver_config.solver_type == "AMGCL"):
            structural_solver.block_builder = True
        else:
            structural_solver.block_builder = False

    return structural_solver
