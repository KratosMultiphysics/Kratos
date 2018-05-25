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
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    
    # add reactions for the displacements
    model_part.AddNodalSolutionStepVariable(REACTION)
    
    # add nodal force variables
    model_part.AddNodalSolutionStepVariable(INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_FORCE)
    #model_part.AddNodalSolutionStepVariable(CONTACT_FORCE)  
     
    # add specific variables for the problem conditions
    #model_part.AddNodalSolutionStepVariable(IMPOSED_DISPLACEMENT)
    #model_part.AddNodalSolutionStepVariable(IMPOSED_ROTATION)
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    #model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)
    #model_part.AddNodalSolutionStepVariable(LINE_LOAD)
    #model_part.AddNodalSolutionStepVariable(SURFACE_LOAD)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    #model_part.AddNodalSolutionStepVariable(POINT_TORQUE)
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(NODAL_MOMENTUM)
    model_part.AddNodalSolutionStepVariable(NODAL_INERTIA)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_AUX)
    model_part.AddNodalSolutionStepVariable(AUX_VELOCITY)
    model_part.AddNodalSolutionStepVariable(AUX_ACCELERATION)
    model_part.AddNodalSolutionStepVariable(AUX_R)
    model_part.AddNodalSolutionStepVariable(AUX_T)
    model_part.AddNodalSolutionStepVariable(AUX_R_VEL)
    model_part.AddNodalSolutionStepVariable(AUX_T_VEL)
    model_part.AddNodalSolutionStepVariable(AUX_R_ACC)
    model_part.AddNodalSolutionStepVariable(AUX_T_ACC)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(NODAL_LUMPED_MASS)
    
    # add for slope with slips
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    

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
                # model_part.AddNodalSolutionStepVariable(PRESSURE)
                model_part.AddNodalSolutionStepVariable(PRESSURE_REACTION);
                model_part.AddNodalSolutionStepVariable(NODAL_MPRESSURE)
                model_part.AddNodalSolutionStepVariable(AUX_PRESSURE)
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

    def __init__(self, model_part1, model_part2, model_part3, new_element, domain_size, geometry_element, number_particle):
        
        Logger.PrintInfo("ParticleSolver", "This Solver is deprecated and could be removed sometime soon.")

        # default settings
        
        self.echo_level = 2
        self.model_part1 = model_part1  #grid_model_part
        self.model_part2 = model_part2  #initial_model_part
        self.model_part3 = model_part3  #mpm_model_part
        
        self.new_element = new_element
        
        self.domain_size = domain_size
        self.geometry_element = geometry_element
        self.number_particle = number_particle

        self.buffer_size = 3 #default buffer_size
        #self.linear_solver = SkylineLUFactorizationSolver()
      
        self.move_mesh_flag = False
        
     

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

       
       
        if(self.domain_size==2):
            self.particle_solver = MPM2D(self.model_part1, self.model_part2, self.model_part3, self.linear_solver, self.new_element, self.move_mesh_flag, self.scheme_type, self.geometry_element, self.number_particle, self.block_builder)
        else:
            self.particle_solver = MPM3D(self.model_part1, self.model_part2, self.model_part3, self.linear_solver, self.new_element, self.move_mesh_flag, self.scheme_type, self.geometry_element,  self.number_particle, self.block_builder)
      
        
        (self.particle_solver).SetEchoLevel(self.echo_level)
        
        # check if everything is assigned correctly
        print('check that everything is correctly assigned')
        self.Check()
        
        print("::[Particle Solver INITIALIZATION]:: -END- ")

    #
    def Solve(self):
        #print("In the solver ")
        (self.particle_solver).Solve()

    #
    def SetEchoLevel(self, level):
        self.echo_level = level
        (self.particle_solver).SetEchoLevel(level)

    #
    
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
  
    

#
#
def CreateSolver(model_part1, model_part2, model_part3, new_element, config, geometry_element, number_particle):

    structural_solver = ParticleSolver(model_part1, model_part2, model_part3, new_element, config.domain_size,geometry_element, number_particle)

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
        
        if(config.linear_solver_config.solver_type == "AMGCL"):
            params = Parameters( """ {
                                    "solver_type" : "AMGCL",
                                    "smoother_type":"damped_jacobi",
                                    "krylov_type": "cg",
                                    "coarsening_type": "aggregation",
                                    "max_iteration": 200,
                                    "provide_coordinates": false,
                                    "gmres_krylov_space_dimension": 100,
                                    "verbosity" : 2,
                                    "tolerance": 1e-7,
                                    "scaling": false,
                                    "block_size": 3,
                                    "use_block_matrices_if_possible" : true,
                                    "coarse_enough" : 50
                                } """)
            structural_solver.linear_solver = linear_solver_factory.ConstructSolver(params)      
            structural_solver.block_builder = True
            
        else:
            structural_solver.linear_solver = linear_solver_factory.ConstructSolver(config.linear_solver_config)
            structural_solver.block_builder = False

    return structural_solver
