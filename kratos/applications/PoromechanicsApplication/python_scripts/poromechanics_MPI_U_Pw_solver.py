from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import poromechanics_U_Pw_solver


def CreateSolver(main_model_part, custom_settings):
    
    return MPIUPwSolver(main_model_part, custom_settings)


class MPIUPwSolver(poromechanics_U_Pw_solver.UPwSolver):

    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "poromechanics_MPI_U_Pw_solver",
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "buffer_size": 2,
            "echo_level": 0,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": true,
            "solution_type": "Quasi-Static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "rayleigh_m": 0.0,
            "rayleigh_k": 0.0,
            "strategy_type": "Newton-Raphson",
            "convergence_criterion": "Displacement_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 15,
            "desired_iterations": 4,
            "max_radius_factor": 20.0,
            "min_radius_factor": 0.5,
            "block_builder": true,
            "nonlocal_damage": false,
            "characteristic_length": 0.05,
            "search_neighbours_step": false,
            "linear_solver_settings":{
                "solver_type": "AmgclMPISolver",
                "tolerance": 1.0e-6,
                "max_iteration": 200,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "None",
                "krylov_type": "fgmres"
            },
            "problem_domain_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "loads_sub_model_part_list": [""],
            "loads_variable_list": [""]
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        # Construct the linear solver
        import trilinos_linear_solver_factory
        self.linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of MPI UPwSolver finished")

    def AddVariables(self):
        
        super(MPIUPwSolver, self).AddVariables()
        
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    def ImportModelPart(self):
        
        # Construct the Trilinos import model part utility
        import trilinos_import_model_part_utility
        TrilinosModelPartImporter = trilinos_import_model_part_utility.TrilinosImportModelPartUtility(self.main_model_part, self.settings)

        # Execute the Metis partitioning and reading
        TrilinosModelPartImporter.ExecutePartitioningAndReading()

        # Create computing_model_part, set constitutive law and buffer size
        self._ExecuteAfterReading()

        # Construct the communicators
        TrilinosModelPartImporter.CreateCommunicators()

    def Initialize(self):
        
        # Construct the communicator
        self.EpetraCommunicator = TrilinosApplication.CreateCommunicator()
        
        # Builder and solver creation
        builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        scheme = self._ConstructScheme(self.settings["scheme_type"].GetString(),
                                         self.settings["solution_type"].GetString())
                
        # Get the convergence criterion
        convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())
                
        # Solver creation
        self.Solver = self._ConstructSolver(builder_and_solver,
                                            scheme,
                                            convergence_criterion,
                                            self.settings["strategy_type"].GetString())

        # Set echo_level
        self.Solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Solver.Check();

        print ("Initialization MPI UPwSolver finished")

    #### Specific internal functions ####

    def _ConstructBuilderAndSolver(self, block_builder):

        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 15
        else:
            guess_row_size = 45

        if(block_builder):
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                   guess_row_size,
                                                                                   self.linear_solver
                                                                                   )
        else:
            builder_and_solver = TrilinosApplication.TrilinosEliminationBuilderAndSolver(self.EpetraCommunicator,
                                                                                         guess_row_size,
                                                                                         self.linear_solver
                                                                                         )

        return builder_and_solver

    def _ConstructScheme(self, scheme_type, solution_type):
        
        if(scheme_type == "Newmark"):
            beta = self.settings["newmark_beta"].GetDouble()
            gamma = self.settings["newmark_gamma"].GetDouble()
            theta = self.settings["newmark_theta"].GetDouble()
            rayleigh_m = self.settings["rayleigh_m"].GetDouble()
            rayleigh_k = self.settings["rayleigh_k"].GetDouble()
            if(solution_type == "Quasi-Static"):
                if(rayleigh_m<1.0e-20 and rayleigh_k<1.0e-20):
                    scheme = KratosPoro.TrilinosNewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    scheme = KratosPoro.TrilinosNewmarkQuasistaticDampedUPwScheme(beta,gamma,theta,rayleigh_m,rayleigh_k)
            else:
                scheme = KratosPoro.TrilinosNewmarkDynamicUPwScheme(beta,gamma,theta,rayleigh_m,rayleigh_k)
        else:
            raise Exception("Apart from Newmark, other scheme_type are not available.")
        
        return scheme

    def _ConstructConvergenceCriterion(self, convergence_criterion):
        
        D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["echo_level"].GetInt()
        
        if(convergence_criterion == "Displacement_criterion"):
            convergence_criterion = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT, self.EpetraCommunicator)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "Residual_criterion"):
            convergence_criterion = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "And_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT, self.EpetraCommunicator)
            Displacement.SetEchoLevel(echo_level)
            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = TrilinosApplication.TrilinosAndCriteria(Residual, Displacement)
        elif(convergence_criterion == "Or_criterion"):
            Displacement = TrilinosApplication.TrilinosDisplacementCriteria(D_RT, D_AT, self.EpetraCommunicator)
            Displacement.SetEchoLevel(echo_level)
            Residual = TrilinosApplication.TrilinosResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = TrilinosApplication.TrilinosOrCriteria(Residual, Displacement)
        
        return convergence_criterion
    
    def _ConstructSolver(self, builder_and_solver, scheme, convergence_criterion, strategy_type):
        
        max_iters = self.settings["max_iteration"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()
                
        if strategy_type == "Newton-Raphson":
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
            solver = TrilinosApplication.TrilinosNewtonRaphsonStrategy(self.main_model_part,
                                                                       scheme,
                                                                       self.linear_solver,
                                                                       convergence_criterion,
                                                                       builder_and_solver,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag
                                                                       )
        else:
            raise Exception("Apart from Newton-Raphson, other strategy_type are not available.")
        
        return solver