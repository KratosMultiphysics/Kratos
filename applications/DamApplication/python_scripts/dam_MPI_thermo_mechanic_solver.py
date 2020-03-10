from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam
from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

from KratosMultiphysics.DamApplication import dam_thermo_mechanic_solver
from KratosMultiphysics.mpi import distributed_import_model_part_utility


def CreateSolver(main_model_part, custom_settings):

    return DamMPIThermoMechanicSolver(main_model_part, custom_settings)


class DamMPIThermoMechanicSolver(dam_thermo_mechanic_solver.DamThermoMechanicSolver):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "dam_MPI_thermo_mechanic_solver",
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "echo_level": 0,
            "buffer_size": 2,
            "reference_temperature" : 10.0,
            "processes_sub_model_part_list": [""],
            "thermal_solver_settings":{
                "echo_level": 0,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": false,
                "move_mesh_flag": true,
                "compute_norm_dx_flag": false,
                "theta_scheme": 1.0,
                "block_builder": true,
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
                "thermal_loads_sub_model_part_list": []
            },
            "mechanical_solver_settings":{
                "echo_level": 0,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": false,
                "move_mesh_flag": true,
                "solution_type": "Quasi-Static",
                "scheme_type": "Newmark",
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
                "body_domain_sub_model_part_list": [],
                "mechanical_loads_sub_model_part_list": [],
                "loads_sub_model_part_list": [],
                "loads_variable_list": []
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Construct the linear solver
        self.thermal_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["thermal_solver_settings"]["linear_solver_settings"])
        self.mechanical_linear_solver = trilinos_linear_solver_factory.ConstructSolver(self.settings["mechanical_solver_settings"]["linear_solver_settings"])

        print("Construction of MPI DamThermoMechanicSolver finished")

    def AddVariables(self):

        super(DamMPIThermoMechanicSolver, self).AddVariables()

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)

    def ImportModelPart(self):

        # Construct the import model part utility
        ModelPartImporter = distributed_import_model_part_utility.DistributedImportModelPartUtility(self.main_model_part, self.settings)

        # Execute the Metis partitioning and reading
        ModelPartImporter.ExecutePartitioningAndReading()

        # Create computing_model_part, set constitutive law and buffer size
        self._ExecuteAfterReading()

        # Construct the communicators
        ModelPartImporter.CreateCommunicators()

    def Initialize(self):

        # Construct the communicator
        self.EpetraCommunicator = TrilinosApplication.CreateCommunicator()

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.REFERENCE_TEMPERATURE, self.settings["reference_temperature"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosConvDiff.THETA, self.settings["thermal_solver_settings"]["theta_scheme"].GetDouble())

        # Get the computing model parts
        self.thermal_computing_model_part = self.main_model_part.GetSubModelPart(self.thermal_model_part_name)
        self.mechanical_computing_model_part = self.GetComputingModelPart()

        # Builder and solver creation
        thermal_builder_and_solver = self._ConstructBuilderAndSolver(self.settings["thermal_solver_settings"]["block_builder"].GetBool(),self.thermal_linear_solver)
        mechanical_builder_and_solver = self._ConstructBuilderAndSolver(self.settings["mechanical_solver_settings"]["block_builder"].GetBool(),self.mechanical_linear_solver)

        # Solution scheme creation
        thermal_scheme = TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme()
        mechanical_scheme = self._ConstructScheme(self.settings["mechanical_solver_settings"]["scheme_type"].GetString(),
                                         self.settings["mechanical_solver_settings"]["solution_type"].GetString())

        # Get the convergence criterion
        convergence_criterion = self._ConstructConvergenceCriterion(self.settings["mechanical_solver_settings"]["convergence_criterion"].GetString())

        # Solver creation (Note: this could be TrilinosResidualBasedLinearStrategy, but there is no such strategy)
        self.Thermal_Solver = TrilinosApplication.TrilinosNewtonRaphsonStrategy(self.thermal_computing_model_part,
                                                                       thermal_scheme,
                                                                       self.thermal_linear_solver,
                                                                       convergence_criterion,
                                                                       thermal_builder_and_solver,
                                                                       self.settings["mechanical_solver_settings"]["max_iteration"].GetInt(),
                                                                       self.settings["thermal_solver_settings"]["compute_reactions"].GetBool(),
                                                                       self.settings["thermal_solver_settings"]["reform_dofs_at_each_step"].GetBool(),
                                                                       self.settings["thermal_solver_settings"]["move_mesh_flag"].GetBool()
                                                                       )
        self.Mechanical_Solver = self._ConstructSolver(mechanical_builder_and_solver,
                                                       mechanical_scheme,
                                                       convergence_criterion,
                                                       self.settings["mechanical_solver_settings"]["strategy_type"].GetString())

        # Set echo_level
        self.Thermal_Solver.SetEchoLevel(self.settings["thermal_solver_settings"]["echo_level"].GetInt())
        self.Mechanical_Solver.SetEchoLevel(self.settings["mechanical_solver_settings"]["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Thermal_Solver.Check()
        self.Mechanical_Solver.Check()

        print ("Initialization MPI DamThermoMechanicSolver finished")

    #### Specific internal functions ####

    def _ConstructBuilderAndSolver(self, block_builder, linear_solver):

        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            guess_row_size = 15
        else:
            guess_row_size = 45

        if(block_builder):
            builder_and_solver = TrilinosApplication.TrilinosBlockBuilderAndSolver(self.EpetraCommunicator,
                                                                                   guess_row_size,
                                                                                   linear_solver
                                                                                   )
        else:
            builder_and_solver = TrilinosApplication.TrilinosEliminationBuilderAndSolver(self.EpetraCommunicator,
                                                                                         guess_row_size,
                                                                                         linear_solver
                                                                                         )

        return builder_and_solver

    def _ConstructScheme(self, scheme_type, solution_type):

        rayleigh_m = self.settings["mechanical_solver_settings"]["rayleigh_m"].GetDouble()
        rayleigh_k = self.settings["mechanical_solver_settings"]["rayleigh_k"].GetDouble()
        self.main_model_part.ProcessInfo.SetValue(KratosStructural.RAYLEIGH_ALPHA, rayleigh_m)
        self.main_model_part.ProcessInfo.SetValue(KratosStructural.RAYLEIGH_BETA, rayleigh_k)
        if(solution_type == "Quasi-Static"):
            if(rayleigh_m<1.0e-20 and rayleigh_k<1.0e-20):
                scheme =  TrilinosApplication.TrilinosResidualBasedIncrementalUpdateStaticScheme()
            else:
                scheme =  KratosDam.TrilinosIncrementalUpdateStaticDampedScheme()
        else:
            if(scheme_type == "Newmark"):
                damp_factor_m = 0.0
            else:
                damp_factor_m = -0.01
            scheme = TrilinosApplication.TrilinosResidualBasedBossakDisplacementScheme(damp_factor_m)

        return scheme

    def _ConstructConvergenceCriterion(self, convergence_criterion):

        D_RT = self.settings["mechanical_solver_settings"]["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["mechanical_solver_settings"]["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["mechanical_solver_settings"]["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["mechanical_solver_settings"]["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["mechanical_solver_settings"]["echo_level"].GetInt()

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

        max_iters = self.settings["mechanical_solver_settings"]["max_iteration"].GetInt()
        compute_reactions = self.settings["mechanical_solver_settings"]["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["mechanical_solver_settings"]["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["mechanical_solver_settings"]["move_mesh_flag"].GetBool()

        if strategy_type == "Newton-Raphson":
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
            solver = TrilinosApplication.TrilinosNewtonRaphsonStrategy(self.mechanical_computing_model_part,
                                                                       scheme,
                                                                       self.mechanical_linear_solver,
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
