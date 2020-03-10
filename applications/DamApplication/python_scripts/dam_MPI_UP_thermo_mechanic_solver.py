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

from KratosMultiphysics.DamApplication import dam_MPI_thermo_mechanic_solver


def CreateSolver(main_model_part, custom_settings):

    return DamMPIUPThermoMechanicSolver(main_model_part, custom_settings)


class DamMPIUPThermoMechanicSolver(dam_MPI_thermo_mechanic_solver.DamMPIThermoMechanicSolver):

    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "dam_MPI_UP_thermo_mechanic_solver",
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
                "solution_type": "Dynamic",
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

        print("Construction of Dam_MPI_UP_ThermoMechanicSolver finished")

    def AddVariables(self):

        super(DamMPIUPThermoMechanicSolver, self).AddVariables()

        ## Fluid variables
        # Add pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        # Add Dynamic pressure variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Dt_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Dt2_PRESSURE)

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            ## Solid dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X,KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y,KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z,KratosMultiphysics.REACTION_Z)
            ## Fluid dofs
            node.AddDof(KratosMultiphysics.PRESSURE)
            ## Thermal dofs
            node.AddDof(KratosMultiphysics.TEMPERATURE)
            # adding VELOCITY as dofs
            node.AddDof(KratosMultiphysics.VELOCITY_X)
            node.AddDof(KratosMultiphysics.VELOCITY_Y)
            node.AddDof(KratosMultiphysics.VELOCITY_Z)
            # adding ACCELERATION as dofs
            node.AddDof(KratosMultiphysics.ACCELERATION_X)
            node.AddDof(KratosMultiphysics.ACCELERATION_Y)
            node.AddDof(KratosMultiphysics.ACCELERATION_Z)

    #### Specific internal functions ####

    def _ConstructScheme(self, scheme_type, solution_type):

        rayleigh_m = self.settings["mechanical_solver_settings"]["rayleigh_m"].GetDouble()
        rayleigh_k = self.settings["mechanical_solver_settings"]["rayleigh_k"].GetDouble()

        beta=0.25
        gamma=0.5
        scheme = KratosDam.TrilinosDamUPScheme(beta,gamma,rayleigh_m,rayleigh_k)

        return scheme
