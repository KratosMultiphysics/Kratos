import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

from KratosMultiphysics.PfemFluidDynamicsApplication import pfem_fluid_solver as BaseSolver


def CreateSolver(model, parameters):
    return PfemFluidThreeStepSolver(model, parameters)

class PfemFluidThreeStepSolver(BaseSolver.PfemFluidSolver):

    def __init__(self, model, parameters):

        super(PfemFluidThreeStepSolver,self).__init__(model, parameters)

        print("Construction of 3-step Pfem Fluid Three Step Solver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
             "solver_type": "pfem_fluid_three_step_solver",
            "model_part_name": "PfemFluidModelPart",
            "physics_type"   : "fluid",
            "domain_size": 2,
            "time_stepping"               : {
                "automatic_time_step" : false,
                "time_step"           : 0.001
            },
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings"           : {
                "materials_filename" : "unknown_name"
            },
            "buffer_size": 3,
            "echo_level": 1,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": true,
            "move_mesh_flag": true,
            "dofs"                : [],
            "stabilization_factor": 1.0,
            "line_search": false,
            "compute_contact_forces": false,
            "block_builder": false,
            "component_wise": false,
            "predictor_corrector": true,
            "time_order": 2,
            "maximum_velocity_iterations": 1,
            "maximum_pressure_iterations": 7,
            "velocity_tolerance": 1e-5,
            "pressure_tolerance": 1e-5,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "provide_coordinates"            : false,
                "scaling"                        : false,
                "smoother_type"                  : "damped_jacobi",
                "krylov_type"                    : "cg",
                "coarsening_type"                : "aggregation",
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "bicgstab",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "preconditioner_type"            : "none",
                "scaling"                        : false
            },
            "solving_strategy_settings":{
               "time_step_prediction_level": 0,
               "max_delta_time": 1.0e-5,
               "fraction_delta_time": 0.9,
               "rayleigh_damping": false,
               "rayleigh_alpha": 0.0,
               "rayleigh_beta" : 0.0
            },
        "bodies_list": [],
        "problem_domain_sub_model_part_list": [],
        "constitutive_laws_list": [],
        "processes_sub_model_part_list": [],
        "constraints_process_list": [],
        "loads_process_list"       : [],
        "output_process_list"      : [],
        "output_configuration"     : {},
        "problem_process_list"     : [],
        "processes"                : {},
        "output_processes"         : {},
        "check_process_list": []
        }""")
        this_defaults.AddMissingParameters(super(PfemFluidThreeStepSolver, cls).GetDefaultParameters())
        return this_defaults


    def Initialize(self):

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        self.fluid_solver = KratosPfemFluid.ThreeStepVPStrategy(self.computing_model_part,
                                                                self.velocity_linear_solver,
                                                                self.pressure_linear_solver,
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["velocity_tolerance"].GetDouble(),
                                                                self.settings["pressure_tolerance"].GetDouble(),
                                                                self.settings["maximum_pressure_iterations"].GetInt(),
                                                                self.settings["time_order"].GetInt(),
                                                                self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION])


        # Set echo_level
        echo_level = self.settings["echo_level"].GetInt()
        self.fluid_solver.SetEchoLevel(echo_level)

        # Self initialize strategy
        self.fluid_solver.Initialize()

        # Check if everything is assigned correctly
        self.fluid_solver.Check() #TODO: This must be done in the Check function

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FRACT_VEL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VOLUME)
