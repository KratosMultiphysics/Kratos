from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver as BaseSolver

def CreateSolver(custom_settings):
    return SegregatedImplicitSolver(custom_settings)

class SegregatedImplicitSolver(BaseSolver.MechanicalSolver):

    def __init__(self, custom_settings):

        # Set defaults and validate custom settings.
        segregated_solver_settings = KratosMultiphysics.Parameters("""
        {     
           "time_integration_settings":{
                "buffer_size": 3
            },
            "solving_strategy_settings":{
                "time_order": 2,
                "maximum_velocity_iterations": 1,
                "maximum_pressure_iterations": 7,
                "velocity_tolerance": 1.0e-4,
                "pressure_tolerance": 1.0e-9,
                "clear_storage": false,
                "move_mesh_flag": true,        
                "reform_dofs_at_each_step": false
            },
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "AMGCL",
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
                "solver_type"                    : "BICGSTABSolver",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "preconditioner_type"            : "None",
                "scaling"                        : false
            }
        } 
        """)

        # Validate and transfer settings
        if( custom_settings.Has("solving_strategy_settings") ):
            self._validate_and_transfer_matching_settings(custom_settings["solving_strategy_settings"], segregated_solver_settings["solving_strategy_settings"])
           
        if( custom_settings.Has("time_integration_settings") ):
            self._validate_and_transfer_matching_settings(custom_settings["time_integration_settings"], segregated_solver_settings["time_integration_settings"])
        else:
            custom_settings.AddValue("time_integration_settings", segregated_solver_settings["time_integration_settings"])
      
        self.segregated_solver_settings = segregated_solver_settings["solving_strategy_settings"]

        self.pressure_solver_settings = segregated_solver_settings["pressure_linear_solver_settings"]
        self.velocity_solver_settings = segregated_solver_settings["velocity_linear_solver_settings"]
        
        # Construct the base solver.
        super(SegregatedImplicitSolver, self).__init__(custom_settings)
            
        print("::[Segregated_Scheme]:: Scheme Ready")


    def _create_solution_scheme(self):
        self._solution_scheme = None
        return self._solution_scheme
        
    def _create_mechanical_solver(self):
        mechanical_solver = self._create_two_step_strategy()
        return mechanical_solver


    def _get_velocity_linear_solver(self):
        if not hasattr(self, '_velocity_linear_solver'):
            self._velocity_linear_solver = self._create_velocity_linear_solver()
        return self._velocity_linear_solver

    def _get_pressure_linear_solver(self):
        if not hasattr(self, '_pressure_linear_solver'):
            self._pressure_linear_solver = self._create_pressure_linear_solver()
        return self._pressure_linear_solver
    
    def _create_velocity_linear_solver(self):
        import linear_solver_factory
        velocity_linear_solver = linear_solver_factory.ConstructSolver(self.velocity_solver_settings)
        return velocity_linear_solver

    def _create_pressure_linear_solver(self):
        import linear_solver_factory
        pressure_linear_solver = linear_solver_factory.ConstructSolver(self.pressure_solver_settings)
        return pressure_linear_solver

    def _create_two_step_strategy(self):
        mechanical_scheme = self._get_solution_scheme()
        velocity_linear_solver = self._get_velocity_linear_solver()
        pressure_linear_solver = self._get_pressure_linear_solver()
        time_order = 2
        mechanical_convergence_criterion = self._get_convergence_criterion()
        builder_and_solver = self._get_builder_and_solver()
        return KratosPfemFluid.TwoStepVPStrategy(self.model_part,
                                                 velocity_linear_solver,
                                                 pressure_linear_solver,
                                                 self.solving_strategy_settings["reform_dofs_at_each_step"].GetBool(),
                                                 self.segregated_solver_settings["velocity_tolerance"].GetDouble(),
                                                 self.segregated_solver_settings["pressure_tolerance"].GetDouble(),
                                                 self.segregated_solver_settings["maximum_pressure_iterations"].GetInt(),
                                                 self.segregated_solver_settings["time_order"].GetInt(),
                                                 self.process_info[KratosMultiphysics.SPACE_DIMENSION])
    def GetVariables(self):
        
        nodal_variables = super(SegregatedImplicitSolver, self).GetVariables()
        nodal_variables = nodal_variables + ['MESH_VELOCITY','NODAL_AREA']
        nodal_variables = nodal_variables + ['FREESURFACE','PRESSURE_VELOCITY','PRESSURE_ACCELERATION']
        nodal_variables = nodal_variables + ['DENSITY','BULK_MODULUS','VISCOSITY','POISSON_RATIO','YOUNG_MODULUS'] #material
        nodal_variables = nodal_variables + ['FLOW_INDEX','YIELD_SHEAR','ADAPTIVE_EXPONENT'] #non-newtonian

        return nodal_variables
    

    def InitializeSolutionStep(self):
        adaptive_time_interval = KratosPfemFluid.AdaptiveTimeIntervalProcess(self.main_model_part,self.echo_level)
        adaptive_time_interval.Execute()

        unactive_peak_elements = False
        unactive_sliver_elements = True
        set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        set_active_flag.Execute()

        #split_elements = KratosPfemFluid.SplitElementsProcess(self.main_model_part,self.echo_level)
        #split_elements.ExecuteInitialize()

    def Predict(self):
        pass

    def FinalizeSolutionStep(self):
        self._get_mechanical_solver().FinalizeSolutionStep()  

        unactive_peak_elements = False
        unactive_sliver_elements = True
        set_active_flag = KratosPfemFluid.SetActiveFlagProcess(self.main_model_part,unactive_peak_elements,unactive_sliver_elements,self.echo_level)
        set_active_flag.ExecuteFinalize()

        #split_elements = KratosPfemFluid.SplitElementsProcess(self.main_model_part,self.echo_level)
        #split_elements.ExecuteFinalize()


