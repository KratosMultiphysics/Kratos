# KratosTopologyOptimizationApplication
#
# License:         BSD License
#                  license: TopologyOptimizationApplication/license.txt
#
# Main authors:    Philipp Hofer, https://github.com/PhiHo-eng
#                  Erich Wehrle, https://github.com/e-dub
#                  based on original file from
#                  Baumgärtner Daniel, https://github.com/dbaumgaertner
#                  Octaviano Malfavón Farías
#                  Eric Gonzales
#
import KratosMultiphysics as km
import KratosMultiphysics.TopologyOptimizationApplication as kto
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.vtk_output_process import VtkOutputProcess
import time


def ConstructOptimizer(opt_model_part, config, analyzer):
    optimizer = SIMPMethod(opt_model_part, config, analyzer)
    return optimizer


class SIMPMethod:
    def __init__(self, opt_model_part, config, analyzer):
        # Set Topology Optimization configurations
        self.config = config

        # For GID output
        self.gid_io = GiDOutputProcess(
            opt_model_part,
            'Topology_Optimization_Results',
            km.Parameters(
                """
                {
                "result_file_configuration": {
                    "gidpost_flags": {
                        "GiDPostMode": "GiD_PostBinary",
                        "WriteDeformedMeshFlag": "WriteUndeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "file_label": "time",
                    "output_control_type": "step",
                    "output_interval": 1.0,
                    "body_output": true,
                    "node_output": false,
                    "skin_output": false,
                    "plane_output": [],
                    "nodal_results": ["DISPLACEMENT","REACTION"],
                    "nodal_nonhistorical_results": [],
                    "nodal_flags_results": [],
                    "gauss_point_results": ["X_PHYS","VON_MISES_STRESS"],
                    "additional_list_files": []
                    }
                }
                """
            ),
        )
        vtk_parameters = km.Parameters(
            """
            {
            "model_part_name": "",
            "output_control_type": "step",
            "output_interval": 1,
            "file_format": "ascii",
            "output_precision": 7,
            "output_sub_model_parts": false,
            "output_path": "vtk_output",
            "save_output_files_in_folder": true,
            "nodal_solution_step_data_variables" : ["DISPLACEMENT","REACTION"],
            "nodal_data_value_variables": [],
            "element_data_value_variables": ["X_PHYS","VON_MISES_STRESS"],
            "condition_data_value_variables": [],
            "gauss_point_variables_extrapolated_to_nodes": []
            }
            """
        )
        vtk_parameters['model_part_name'].SetString(opt_model_part.Name)
        self.vtk_io = VtkOutputProcess(
            opt_model_part.GetModel(), vtk_parameters
        )

        # Set analyzer
        self.analyzer = analyzer

        # Set response functions
        self.objectives = config['objectives'].items()
        self.constraints = config['constraints'].items()

        km.Logger.Print(
            '\n::[Initializing Topology Optimization Application]::'
        )
        km.Logger.Print('  The following objectives are defined:')
        for func_id, obj_settings in config['objectives'].items():
            obj_settings['grad'].GetString()
            km.Logger.Print(
                '   ',
                func_id,
                "-> 'grad' : ",
                obj_settings['grad'].GetString(),
                '\n',
            )

        km.Logger.Print('  The following constraints are defined:')
        for func_id, const_settings in config['constraints'].items():
            const_settings['grad'].GetString()
            km.Logger.Print(
                '   ',
                func_id,
                "-> 'type' : ",
                const_settings['type'].GetString(),
                ", 'grad' : ",
                const_settings['grad'].GetString(),
                '\n',
            )

        # Create controller object
        self.controller = Controller(config)

        # Model parameters
        self.opt_model_part = opt_model_part

        # Initialize element variables
        for element_i in opt_model_part.Elements:
            element_i.SetValue(kto.PENAL, config['penalty'].GetInt())
            element_i.SetValue(
                kto.MAT_INTERP, config['material_interpolation'].GetString()
            )
            element_i.SetValue(
                kto.X_PHYS, config['initial_volume_fraction'].GetDouble()
            )
            element_i.SetValue(
                kto.X_PHYS_OLD, config['initial_volume_fraction'].GetDouble()
            )
            element_i.SetValue(
                km.YOUNG_MODULUS,
                opt_model_part.GetProperties()[
                    config['simp_property'].GetInt()
                ].GetValue(km.YOUNG_MODULUS),
            )
            elemental_volume = element_i.GetGeometry().DomainSize()
            element_i.SetValue(kto.INITIAL_ELEMENT_SIZE, elemental_volume)

        # Only happens if continuation strategy is activated (Initialization of penalty factor)
        if self.config['continuation_strategy'] == 1:
            for element_i in self.opt_model_part.Elements:
                element_i.SetValue(kto.PENAL, 1)

        # Add toolbox for topology filtering utilities
        self.filter_utils = kto.TopologyFilteringUtilities(
            opt_model_part,
            self.config['filter_radius'].GetDouble(),
            self.config['max_elements_in_filter_radius'].GetInt(),
        )

        # Add toolbox for topology updating utilities
        self.design_update_utils = kto.TopologyUpdatingUtilities(
            opt_model_part
        )

        # Add toolbox for I/O
        self.io_utils = kto.IOUtilities()

    def optimize(self):
        km.Logger.Print('\n>' + '-' * 45)
        km.Logger.Print('> Starting topology optimization')
        km.Logger.Print('>' + '-' * 45 + '\n')

        # Start timer and assign to object such that total time of opt may be measured at each step
        self.opt_start_time = time.time()

        # Initialize the design output in GiD format and print initial 0 state
        self.gid_io.ExecuteInitialize()
        self.gid_io.ExecuteBeforeSolutionLoop()
        self.vtk_io.ExecuteInitialize()
        self.vtk_io.ExecuteBeforeSolutionLoop()

        # Call for the specified optimization algorithm
        if self.config['optimization_algorithm'].GetString() == 'oc_algorithm':
            self.start_oc_algorithm()

        else:
            raise TypeError(
                'Specified optimization_algorithm not implemented!'
            )

        # Finalize the design output in GiD and vtk formats
        self.gid_io.PrintOutput()
        self.gid_io.ExecuteFinalizeSolutionStep()
        self.gid_io.ExecuteFinalize()

        self.vtk_io.PrintOutput()
        self.vtk_io.ExecuteFinalizeSolutionStep()
        self.vtk_io.ExecuteFinalize()

        # Stop timer
        opt_end_time = time.time()

        km.Logger.Print('\n>' + '-' * 45)
        km.Logger.Print('> Topology optimization complete')
        # km.Logger.Print('> Iteration: ', opt_itr)
        km.Logger.Print(
            '> Duration: ', round(opt_end_time - self.opt_start_time, 1), ' s!'
        )
        km.Logger.Print('>' + '-' * 45 + '\n')

    # Topology Optimization with OC
    def start_oc_algorithm(self):
        # Get Id of objective & constraint
        only_F_id = None
        only_C_id = None
        for F_id, empty_id_f in self.objectives:
            only_F_id = F_id
            break
        for C_id, empty_id_c in self.constraints:
            only_C_id = C_id
            break

        # Initialize variables for comparison purposes in Topology Optimization Tool
        pmax = self.config[
            'penalty'
        ].GetInt()  # Maximum penalty value used for continuation strategy
        Obj_Function = None
        Obj_Function_old = None
        Obj_Function_initial = None
        Obj_Function_relative_change = None
        Obj_Function_absolute_change = None

        # Print the Topology Optimization Settings that will be used in the program
        km.Logger.Print('  Topology optimization settings:')
        km.Logger.Print(
            '  Algorithm:         ',
            self.config['optimization_algorithm'].GetString(),
        )
        km.Logger.Print(
            '  Material interp.:  ',
            self.config['material_interpolation'].GetString(),
        )
        km.Logger.Print(
            '  Penalty factor:    ', self.config['penalty'].GetInt()
        )
        km.Logger.Print(
            '  Emin:              ',
            self.opt_model_part.GetProperties()[
                self.config['simp_property'].GetInt()
            ].GetValue(kto.YOUNGS_MODULUS_MIN),
        )
        km.Logger.Print(
            '  Filter radius:     ', self.config['filter_radius'].GetDouble()
        )
        km.Logger.Print(
            '  Relative tolerance:',
            self.config['relative_tolerance'].GetDouble(),
        )
        km.Logger.Print(
            '  Volume fraction:   ',
            self.config['initial_volume_fraction'].GetDouble(),
        )
        km.Logger.Print(
            '  Maximum iterations:', self.config['max_opt_iterations'].GetInt()
        )

        if (
            self.config['restart_write_frequency'].GetInt()
            < self.config['max_opt_iterations'].GetInt()
        ):
            if self.config['restart_write_frequency'].GetInt() == 1:
                km.Logger.Print('  Restart file every iteration')
            elif self.config['restart_write_frequency'].GetInt() > 1:
                km.Logger.Print(
                    '  Restart file:       every',
                    self.config['restart_write_frequency'].GetInt(),
                    'iterations',
                )
            else:
                km.Logger.Print(
                    '  No restart file created during the simulation'
                )
        else:
            km.Logger.Print('  No restart file created during the simulation')

        # Start optimization loop
        for opt_itr in range(
            1, self.config['max_opt_iterations'].GetInt() + 1
        ):
            # Some output
            km.Logger.Print('\n>' + '-' * 45)
            km.Logger.Print('> Topology optimization iteration ', opt_itr)
            km.Logger.Print('>' + '-' * 45 + '\n')

            # Start measuring time needed for current optimization step
            start_time = time.time()

            # Initialize response container
            response = self.controller.create_response_container()

            # Set controller to evaluate objective & constraint
            self.controller.initialize_controls()
            self.controller.get_controls()[only_F_id]['calc_func'] = 1
            self.controller.get_controls()[only_C_id]['calc_func'] = 1

            # Set to evaluate objective & constraint gradient if provided
            if (
                self.config['objectives'][only_F_id]['grad'].GetString()
                == 'provided'
            ):
                self.controller.get_controls()[only_F_id]['calc_grad'] = 1
            if (
                self.config['constraints'][only_C_id]['grad'].GetString()
                == 'provided'
            ):
                self.controller.get_controls()[only_C_id]['calc_grad'] = 1

            # RUN FEM: Call analyzer with current X to compute response (global_strain_energy, dcdx)
            self.analyzer(self.controller.get_controls(), response, opt_itr)

            # Filter sensitivities
            km.Logger.Print('\n[TopOpt]:   ::[Filter Sensitivities]::')
            self.filter_utils.ApplyFilterSensitivity(
                self.config['filter_type'].GetString(),
                self.config['filter_kernel'].GetString(),
            )

            # Update design variables ( densities )  --> new X by:
            km.Logger.Print('\n[TopOpt]    ::[Update Densities with OC]::')
            self.design_update_utils.UpdateDensitiesUsingOCMethod(
                self.config['optimization_algorithm'].GetString(),
                self.config['initial_volume_fraction'].GetDouble(),
                self.config['grey_scale_filter'].GetInt(),
                opt_itr,
                self.config['q_max'].GetDouble(),
            )

            if self.config['density_filter'].GetString() == 'density':
                km.Logger.Print('\n[TopOpt]   ::[Filter Densities]::')
                self.filter_utils.ApplyFilterDensity(
                    self.config['density_filter'].GetString(),
                    self.config['filter_kernel'].GetString(),
                )

            # Print of results
            km.Logger.Print('\n[TopOpt]:   ::[RESULTS]::')
            Obj_Function = response[only_F_id]['func']
            C_Function = response[only_C_id]['func']

            km.Logger.Print(
                '  Objective value:              ',
                '{:.6f}'.format(Obj_Function),
            )
            km.Logger.Print(
                '  Constraint value:             ',
                '{:.6f}'.format(C_Function),
            )

            if opt_itr == 1:
                Obj_Function_initial = Obj_Function

            if opt_itr > 1:
                Obj_Function_relative_change = (
                    Obj_Function - Obj_Function_old
                ) / Obj_Function_initial
                km.Logger.Print(
                    '  Relative change in objective: ',
                    '{:.9f}'.format(Obj_Function_relative_change),
                )

                Obj_Function_absolute_change = (
                    Obj_Function - Obj_Function_initial
                ) / Obj_Function_initial
                km.Logger.Print(
                    '  Absolute change in objective: ',
                    '{:.9f}'.format(Obj_Function_absolute_change),
                )

            Obj_Function_old = Obj_Function

            # Write design in GiD format
            self.gid_io.ExecuteFinalizeSolutionStep()

            # Continuation Strategy
            if self.config['continuation_strategy'].GetInt() == 1:
                km.Logger.Print('  Continuation strategy:         active')
                if opt_itr < 20:
                    for element_i in self.opt_model_part.Elements:
                        element_i.SetValue(kto.PENAL, 1)
                else:
                    for element_i in self.opt_model_part.Elements:
                        element_i.SetValue(
                            kto.PENAL,
                            min(pmax, 1.02 * element_i.GetValue(kto.PENAL)),
                        )
            else:
                km.Logger.Print('  Continuation strategy:         inactive')

            # Write restart file every selected number of iterations
            restart_filename = (
                self.config['restart_output_file']
                .GetString()
                .replace('.mdpa', '_' + str(opt_itr) + '.mdpa')
            )
            if self.config['restart_write_frequency'].GetInt() > 0:
                if (
                    opt_itr % self.config['restart_write_frequency'].GetInt()
                    == False
                ):
                    km.Logger.Print('\n::[Restart File]::')
                    km.Logger.Print('  Saving file at iteration', opt_itr)
                    self.io_utils.SaveOptimizationResults(
                        self.config['restart_input_file'].GetString(),
                        self.opt_model_part,
                        restart_filename,
                    )

            # Check convergence
            if opt_itr > 1:
                # Check if maximum iterations were reached
                if opt_itr == self.config['max_opt_iterations'].GetInt():
                    end_time = time.time()
                    km.Logger.Print(
                        '  Iteration time:               ',
                        round(end_time - start_time, 1),
                        's',
                    )
                    km.Logger.Print(
                        '  Elapsed time:                 ',
                        round(end_time - self.opt_start_time, 1),
                        's',
                    )
                    km.Logger.Print('\n  Maximum iterations reached!')
                    self.io_utils.SaveOptimizationResults(
                        self.config['restart_input_file'].GetString(),
                        self.opt_model_part,
                        restart_filename,
                    )
                    break

                # Check for relative tolerance
                if (
                    abs(Obj_Function_relative_change)
                    < self.config['relative_tolerance'].GetDouble()
                ):
                    end_time = time.time()
                    km.Logger.Print(
                        '  Iteration time:               ',
                        round(end_time - start_time, 1),
                        's',
                    )
                    km.Logger.Print(
                        '  Elapsed time:                 ',
                        round(end_time - self.opt_start_time, 1),
                        's',
                    )
                    km.Logger.Print(
                        '\n  Converged to relative change in objective: ',
                        self.config['relative_tolerance'].GetDouble(),
                    )
                    self.io_utils.SaveOptimizationResults(
                        self.config['restart_input_file'].GetString(),
                        self.opt_model_part,
                        restart_filename,
                    )

                    break

            # Set X_PHYS_OLD to update the value for the next simulation's "change percentage"
            for element_i in self.opt_model_part.Elements:
                element_i.SetValue(
                    kto.X_PHYS_OLD, element_i.GetValue(kto.X_PHYS)
                )

            # Take time needed for current optimization step
            end_time = time.time()
            km.Logger.Print(
                '  Iteration time:               ',
                round(end_time - start_time, 1),
                's',
            )
            km.Logger.Print(
                '  Elapsed time:                 ',
                round(end_time - self.opt_start_time, 1),
                's',
            )


class Controller:
    def __init__(self, config):
        # Create and initialize controller
        self.controls = {}
        for func_id, empty in config['objectives'].items():
            km.Logger.Print(' Print the controller input   ', func_id, '\n')
            self.controls[func_id] = {'calc_func': 0, 'calc_grad': 0}
            km.Logger.Print(
                ' Print the controller output  ', self.controls[func_id], '\n'
            )

        for func_id, empty in config['constraints'].items():
            self.controls[func_id] = {'calc_func': 0, 'calc_grad': 0}

        # Initialize response container to provide storage for any response
        self.response_container = {}
        for func_id, empty in config['objectives'].items():
            self.response_container[func_id] = {'func': None, 'grad': None}
        for func_id, empty in config['constraints'].items():
            self.response_container[func_id] = {'func': None, 'grad': None}

    def initialize_controls(self):
        # Sets
        for func_id in self.controls:
            self.controls[func_id] = {'calc_func': 0, 'calc_grad': 0}

    def get_controls(self):
        return self.controls

    def create_response_container(self):
        # Create and initialize container to store any response defined
        for func_id in self.response_container:
            self.response_container[func_id] = {'func': None, 'grad': None}

        # Return container
        return self.response_container
