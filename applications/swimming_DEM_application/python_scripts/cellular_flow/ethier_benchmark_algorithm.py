from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_algorithm
import swimming_DEM_procedures as SDP
import math
import numpy as np
BaseAlgorithm = swimming_DEM_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)

    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)
        Add = self.pp.CFD_DEM.AddEmptyValue
        Add("pressure_grad_recovery_type")
        Add("size_parameter").SetDouble(1)
        Add("store_full_gradient_option").SetBool(True)
        Add("regular_mesh_option").SetBool(True)
        Add("mesh_tag").SetString('')
        Add("print_VECTORIAL_ERROR_option").SetBool(True)

    def SetCustomBetaParameters(self, custom_parameters):
        BaseAlgorithm.SetCustomBetaParameters(self, custom_parameters)
        self.pp.CFD_DEM.size_parameter = self.pp.CFD_DEM["size_parameter"].GetDouble()
        # Creating a code for the used input variables
        self.run_code = '_ndiv_' \
                      + str(self.pp.CFD_DEM["size_parameter"].GetDouble()) \
                      + '_mat_deriv_type_' \
                      + str(self.pp.CFD_DEM["material_acceleration_calculation_type"].GetInt()) \
                      + '_lapl_type_' \
                      + str(self.pp.CFD_DEM["laplacian_calculation_type"].GetInt())

    def ReadFluidModelParts(self):
        problem_name = self.pp.problem_name.replace('Fluid', '')
        is_regular_mesh = self.pp.CFD_DEM["regular_mesh_option"].GetBool()
        tag =  self.pp.CFD_DEM["mesh_tag"].GetString()

        if is_regular_mesh:
            mdpa_name = problem_name + '_ndiv_' + str(int(self.pp.CFD_DEM.size_parameter)) + tag + 'Fluid'
        else:
            mdpa_name = problem_name + '_h_' + str(self.pp.CFD_DEM.size_parameter) + 'Fluid'

        model_part_io_fluid = ModelPartIO(mdpa_name)
        model_part_io_fluid.ReadModelPart(self.fluid_solution.fluid_model_part)

    def AddExtraVariables(self, run_code = ''):
        BaseAlgorithm.AddExtraVariables(self, self.run_code)

    def GetParticlesResultsCounter(self):
        return SDP.Counter()

    def GetPrintCounterUpdatedDEM(self):
        return SDP.Counter(is_dead=True)

    def GetFieldUtility(self):
        a = math.pi / 4
        d = math.pi / 2 * 0

        self.flow_field = EthierFlowField(a, d)
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def GetRecoveryCounter(self):
        return SDP.Counter(1, 1, self.pp.CFD_DEM["coupling_level_type"].GetInt() or self.pp.CFD_DEM.print_PRESSURE_GRADIENT_option)

    def GetRunCode(self):
        return self.run_code

    def PerformZeroStepInitializations(self):
        self.mat_deriv_errors = []
        self.laplacian_errors = []
        self.current_mat_deriv_errors = np.zeros(2)
        self.current_laplacian_errors = np.zeros(2)

        for node in self.fluid_model_part.Nodes:
            vel= Vector(3)
            coor = Vector([node.X, node.Y, node.Z])
            self.flow_field.Evaluate(0.0, coor, vel, 0)
            node.SetSolutionStepValue(VELOCITY, vel)

    def FluidSolve(self, time = 'None', solve_system = True):

        for node in self.fluid_model_part.Nodes:
            vel= Vector(3)
            coor = Vector([node.X, node.Y, node.Z])
            self.flow_field.Evaluate(time, coor, vel, 0)
            node.SetSolutionStepValue(VELOCITY, vel)

        self.CalculateRecoveryErrors(time)

    def CalculateRecoveryErrors(self, time):
        L2_norm_mat_deriv = 0.
        L2_norm_mat_deriv_error = 0.
        L2_norm_laplacian = 0.
        L2_norm_laplacian_error = 0.
        max_mat_deriv_error = 0.
        max_laplacian_error = 0.

        calc_mat_deriv = np.zeros(3)
        calc_laplacian = np.zeros(3)
        mat_deriv= Vector(3)
        laplacian= Vector(3)
        total_volume = 0.

        for node in self.fluid_model_part.Nodes:
            nodal_volume = node.GetSolutionStepValue(NODAL_AREA)
            total_volume += nodal_volume
            coor = Vector([node.X, node.Y, node.Z])

            self.flow_field.CalculateMaterialAcceleration(time, coor, mat_deriv, 0)
            self.flow_field.CalculateLaplacian(time, coor, laplacian, 0)
            calc_mat_deriv = node.GetSolutionStepValue(MATERIAL_ACCELERATION)
            calc_laplacian = node.GetSolutionStepValue(VELOCITY_LAPLACIAN)

            module_mat_deriv_squared = sum(x**2 for x in mat_deriv)
            module_laplacian_squared = sum(x**2 for x in laplacian)
            L2_norm_mat_deriv += module_mat_deriv_squared * nodal_volume
            L2_norm_laplacian += module_laplacian_squared * nodal_volume
            diff_mat_deriv = calc_mat_deriv - mat_deriv
            diff_laplacian = calc_laplacian - laplacian
            node.SetSolutionStepValue(VECTORIAL_ERROR, Vector(list(diff_mat_deriv)))
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN, Vector(list(diff_laplacian)))
            module_mat_deriv_error_squared = sum(x**2 for x in diff_mat_deriv)
            module_laplacian_error_squared = sum(x**2 for x in diff_laplacian)
            L2_norm_mat_deriv_error += module_mat_deriv_error_squared * nodal_volume
            L2_norm_laplacian_error += module_laplacian_error_squared * nodal_volume
            max_mat_deriv_error = max(max_mat_deriv_error, module_mat_deriv_error_squared)
            max_laplacian_error = max(max_laplacian_error, module_laplacian_error_squared)

        L2_norm_mat_deriv **= 0.5
        L2_norm_mat_deriv /= total_volume ** 0.5
        L2_norm_mat_deriv_error **= 0.5
        L2_norm_mat_deriv_error /= total_volume ** 0.5
        L2_norm_laplacian **= 0.5
        L2_norm_laplacian /= total_volume ** 0.5
        L2_norm_laplacian_error **= 0.5
        L2_norm_laplacian_error /= total_volume ** 0.5
        max_mat_deriv_error ** 0.5
        max_laplacian_error ** 0.5

        if L2_norm_mat_deriv > 0 and L2_norm_laplacian > 0:
            SDP.MultiplyNodalVariableByFactor(self.fluid_model_part, VECTORIAL_ERROR, 1.0 / L2_norm_mat_deriv)
            SDP.MultiplyNodalVariableByFactor(self.fluid_model_part, VELOCITY_LAPLACIAN, 1.0 / L2_norm_laplacian)
            self.current_mat_deriv_errors[0] = L2_norm_mat_deriv_error / L2_norm_mat_deriv
            self.current_mat_deriv_errors[1] = max_mat_deriv_error / L2_norm_mat_deriv
            self.current_laplacian_errors[0] = L2_norm_laplacian_error / L2_norm_laplacian
            self.current_laplacian_errors[1] = max_laplacian_error / L2_norm_laplacian
            self.mat_deriv_errors.append(self.current_mat_deriv_errors)
            self.laplacian_errors.append(self.current_laplacian_errors)
            #print('mat_deriv: min, max, avg, ', mat_deriv_averager.GetCurrentData())
            #print('laplacian: min, max, avg, ', laplacian_averager.GetCurrentData())
            text_width = 40
            print('\n' + '-.' * text_width)
            print('L2 error for the material derivative'.ljust(text_width), self.current_mat_deriv_errors[0])
            print('max error for the material derivative'.ljust(text_width), self.current_mat_deriv_errors[1])
            print('L2 error for the laplacian'.ljust(text_width), self.current_laplacian_errors[0])
            print('max error for the laplacian'.ljust(text_width), self.current_laplacian_errors[1])
            print('-.' * text_width + '\n')

    def PerformFinalOperations(self, time = None):

        if not os.path.exists('../errors_recorded'):
            os.makedirs('../errors_recorded')

        import h5py

        file_name = self.main_path + '/errors_recorded/recovery_errors.hdf5'
        # with h5py.File(self.file_name, 'r+') as f:
        #     f.create_dataset('material_derivative', shape = self.shape, dtype = np.float32)
        size_parameter = self.pp.CFD_DEM.size_parameter
        is_regular_mesh = self.pp.CFD_DEM["regular_mesh_option"].GetBool()

        with h5py.File(file_name) as f:
            mat_deriv_grp = f.require_group('material derivative')
            mat_deriv_mthd_group = mat_deriv_grp.require_group('method = ' + str(self.pp.CFD_DEM["material_acceleration_calculation_type"].GetInt()))
            laplacian_grp = f.require_group('laplacian')
            laplacian_mthd_group = laplacian_grp.require_group('method = ' + str(self.pp.CFD_DEM["laplacian_calculation_type"].GetInt()))

            if is_regular_mesh:
                mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('regular mesh')
                mesh_grp_laplacian = laplacian_mthd_group.require_group('regular mesh')
                dset_mat_deriv = mesh_grp_mat_deriv.require_dataset('n_div = ' + str(size_parameter),
                                                                    (2, ),
                                                                    dtype=np.float64)
                dset_laplacian = mesh_grp_laplacian.require_dataset('n_div = ' + str(size_parameter),
                                                                    (2,),
                                                                    dtype=np.float64)
            else:
                mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('irregular mesh')
                mesh_grp_laplacian = laplacian_mthd_group.require_group('irregular mesh')
                dset_mat_deriv = mesh_grp_mat_deriv.require_dataset('h = ' + str(size_parameter),
                                                                    (2, ),
                                                                    dtype=np.float64)
                dset_laplacian = mesh_grp_laplacian.require_dataset('h = ' + str(size_parameter),
                                                                    (2, ),
                                                                    dtype = np.float64)
            dset_mat_deriv.attrs['mesh size'] =  min(float(size_parameter), 1.0 / float(size_parameter))
            dset_mat_deriv[:] = self.current_mat_deriv_errors[:]
            dset_laplacian[:] = self.current_laplacian_errors[:]

        if is_regular_mesh:
            size_parameter_name = '_ndiv_'
        else:
            size_parameter_name = '_h_'
        with open('../errors_recorded/mat_deriv_errors'
                  + size_parameter_name
                  + str(self.pp.CFD_DEM.size_parameter)
                  + '_type_'
                  + str(self.pp.CFD_DEM["material_acceleration_calculation_type"].GetInt())
                  + '.txt', 'w') as mat_errors_file:
            for error in self.mat_deriv_errors:
                mat_errors_file.write(str(error) + '\n')
        with open('../errors_recorded/laplacian_errors'
                  + size_parameter_name
                  + str(size_parameter)
                  + '_type_'
                  + str(self.pp.CFD_DEM["laplacian_calculation_type"].GetInt())
                  + '.txt', 'w') as laplacian_errors_file:
            for error in self.laplacian_errors:
                laplacian_errors_file.write(str(error) + '\n')
        sys.stdout.flush()
        sys.stdout.path_to_console_out_file
        os.rename(sys.stdout.path_to_console_out_file,
                  self.post_path +
                  '/' + sys.stdout.console_output_file_name)

        dir_name = self.post_path + '_FINISHED_AT_t=' + str(round(time, int(math.log(10. / time))))

        if os.path.isdir(dir_name):
            import shutil
            shutil.rmtree(dir_name)
        os.rename(self.post_path, dir_name)
