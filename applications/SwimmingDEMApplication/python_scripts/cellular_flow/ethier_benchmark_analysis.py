import KratosMultiphysics as Kratos
from KratosMultiphysics import Vector, Parameters, ModelPartIO
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import sys
import os
import math
import numpy as np
import time as timer
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis as swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis


class EthierBenchmarkAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)
        self.SetCustomBetaParameters()

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        Add = self.project_parameters.AddEmptyValue
        Add("field_identifier").SetString('ethier')
        Add("pressure_grad_recovery_type")
        Add("size_parameter").SetDouble(1)
        Add("store_full_gradient_option").SetBool(True)
        Add("regular_mesh_option").SetBool(True)
        Add("mesh_tag").SetString('Kratos')
        Add("print_VELOCITY_LAPLACIAN_option").SetBool(True)
        Add("print_VECTORIAL_ERROR_option").SetBool(True)

    def SetCustomBetaParameters(self):
        self.field_identifier = self.project_parameters["field_identifier"].GetString()
        self.mesh_tag = self.project_parameters["mesh_tag"].GetString()
        # Creating a code for the used input variables
        self.run_code = '_ndiv_' \
                      + str(self.project_parameters["size_parameter"].GetDouble()) \
                      + '_mat_deriv_type_' \
                      + str(self.project_parameters["material_acceleration_calculation_type"].GetInt()) \
                      + '_lapl_type_' \
                      + str(self.project_parameters["laplacian_calculation_type"].GetInt()) \
                      + '_' + self.field_identifier \
                      + '_' + self.mesh_tag

    def ReadFluidModelParts(self):
        problem_name = self.pp.problem_name.replace('Fluid', '')
        is_regular_mesh = self.project_parameters["regular_mesh_option"].GetBool()
        tag =  self.project_parameters["mesh_tag"].GetString()

        if is_regular_mesh:
            mdpa_name = problem_name + '_ndiv_' + str(int(self.project_parameters["size_parameter"].GetDouble())) + tag + 'Fluid'
        else:
            mdpa_name = problem_name + '_h_' + str(self.project_parameters["size_parameter"].GetDouble()) + 'Fluid'

        model_part_io_fluid = ModelPartIO(mdpa_name)
        model_part_io_fluid.ReadModelPart(self._GetFluidAnalysis().fluid_model_part)

    def AddExtraVariables(self, run_code = ''):
        BaseAnalysis.AddExtraVariables(self, self.run_code)

    def GetParticlesResultsCounter(self):
        return SDP.Counter()

    def GetPrintCounter(self):
        return SDP.Counter(is_dead=True)

    def GetFieldUtility(self):
        a = math.pi / 4
        d = math.pi / 2

        self.flow_field = SDEM.EthierFlowField(a, d)
        space_time_set = SDEM.SpaceTimeSet()
        self.field_utility = SDEM.FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def GetRecoveryCounter(self):
        return SDP.Counter(1, 1, self.project_parameters["coupling"]["coupling_level_type"].GetInt() or self.project_parameters.print_PRESSURE_GRADIENT_option)

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
            node.SetSolutionStepValue(Kratos.VELOCITY, vel)

    def RecoverDerivatives(self):
        t0 = timer.clock()
        self.recovery.Recover()
        t1 = timer.clock()
        self.recovery_time = t1 - t0
        self.CalculateRecoveryErrors(0.0)

    def FluidSolve(self, time = 'None', solve_system = True):
        self.CalculateRecoveryErrors(time)

    def GetFirstStepForFluidComputation(self):
        return 1

    def CalculateRecoveryErrors(self, time):
        L2_norm_mat_deriv = 0.
        L2_norm_mat_deriv_error = 0.
        L2_norm_laplacian = 0.
        L2_norm_laplacian_error = 0.
        max_mat_deriv_error = 0.
        max_laplacian_error = 0.
        total_volume = 0.

        calc_mat_deriv = np.zeros(3)
        calc_laplacian = np.zeros(3)
        mat_deriv= Vector(3)
        laplacian= Vector(3)

        for node in self.fluid_model_part.Nodes:
            nodal_volume = node.GetSolutionStepValue(Kratos.NODAL_AREA)
            total_volume += nodal_volume
            coor = Vector([node.X, node.Y, node.Z])

            self.flow_field.CalculateConvectiveDerivative(0., coor, mat_deriv, 0)
            self.flow_field.CalculateLaplacian(0., coor, laplacian, 0)
            calc_mat_deriv = node.GetSolutionStepValue(Kratos.MATERIAL_ACCELERATION)
            calc_laplacian = node.GetSolutionStepValue(Kratos.VELOCITY_LAPLACIAN)

            module_mat_deriv_squared = sum(x**2 for x in mat_deriv)
            module_laplacian_squared = sum(x**2 for x in laplacian)
            L2_norm_mat_deriv += module_mat_deriv_squared * nodal_volume
            L2_norm_laplacian += module_laplacian_squared * nodal_volume
            diff_mat_deriv = calc_mat_deriv - mat_deriv
            diff_laplacian = calc_laplacian - laplacian
            node.SetSolutionStepValue(Kratos.VECTORIAL_ERROR, Vector(list(diff_mat_deriv)))
            node.SetSolutionStepValue(Kratos.VECTORIAL_ERROR_1, Vector(list(diff_laplacian)))
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
        max_mat_deriv_error **= 0.5
        max_laplacian_error **= 0.5

        if L2_norm_mat_deriv > 0 and L2_norm_laplacian > 0:
            SDP.MultiplyNodalVariableByFactor(self.fluid_model_part, Kratos.VECTORIAL_ERROR, 1.0 / L2_norm_mat_deriv)
            SDP.MultiplyNodalVariableByFactor(self.fluid_model_part, Kratos.VECTORIAL_ERROR_1, 1.0 / L2_norm_laplacian)
            self.current_mat_deriv_errors[0] = L2_norm_mat_deriv_error / L2_norm_mat_deriv
            self.current_mat_deriv_errors[1] = max_mat_deriv_error / L2_norm_mat_deriv
            self.current_laplacian_errors[0] = L2_norm_laplacian_error / L2_norm_laplacian
            self.current_laplacian_errors[1] = max_laplacian_error / L2_norm_laplacian
            self.mat_deriv_errors.append(self.current_mat_deriv_errors)
            self.laplacian_errors.append(self.current_laplacian_errors)

            text_width = 40
            Kratos.Logger.PrintInfo("SwimmingDEM", '\n' + '-.' * text_width)
            Kratos.Logger.PrintInfo("SwimmingDEM", 'L2 error for the material derivative'.ljust(text_width), self.current_mat_deriv_errors[0])
            Kratos.Logger.PrintInfo("SwimmingDEM", 'max error for the material derivative'.ljust(text_width), self.current_mat_deriv_errors[1])
            Kratos.Logger.PrintInfo("SwimmingDEM", 'L2 error for the laplacian'.ljust(text_width), self.current_laplacian_errors[0])
            Kratos.Logger.PrintInfo("SwimmingDEM", 'max error for the laplacian'.ljust(text_width), self.current_laplacian_errors[1])
            Kratos.Logger.PrintInfo('-.' * text_width + '\n')

    def PerformFinalOperations(self, time = None):

        if not os.path.exists('../errors_recorded'):
            os.makedirs('../errors_recorded')

        import h5py

        file_name = self.main_path + '/errors_recorded/recovery_errors.hdf5'
        # with h5py.File(self.file_name, 'r+') as f:
        #     f.create_dataset('material_derivative', shape = self.shape, dtype = np.float32)
        size_parameter = self.project_parameters["size_parameter"].GetInt()
        is_regular_mesh = self.project_parameters["regular_mesh_option"].GetBool()

        with h5py.File(file_name) as f:
            field_identifier = self.project_parameters["field_identifier"].GetString()
            field_group = f.require_group(field_identifier)
            mat_deriv_grp = field_group.require_group('material derivative')
            mat_deriv_mthd_group = mat_deriv_grp.require_group('method = ' + str(self.project_parameters["material_acceleration_calculation_type"].GetInt()))
            laplacian_grp = field_group.require_group('laplacian')
            laplacian_mthd_group = laplacian_grp.require_group('method = ' + str(self.project_parameters["laplacian_calculation_type"].GetInt()))

            if is_regular_mesh:
                mesh_tag = self.project_parameters["mesh_tag"].GetString()
                mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('regular mesh (' + mesh_tag + ')')
                mesh_grp_laplacian = laplacian_mthd_group.require_group('regular mesh (' + mesh_tag + ')')
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
            dset_mat_deriv.attrs['recovery time'] = self.recovery_time
            dset_laplacian.attrs['recovery time'] = self.recovery_time
            dset_mat_deriv.attrs['mesh size'] =  min(float(size_parameter), 1.0 / float(size_parameter))
            dset_laplacian.attrs['mesh size'] =  min(float(size_parameter), 1.0 / float(size_parameter))
            dset_mat_deriv[:] = self.current_mat_deriv_errors[:]
            dset_laplacian[:] = self.current_laplacian_errors[:]

        sys.stdout.flush()
        os.rename(sys.stdout.path_to_console_out_file,
                  self.post_path +
                  '/' + sys.stdout.console_output_file_name)

        if time <= 0:
            time_string = '0.0'
        else:
            time_string = str(round(time, int(math.log(10. / time))))
        dir_name = self.post_path + '_FINISHED_AT_t=' + time_string

        if os.path.isdir(dir_name):
            import shutil
            shutil.rmtree(dir_name)
        os.rename(self.post_path, dir_name)
