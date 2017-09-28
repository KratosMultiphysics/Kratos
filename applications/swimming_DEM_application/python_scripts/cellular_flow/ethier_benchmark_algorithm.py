from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_algorithm
import swimming_DEM_procedures as SDP
import math
import numpy as np
BaseAlgorithm = swimming_DEM_algorithm.Algorithm

def num_type(value):
  try:
    int(str(value))
    return 'int'
  except:
    try:
        float(str(value))
        return 'float'
    except:
        return None

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)

    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)
        self.pp.CFD_DEM.AddEmptyValue("pressure_grad_recovery_type")
        self.pp.CFD_DEM.AddEmptyValue("size_parameter").SetInt(1)
        self.pp.CFD_DEM.AddEmptyValue("store_full_gradient_option").SetBool(True)

    def SetCustomBetaParameters(self, custom_parameters):
        BaseAlgorithm.SetCustomBetaParameters(self, custom_parameters)
        self.pp.CFD_DEM.size_parameter = self.pp.CFD_DEM["size_parameter"].GetInt()
        # Creating a code for the used input variables
        self.run_code = '_ndiv_' + str(self.pp.CFD_DEM["size_parameter"].GetDouble()) \
                                 + '_mat_deriv_type_' \
                                 + str(self.pp.CFD_DEM["material_acceleration_calculation_type"].GetInt()) \
                                 + '_lapl_type_' \
                                 + str(self.pp.CFD_DEM["laplacian_calculation_type"].GetInt())

    def ReadFluidModelParts(self):
        os.chdir(self.main_path)
        if num_type(self.pp.CFD_DEM.size_parameter) == 'int':
            model_part_io_fluid = ModelPartIO(self.pp.problem_name.replace('ethier', 'ethier_ndiv_' + str(self.pp.CFD_DEM.size_parameter)))
        elif num_type(self.pp.CFD_DEM.size_parameter) == 'float':
            model_part_io_fluid = ModelPartIO(self.pp.problem_name.replace('ethier', 'ethier_h_' + str(self.pp.CFD_DEM.size_parameter)))
        model_part_io_fluid.ReadModelPart(self.fluid_algorithm.fluid_model_part)

    def AddExtraVariables(self, run_code = ''):
        BaseAlgorithm.AddExtraVariables(self, self.run_code)

    def GetFieldUtility(self):
        a = math.pi / 4
        d = math.pi / 2*0

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
        fluid_model_part = self.all_model_parts.Get('FluidPart')
        for node in fluid_model_part.Nodes:
            vel = Vector(3)
            coor = Vector(3)
            coor[0]=node.X
            coor[1]=node.Y
            coor[2]=node.Z
            self.flow_field.Evaluate(0.0, coor, vel, 0)
            node.SetSolutionStepValue(VELOCITY_X, vel[0])
            node.SetSolutionStepValue(VELOCITY_Y, vel[1])
            node.SetSolutionStepValue(VELOCITY_Z, vel[2])

    def FluidSolve(self, time = 'None'):
        fluid_model_part = self.all_model_parts.Get('FluidPart')

        for node in fluid_model_part.Nodes:
            vel= Vector(3)
            coor= Vector(3)
            coor[0]=node.X
            coor[1]=node.Y
            coor[2]=node.Z
            self.flow_field.Evaluate(time,coor,vel,0)
            node.SetSolutionStepValue(VELOCITY_X, vel[0])
            node.SetSolutionStepValue(VELOCITY_Y, vel[1])
            node.SetSolutionStepValue(VELOCITY_Z, vel[2])
            #node.SetSolutionStepValue(VELOCITY_X, 7*node.X**2 + 6 * node.Y + 19)
            #node.SetSolutionStepValue(VELOCITY_Y,  9 *node.X**2 - 8 * node.Y**2 +30*node.X)
            #node.SetSolutionStepValue(VELOCITY_Z, 0.0)
        self.CalculateRecoveryErrors(time)

    def CalculateRecoveryErrors(self, time):
        error_mat_deriv = 0.
        max_error_mat_deriv = - float('inf')
        error_laplacian = 0.
        max_error_laplacian = - float('inf')

        total_volume = 0.
        mat_deriv_average = Vector(3)
        laplacian_average = Vector(3)
        for k in range(3):
            mat_deriv_average[k] = 0.0
            laplacian_average[k] = 0.0
        norm_mat_deriv_average = 0.
        norm_laplacian_average = 0.

        module_mat_deriv = 0.
        module_laplacian = 0.
        fluid_model_part = self.all_model_parts.Get('FluidPart')

        for i_node, node in enumerate(fluid_model_part.Nodes):
            calc_mat_deriv = [0.] * 3
            calc_laplacian = [0.] * 3
            mat_deriv= Vector(3)
            laplacian= Vector(3)
            coor= Vector(3)
            coor[0]=node.X
            coor[1]=node.Y
            coor[2]=node.Z
            self.flow_field.CalculateMaterialAcceleration(time, coor, mat_deriv, 0)
            self.flow_field.CalculateLaplacian(time, coor, laplacian, 0)
            calc_mat_deriv[0] = node.GetSolutionStepValue(MATERIAL_ACCELERATION_X)
            calc_mat_deriv[1] = node.GetSolutionStepValue(MATERIAL_ACCELERATION_Y)
            calc_mat_deriv[2] = node.GetSolutionStepValue(MATERIAL_ACCELERATION_Z)
            calc_laplacian[0] = node.GetSolutionStepValue(VELOCITY_LAPLACIAN_X)
            calc_laplacian[1] = node.GetSolutionStepValue(VELOCITY_LAPLACIAN_Y)
            calc_laplacian[2] = node.GetSolutionStepValue(VELOCITY_LAPLACIAN_Z)
            module_mat_deriv += math.sqrt(calc_mat_deriv[0] ** 2 + calc_mat_deriv[1] ** 2 + calc_mat_deriv[2] ** 2)
            module_laplacian += math.sqrt(calc_laplacian[0] ** 2 + calc_laplacian[1] ** 2 + calc_laplacian[2] ** 2)
            #module_mat_deriv = max(math.sqrt(mat_deriv[0] ** 2 + mat_deriv[1] ** 2 + mat_deriv[2] ** 2), 1e-8)
            #module_laplacian = max(math.sqrt(laplacian[0] ** 2 + laplacian[1] ** 2 + laplacian[2] ** 2), 1e-8)
            #laplacian[0] = 14
            #laplacian[1] = 2
            #laplacian[2] = 0
            #nodal_volume = node.GetSolutionStepValue(NODAL_AREA)
            #total_volume += nodal_volume
            current_error = SDP.NormOfDifference(calc_mat_deriv, mat_deriv)
            error_mat_deriv += current_error
            max_error_mat_deriv = max(max_error_mat_deriv, current_error)
            current_error = SDP.NormOfDifference(calc_laplacian, laplacian)
            error_laplacian += current_error
            max_error_laplacian = max(max_error_laplacian, current_error)
            diff_mat_deriv = [calc_mat_deriv[i] - mat_deriv[i] for i in range(len(calc_mat_deriv))]
            diff_laplacian = [calc_laplacian[i] - laplacian[i] for i in range(len(calc_laplacian))]
            #mat_deriv_averager.SDP.Norm(diff_mat_deriv)
            #laplacian_averager.SDP.Norm(diff_laplacian)
            #for k in range(3):
                #mat_deriv_average[k] += mat_deriv[k]
                #laplacian_average[k] += laplacian[k]
            norm_mat_deriv_average += SDP.Norm(mat_deriv)
            norm_laplacian_average += SDP.Norm(laplacian)

            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_RATE_X, calc_mat_deriv[0] - mat_deriv[0])
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_RATE_Y, calc_mat_deriv[1] - mat_deriv[1])
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_RATE_Z, calc_mat_deriv[2] - mat_deriv[2])
            #node.SetSolutionStepValue(MATERIAL_ACCELERATION_X, mat_deriv[0])
            #node.SetSolutionStepValue(MATERIAL_ACCELERATION_Y, mat_deriv[1])
            #node.SetSolutionStepValue(MATERIAL_ACCELERATION_Z, mat_deriv[2])


            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, calc_laplacian[0] - laplacian[0])
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, calc_laplacian[1] - laplacian[1])
            node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, calc_laplacian[2] - laplacian[2])
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, calc_laplacian_0)
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, calc_laplacian_1)
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, calc_laplacian_2)
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_X, laplacian[0])
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Y, laplacian[1])
            #node.SetSolutionStepValue(VELOCITY_LAPLACIAN_Z, laplacian[2])

        module_mat_deriv /= len(fluid_model_part.Nodes)
        module_laplacian /= len(fluid_model_part.Nodes)
        SDP.MultiplyNodalVariableByFactor(fluid_model_part, VELOCITY_LAPLACIAN_RATE, 1.0 / module_mat_deriv)
        SDP.MultiplyNodalVariableByFactor(fluid_model_part, VELOCITY_LAPLACIAN, 1.0 / module_laplacian)

        if norm_mat_deriv_average > 0. and norm_laplacian_average > 0:
            self.current_mat_deriv_errors[0] = error_mat_deriv / norm_mat_deriv_average
            self.current_mat_deriv_errors[1] = max_error_mat_deriv / norm_mat_deriv_average * len(fluid_model_part.Nodes)
            self.current_laplacian_errors[0] = error_laplacian / norm_laplacian_average
            self.current_laplacian_errors[1] = max_error_laplacian / norm_laplacian_average * len(fluid_model_part.Nodes)
            self.mat_deriv_errors.append(self.current_mat_deriv_errors)
            self.laplacian_errors.append(self.current_laplacian_errors)
            #print('mat_deriv: min, max, avg, ', mat_deriv_averager.GetCurrentData())
            #print('laplacian: min, max, avg, ', laplacian_averager.GetCurrentData())
            print('rel_error_mat_deriv', error_mat_deriv / norm_mat_deriv_average)
            print('rel_error_laplacian', error_laplacian / norm_laplacian_average)

    def PerformFinalOperations(self, time = None):

        if not os.path.exists('../errors_recorded'):
            os.makedirs('../errors_recorded')

        import h5py

        file_name = self.main_path + '/errors_recorded/recovery_errors.hdf5'
        # with h5py.File(self.file_name, 'r+') as f:
        #     f.create_dataset('material_derivative', shape = self.shape, dtype = np.float32)
        size_parameter = self.pp.CFD_DEM.size_parameter
        with h5py.File(file_name) as f:
            mat_deriv_grp = f.require_group('material derivative')
            mat_deriv_mthd_group = mat_deriv_grp.require_group('method = ' + str(self.pp.CFD_DEM["material_acceleration_calculation_type"].GetInt()))
            laplacian_grp = f.require_group('laplacian')
            laplacian_mthd_group = laplacian_grp.require_group('method = ' + str(self.pp.CFD_DEM["laplacian_calculation_type"].GetInt()))

            if num_type(size_parameter) == 'int':
                mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('regular mesh')
                mesh_grp_laplacian = laplacian_mthd_group.require_group('regular mesh')
                dset_mat_deriv = mesh_grp_mat_deriv.require_dataset('n_div = ' + str(size_parameter), (2,), dtype = np.float64)
                dset_laplacian = mesh_grp_laplacian.require_dataset('n_div = ' + str(size_parameter), (2,), dtype = np.float64)
            elif num_type(size_parameter) == 'float':
                mesh_grp_mat_deriv = mat_deriv_mthd_group.require_group('irregular mesh')
                mesh_grp_laplacian = laplacian_mthd_group.require_group('irregular mesh')
                dset_mat_deriv = mesh_grp_mat_deriv.require_dataset('h = ' + str(size_parameter), (2,), dtype = np.float64)
                dset_laplacian = mesh_grp_laplacian.require_dataset('h = ' + str(size_parameter), (2,), dtype = np.float64)
            dset_mat_deriv.attrs['mesh size'] =  min(float(size_parameter), 1.0 / float(size_parameter))
            dset_mat_deriv[0] = self.current_mat_deriv_errors[0]
            dset_mat_deriv[1] = self.current_mat_deriv_errors[1]
            dset_laplacian[0] = self.current_laplacian_errors[0]
            dset_laplacian[1] = self.current_laplacian_errors[1]

        if num_type(size_parameter) == 'int':
            size_parameter_name = '_ndiv_'
        elif num_type(size_parameter) == 'float':
            size_parameter_name = '_h_'
        with open('../errors_recorded/mat_deriv_errors' + size_parameter_name + str(self.pp.CFD_DEM.size_parameter) + '_type_' + str(self.pp.CFD_DEM["material_acceleration_calculation_type"].GetInt()) + '.txt', 'w') as mat_errors_file:
            for error in self.mat_deriv_errors:
                mat_errors_file.write(str(error) + '\n')
        with open('../errors_recorded/laplacian_errors' + size_parameter_name + str(size_parameter) + '_type_' + str(self.pp.CFD_DEM["laplacian_calculation_type"].GetInt()) + '.txt', 'w') as laplacian_errors_file:
            for error in self.laplacian_errors:
                laplacian_errors_file.write(str(error) + '\n')
        sys.stdout.flush()
        os.chdir(self.main_path)
        sys.stdout.path_to_console_out_file
        os.rename(sys.stdout.path_to_console_out_file, self.post_path + '/' + sys.stdout.console_output_file_name)
        empty = False
        add_to_name = ''
        count = 0
        dir_name = self.post_path + '_FINISHED_AT_t=' + str(round(time, int(math.log(10. / time))))
        #dir_name_it = dir_name
        #while os.path.isdir(dir_name_it):
        #dir_name_it = dir_name + '_version_' + str(count + 1)
        #count += 1
        if os.path.isdir(dir_name):
            import shutil
            shutil.rmtree(dir_name)
            print(dir_name)
            print(self.post_path)
        os.rename(self.post_path, dir_name)
