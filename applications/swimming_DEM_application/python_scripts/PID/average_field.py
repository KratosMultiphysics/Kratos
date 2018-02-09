import math
import numpy as np
import h5py

def GetLastStationaryFieldId(hdf5_file):
    max_Id = - 1
    valid_entries = (entry for entry in hdf5_file if 'stationary_field' in entry)
    for entry in valid_entries:
        max_Id = max(max_Id, hdf5_file['/' + entry].attrs['Id'])
        print('max_Id', max_Id)

    return max_Id

def CreateGroup(file_or_group, name, overwrite_previous=True):
    if overwrite_previous:
        if name in file_or_group:
            file_or_group['/'].__delitem__(name)
        Id = 0
    else:
        Id = GetLastStationaryFieldId(file_or_group) + 1
        if Id:
            name += '_' + str(Id)
        print(name)

    group = file_or_group.create_group(name)
    group.attrs['Id'] = Id

    return group

class Rotator:
    def __init__(self,
                 rotation_axis_initial_point,
                 rotation_axis_final_point,
                 angular_velocity_module):
        self.a_init = np.array(rotation_axis_initial_point)
        self.a_final = np.array(rotation_axis_final_point)
        self.omega = angular_velocity_module
        self.axis = np.array(self.Normalize(self.a_final - self.a_init))
        self.CalculateRodriguesMatrices(self.axis)

    def CalculateRodriguesMatrices(self, axis):
        self.I = np.identity(3)
        self.UU = np.array([a * axis for a in axis])
        self.Ux = np.array([[0, - axis[2], axis[1]],
                            [axis[2], 0., -axis[0]],
                            [-axis[1], axis[0], 0.]])

    def GetRotationMatrices(self, time):
        sin = math.sin(self.omega * time)
        cos = math.cos(self.omega * time)

        # Rotation matrix
        R = cos * self.I + sin * self.Ux + (1.0 - cos) * self.UU

        # Rotation matrix' (derivative ofname R with respect to time)
        Rp = - self.omega * sin * self.I + self.omega * cos * self.Ux + self.omega * sin * self.UU

        return R, Rp

    def UndoRotationOfVectors(self, time, list_of_vectors):
        R, Rp = self.GetRotationMatrices(time)
        R_inv = np.linalg.inv(R)
        list_of_vectors = np.dot(R_inv, list_of_vectors)

    def Normalize(self, v):
        mod_2 = sum([x ** 2 for x in v])

        if mod_2 == 0:
            return v
        else:
            mod_inv = 1.0 / math.sqrt(mod_2)
            return mod_inv * v

class Averager:
    def __init__(self,
                 rotation_axis_initial_point,
                 rotation_axis_final_point,
                 angular_velocity_module,
                 dataset_name,
                 original_file_name,
                 original_file_path,
                 initial_time = 0.0,
                 final_time = float('inf'),
                 steps_per_average_step = 1,
                 calculate_standard_deviations = True,
                 normalize_standard_deviation = True,
                 overwrite_previous=True):
        self.initial_time = initial_time
        self.final_time = final_time
        self.steps_per_average_step = steps_per_average_step
        self.calculate_standard_deviations = calculate_standard_deviations
        self.dataset_name = dataset_name
        self.original_file_name = original_file_name
        self.original_file_path = original_file_path + '/' + original_file_name
        self.normalize_standard_deviation = normalize_standard_deviation
        self.overwrite_previous = overwrite_previous
        self.rotator = Rotator(rotation_axis_initial_point, rotation_axis_final_point, angular_velocity_module)

    def IsRelevant(self, i, time_str):
        if i % self.steps_per_average_step != 0:
            return False
        else:
            time = float(time_str)
            return time >= self.initial_time and time < self.final_time

    def GetFilePath(self):
        return self.averaged_field_file_name

    def PerformAverage(self, reference_time=0.0):
        self.reference_time = reference_time

        with h5py.File(self.original_file_path, 'r') as f:
            self.times_str = [t for t in f if 'time' in f['/' + t].attrs]
            self.n_nodes = f['/nodes'][:, 0].size

            p = np.array(self.n_nodes)
            vel = np.zeros((3, self.n_nodes))
            field = np.zeros((4, self.n_nodes))

            relevant_time_strings = [t for i, t in enumerate(self.times_str) if self.IsRelevant(i, t)]

            for i_sample, time_str in enumerate(relevant_time_strings, start=1):
                p = f['/' + time_str + '/p']
                vel[0] = f['/' + time_str + '/vx']
                vel[1] = f['/' + time_str + '/vy']
                vel[2] = f['/' + time_str + '/vz']
                self.rotator.UndoRotationOfVectors(float(time_str) - self.reference_time, vel)
                field[:3] += vel
                field[3] += p
                print('averaging step ', i_sample, '...')

            field /= i_sample

            if self.calculate_standard_deviations:
                std_dev = np.zeros((2, self.n_nodes))

                for i_sample, time_str in enumerate(relevant_time_strings, start=1):
                    p = f['/' + time_str + '/p']
                    vel[0] = f['/' + time_str + '/vx']
                    vel[1] = f['/' + time_str + '/vy']
                    vel[2] = f['/' + time_str + '/vz']
                    self.rotator.UndoRotationOfVectors(float(time_str) - self.reference_time, vel)
                    for i in range(3):
                        std_dev[0] += (field[i] - vel[i]) ** 2
                    std_dev[1] += (field[3] - p) ** 2
                    print('calculating standard deviations step ', i_sample, '...')


            std_dev **= 0.5
            std_dev /= i_sample

            if self.normalize_standard_deviation:
                norms = np.sqrt(field[0] ** 2 + field[1] ** 2 + field[2] ** 2)

                with np.errstate(divide='ignore', invalid='ignore'):
                    std_dev[0] = np.where(norms > 0.0, std_dev[0] / norms, 0.0)
                    std_dev[1] = np.where(abs(field[3]) > 0., std_dev[1] / abs(field[3]), 0.0)

            v_stdv_modulus = sum(std_dev[0]) / self.n_nodes
            p_stdv_modulus = sum(std_dev[1]) / self.n_nodes

            self.averaged_field_file_name = self.original_file_path.replace('.hdf5', '') + '_averaged.hdf5'

            with h5py.File(self.averaged_field_file_name) as f:
                f.attrs['initial_time'] = self.initial_time
                f.attrs['final_time'] = float(time_str)
                f.attrs['number_of_samples'] = i_sample
                stat_group = CreateGroup(file_or_group=f,
                                         name=self.dataset_name,
                                         overwrite_previous=self.overwrite_previous)
                stat_group.create_dataset('vx', data = field[0])
                stat_group.create_dataset('vy', data = field[1])
                stat_group.create_dataset('vz', data = field[2])
                stat_group.create_dataset('p', data = field[3])

                if self.calculate_standard_deviations:
                    stat_group.attrs['reference_time'] = self.reference_time
                    stat_group.attrs['velocity_standard_deviation'] = v_stdv_modulus
                    stat_group.attrs['pressure_standard_deviation'] = p_stdv_modulus
                    stat_group.create_dataset('v_mod_stdv', data = std_dev[0])
                    stat_group.create_dataset('p_stdv', data = std_dev[1])

            print('velocity modulus mean standard deviation:', v_stdv_modulus)
            print('pressure mean standard deviation:', p_stdv_modulus)


if __name__ == '__main__':
    import os
    averager = Averager(rotation_axis_initial_point = [0., 0., 0.],
                        rotation_axis_final_point = [0., 0., 1.],
                        angular_velocity_module = - 2 * math.pi,
                        dataset_name = 'stationary_field',
                        original_file_name = 'sample.hdf5',
                        original_file_path = os.getcwd(),
                        initial_time = 0.04,
                        steps_per_average_step = 50,
                        calculate_standard_deviations = True,
                        normalize_standard_deviation = True,
                        overwrite_previous=False)
    averager.PerformAverage()
