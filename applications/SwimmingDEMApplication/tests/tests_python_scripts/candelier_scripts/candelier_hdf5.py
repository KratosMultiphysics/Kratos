import h5py
import numpy as np
import math
import candelier
import candelier_parameters as ch_pp
import parameters_tools as PT

class ResultsCandelier:
    def __init__(self, project_parameters, path):
        self.sim = candelier.AnalyticSimulator(ch_pp)
        self.sim.CalculateNonDimensionalVars()
        self.path = path + '/candelier_results.h5py'
        self.reading_index = 0
        self.times = []
        self.errors = []
        self.do_include_history_force = (PT.RecursiveFindParametersWithCondition(
                                         project_parameters["properties"], 'history_force_parameters',
                                         condition=lambda value: value['name'].GetString() != 'default'))

        ch_pp.include_history_force = self.do_include_history_force
        self.dt = project_parameters["time_stepping"]["time_step"].GetDouble()

        if self.do_include_history_force: #TODO: extend to multiple properties
            for prop in project_parameters["properties"].values():
                self.history_force_parameters =  prop["hydrodynamic_law_parameters"]["history_force_parameters"]
                break

            self.N_q = self.history_force_parameters["time_steps_per_quadrature_step"].GetInt()
            self.quadrature_order = self.history_force_parameters["quadrature_order"].GetInt()

        if not self.history_force_parameters["mae_parameters"]['do_use_mae'].GetBool():
            self.method = 'Daitche'
        else:
            self.method = 'Hinsberg'
            self.m = self.history_force_parameters["mae_parameters"]["m"].GetInt()
            self.t_w = self.history_force_parameters["mae_parameters"]["window_time_interval"].GetDouble()

        self.result_code = self.method + '_dt=' + str(self.dt) + '_Nq=' + str(self.N_q) + '_quadrature_order=' + str(self.quadrature_order)
        if self.method == 'Hinsberg':
            self.result_code += '_tw=' + str(self.t_w) + '_m=' + str(self.m)

        with h5py.File(self.path) as f:
            result = f.require_group(self.result_code)
            result.attrs['method'] = self.method
            result.attrs['dt'] = self.dt
            result.attrs['N_q'] = self.N_q
            result.attrs['quadrature_order'] = self.quadrature_order

            if self.method == 'Hinsberg':
                result.attrs['t_w'] = self.t_w
                result.attrs['m'] = self.m

    def WriteToHDF5(self):
        with h5py.File(self.path, 'r+') as f:
            times = np.array(self.times)
            errors = np.array(self.errors)
            shape = (len(times),)
            f[self.result_code].require_dataset(name = 'time', data = times, shape = shape, dtype = 'float64')
            f[self.result_code].require_dataset(name = 'E(t)', data = errors, shape = shape, dtype = 'float64')

    def CalculateError(self, time, coor_calculated):
        coor_theory = np.zeros(3)
        vel_theory = np.zeros(3)
        self.sim.CalculatePosition(coor_theory, time * ch_pp.omega, vel_theory)
        coor_theory *= ch_pp.R
        try:
            r_inv = 1.0 / sum([x ** 2 for x in coor_theory[:-1]])
        except:
            r_inv = 0

        coor_calculated[2] = coor_theory[2]
        error = math.sqrt(r_inv * sum([(x - y) ** 2 for (x, y) in zip(coor_theory[:-1], coor_calculated[:-1])]))
        return error

    # def MakeReading(self, time, coor_calculated):
    #     with h5py.File(self.path, 'r+') as f:
    #         entry = str(self.reading_index).zfill(8)
    #         f[self.result_code].require_group(entry)
    #         f[self.result_code + '/' + entry].attrs['time'] = time
    #         f[self.result_code + '/' + entry].attrs['E(time)'] = self.CalculateError(time, coor_calculated)
    #         self.reading_index += 1

    def MakeReading(self, time, coor_calculated):
        self.times.append(time)
        self.errors.append(self.CalculateError(time, coor_calculated))
