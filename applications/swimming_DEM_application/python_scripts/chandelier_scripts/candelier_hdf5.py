import h5py
import numpy as np
import math
from KratosMultiphysics import *
import chandelier_parameters as ch_pp
import chandelier as ch

class ResultsCandelier:
    def __init__(self, pp, path):
        self.sim = ch.AnalyticSimulator(ch_pp)
        self.sim.CalculateNonDimensionalVars()
        self.path = path + '/candelier_results.h5py'
        self.dt = pp.CFD_DEM["MaxTimeStep"].GetDouble()
        self.N_q = pp.CFD_DEM["time_steps_per_quadrature_step"].GetInt()
        self.quadrature_order = pp.CFD_DEM["quadrature_order"].GetInt()
        self.reading_index = 0
        self.times = []
        self.errors = []
        ch_pp.include_history_force = bool(pp.CFD_DEM["basset_force_type"].GetInt())

        if pp.CFD_DEM["basset_force_type"].GetInt() == 2:
            self.method = 'Daitche'
        else:
            self.method = 'Hinsberg'
            self.m = pp.CFD_DEM["number_of_exponentials"].GetInt()
            self.t_w = pp.CFD_DEM["time_window"].GetDouble()

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
