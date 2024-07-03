import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import h5py
from fractions import Fraction
from matplotlib import colors
from decimal import Decimal


plt.rcParams.update({'font.size': 13.5})
pd.set_option("display.max_columns", None)

class Curve:
    def __init__(self, filename, norm):
        self.filename = filename
        self.norm = norm
        self.data_base = self.CreateDataBase(self.filename)

    def CreateDataBase(self, filename):
        hdata = []
        v_df = []
        p_df = []
        projection_type = []
        model_type = []
        subscale_type = []
        delta_time = []
        lowest_alpha = []
        convergence = []
        n_case = []
        relaxation = []
        damkohler_number = []
        reynolds_number = []
        results = pd.DataFrame()
        f = h5py.File(filename, 'r')
        keys = sorted([int(key) for key in f.keys()])
        for key in keys:
            groups = f[str(key)]
            v_error = groups[('V_' + self.norm + '_ERROR')]
            p_error = groups[('P_' + self.norm + '_ERROR')]
            n_iter = groups[('N_ITERATIONS')]
            #element_size = 2*float(groups.attrs['element_size'])
            if key == 0:
                element_size = np.sqrt(2)
            if key == 1:
                element_size = np.sqrt(2)/2
            if key == 2:
                element_size = np.sqrt(2)/4
            proj_data = groups.attrs['projection_type']
            model_data = groups.attrs['model_type']
            sub_data = groups.attrs['subscale_type']
            dt_data = groups.attrs['delta_time']
            relaxation_data = groups.attrs['relaxation_alpha']
            max_iter_data = float(groups.attrs['max_iteration'])
            lowest_alpha_data = float(groups.attrs['lowest_alpha'])
            damkohler_number_data = float(groups.attrs['damkohler_number'])
            reynolds_number_data = float(groups.attrs['reynolds_number'])
            n_case.append(key)
            hdata.append(round(element_size,4))
            projection_type.append(proj_data)
            model_type.append(model_data)
            subscale_type.append(sub_data)
            delta_time.append(round(float(dt_data),8))
            lowest_alpha.append(round(lowest_alpha_data,15))
            damkohler_number.append(round(damkohler_number_data,15))
            reynolds_number.append(round(reynolds_number_data,8))
            relaxation.append(relaxation_data)
            array_v_error = np.array(v_error)
            array_p_error = np.array(p_error)
            n_iterations = np.array(n_iter)
            vdata = array_v_error[0:]
            pdata = array_p_error[0:]
            vdata_avg = vdata[-1]
            pdata_avg = pdata[-1]
            if np.any(n_iterations == max_iter_data):
                if vdata[-1] >= 1e10 or np.isnan(vdata[-1]):
                    convergence.append('not_converged')
                else:
                    convergence.append('may_converge')
            else:
                convergence.append('converges')
            v_df.append(vdata_avg)
            p_df.append(pdata_avg)
        v_array = (np.array(v_df)).T
        p_array = (np.array(p_df)).T
        results['sp_data_v'] = v_array
        results['sp_data_p'] = p_array
        results['element_size'] = (np.array(hdata)).T
        results['projection_type'] = (np.array(projection_type)).T
        results['model_type'] = (np.array(model_type)).T
        results['subscale_type'] = (np.array(subscale_type)).T
        results['delta_time'] = (np.array(delta_time)).T
        results['lowest_alpha'] = (np.array(lowest_alpha)).T
        results['convergence'] = (np.array(convergence)).T
        results['Da'] = (np.array(damkohler_number)).T
        results['case_number'] = (np.array(n_case)).T
        results['relaxation'] = (np.array(relaxation)).T
        results['Re'] = (np.array(reynolds_number)).T

        return results.sort_values(['element_size'], ascending=False).reset_index(drop = True)

class Data:

    def __init__(self, data_base, norm, dof, set_index):
        self.initial_data_base = data_base
        self.norm_type = norm
        self.dof = dof
        db = self.initial_data_base.loc[(self.initial_data_base['subscale_type'] == 'quasi-static')].reset_index(drop = True)
        self.data_list = self.ClasifyDataBase(db)
        self.result_db = pd.DataFrame()
        for datb in self.data_list:
            if not set_index:
                self.result_db['case'] = self.ExtractProblem(datb)
                set_index = True
            column_name_min_error = self.dof + ' ' + self.norm_type + ' ' + datb['projection_type'].iloc[0] + ' ' + 'min'
            column_name_slope = self.dof + ' ' + self.norm_type + ' ' + datb['projection_type'].iloc[0] + ' ' + 'slope'
            self.result_db[column_name_slope] = self.ExtractConvergenceSlope(datb)
            self.result_db[column_name_min_error] = self.ExtractMinimumError(datb)

    def ClasifyDataBase(self, db):
        db_1 = db.loc[(db['projection_type'] == 'ASGS') & (db['model_type'] == 'Drew model') & (db['convergence'] == 'converges')].sort_values(['element_size'], ascending=False).reset_index(drop = True)
        db_2 = db.loc[(db['projection_type'] == 'OSS') & (db['model_type'] == 'Drew model') & (db['convergence'] == 'converges')].sort_values(['element_size'], ascending=False).reset_index(drop = True)
        if db_1.empty:
            return [db_2]
        elif db_2.empty:
            return [db_1]
        else:
            return [db_1,db_2]

    def ExtractMinimumError(self,db):
        lowest_alpha_set = sorted(self.initial_data_base['lowest_alpha'].unique().tolist(),reverse=True)
        damkohler_set = sorted(self.initial_data_base['Da'].unique().tolist(),reverse=False)
        reynolds_set = sorted(self.initial_data_base['Re'].unique().tolist(),reverse=False)
        aux = []

        #print('Minimum ' + self.norm_type + ' ' + self.dof + ' error values')
        for d in damkohler_set:
            for r in reynolds_set:
                for lowest_alpha in lowest_alpha_set:
                    if not db.empty:
                        sdb = db.loc[(db['Da'] == d) & (db['lowest_alpha'] == lowest_alpha) & (db['Re'] == r)].sort_values(['element_size'], ascending=False)
                        #print(r,d,lowest_alpha)
                        if not sdb.empty:
                            if self.dof == 'v':
                                aux.append(min(sdb['sp_data_v']))
                            else:
                                aux.append(min(sdb['sp_data_p']))
        return aux


    def ExtractConvergenceSlope(self,db):
        lowest_alpha_set = sorted(self.initial_data_base['lowest_alpha'].unique().tolist(),reverse=True)
        damkohler_set = sorted(self.initial_data_base['Da'].unique().tolist(),reverse=False)
        reynolds_set = sorted(self.initial_data_base['Re'].unique().tolist(),reverse=False)
        aux = []
        for d in damkohler_set:
            for r in reynolds_set:
                for lowest_alpha in lowest_alpha_set:
                    if not db.empty:
                        sdb = db.loc[(db['Da'] == d) & (db['lowest_alpha'] == lowest_alpha) & (db['Re'] == r)].sort_values(['element_size'], ascending=False)
                        #print(r,d,lowest_alpha)
                        if not sdb.empty:
                            if self.dof == 'v':
                                aux.append((np.log(sdb['sp_data_v'].iloc[-2]) - np.log(sdb['sp_data_v'].iloc[-1]))/(np.log(sdb['element_size'].iloc[-2]) - np.log(sdb['element_size'].iloc[-1])))
                            else:
                                aux.append((np.log(sdb['sp_data_p'].iloc[-2]) - np.log(sdb['sp_data_p'].iloc[-1]))/(np.log(sdb['element_size'].iloc[-2]) - np.log(sdb['element_size'].iloc[-1])))
        return aux

    def ExtractProblem(self,db):
        lowest_alpha_set = sorted(self.initial_data_base['lowest_alpha'].unique().tolist(),reverse=True)
        damkohler_set = sorted(self.initial_data_base['Da'].unique().tolist(),reverse=False)
        reynolds_set = sorted(self.initial_data_base['Re'].unique().tolist(),reverse=False)
        aux = []
        for d in damkohler_set:
            for r in reynolds_set:
                for lowest_alpha in lowest_alpha_set:
                    if not db.empty:
                        sdb = db.loc[(db['Da'] == d) & (db['lowest_alpha'] == lowest_alpha) & (db['Re'] == r)].sort_values(['element_size'], ascending=False)
                        case = '(' +str(r) + ', ' + str(d) + ', ' + str(lowest_alpha) + ')'
                        #print(r,d,lowest_alpha)
                        if not sdb.empty:
                            aux.append(case)

        return aux


norms = ['L2', 'H1']
dofs = ['v', 'p']
appended_data = []
set_index = False
for norm in norms:
    for dof in dofs:
        curve1 = Curve('sp_data.hdf5', norm)
        data = Data(curve1.data_base, norm, dof, set_index)
        appended_data.append(data.result_db)
        set_index = True

data = pd.concat(appended_data,axis=1)
print(data)
to_store = data.set_index('case')
to_store.to_excel('table.xlsx')
to_store.to_hdf('table.h5', key='data', format='table', data_columns=True)


