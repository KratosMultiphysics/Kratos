# import modules
import numpy as np
import os
import re
from scipy.stats import tmean, tstd, skew, kurtosis
from matplotlib.pylab import*
from scipy.spatial import distance
import json


class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class TimeHistoryData:

    def __init__(self, folder_name, file_name, data_index):

        with WorkFolderScope(folder_name):
            file_name = file_name + "_" + str(data_index) + ".dat"

            # read quantities
            with open(file_name, 'r') as data:
                data.seek(0)
                first_line = data.readline()
                second_line = data.readline()

            self.position = [float(item)
                            for item in (re.findall(r"[-+]?\d*\.\d+|\d+", first_line))]
            self.variable_names = re.findall('\w+', second_line)

            self.variable_results = {}
            for variable in range(len(self.variable_names)):
                self.variable_results[self.variable_names[variable]] = np.genfromtxt(
                    file_name, skip_header=2, skip_footer=0, usecols=(variable,))

            self.get_velocity_magnitude()


    def get_velocity_magnitude(self):
        vel_x = self.variable_results["TIME_AVERAGED_VELOCITY_X"]
        vel_y = self.variable_results["TIME_AVERAGED_VELOCITY_Y"]
        vel_z = self.variable_results["TIME_AVERAGED_VELOCITY_Z"]
        vel = np.sqrt(vel_x * vel_x + vel_y * vel_y + vel_z * vel_z)
        self.variable_results["TIME_AVERAGED_VELOCITY"] = vel


class PostProcessing:

    def __init__(self, case_name, case_folder_name, parameter_file):
        self.case_name = case_name
        self.result_folder_name = "."
        self.fig_folder_name = "."
        self.file_name = "time_averaged"
        self.ref_folder_name = "."
        self.ref_file_name = "time_accurate"

        with WorkFolderScope(case_folder_name):
            with open(parameter_file) as parameter_file:
                self.project_parameters = json.load(parameter_file)

        self.num_of_sampling_points = self.project_parameters["auxiliar_process_list"][0]["Parameters"]["sampling_points"]
        self.data = []
        self.ref_data = []
        for i in range(1, self.num_of_sampling_points + 1):
            self.data.append(TimeHistoryData(self.result_folder_name, self.file_name, i))
            self.ref_data.append(TimeHistoryData(self.ref_folder_name, self.ref_file_name, i))


    def get_space_data(self, t):
        dt = self.project_parameters["solver_settings"]["time_stepping"]["time_step"]
        time_index = int(t / dt) - 1
        velocity, ref_velocity = [], []
        velocity_x, ref_velocity_x = [], []
        velocity_y, ref_velocity_y = [], []
        pressure, ref_pressure = [], []
        for i in range(0, self.num_of_sampling_points):
            velocity.append(self.data[i].variable_results["TIME_AVERAGED_VELOCITY"][time_index])
            ref_velocity.append(self.ref_data[i].variable_results["TIME_AVERAGED_VELOCITY"][time_index])
            velocity_x.append(self.data[i].variable_results["TIME_AVERAGED_VELOCITY_X"][time_index])
            ref_velocity_x.append(self.ref_data[i].variable_results["TIME_AVERAGED_VELOCITY_X"][time_index])
            velocity_y.append(self.data[i].variable_results["TIME_AVERAGED_VELOCITY_Y"][time_index])
            ref_velocity_y.append(self.ref_data[i].variable_results["TIME_AVERAGED_VELOCITY_Y"][time_index])
            pressure.append(self.data[i].variable_results["TIME_AVERAGED_PRESSURE"][time_index])
            ref_pressure.append(self.ref_data[i].variable_results["TIME_AVERAGED_PRESSURE"][time_index])
        return velocity, ref_velocity, pressure, ref_pressure, velocity_x, ref_velocity_x,  velocity_y, ref_velocity_y


    def plot(self):
        from itertools import cycle
        cycol = cycle('bgrcmk')
        t_end = int(self.project_parameters["problem_data"]["end_time"]) + 1

        step = int(t_end / 5);
        for t in range(0, t_end, step):
            color = next(cycol)
            v, v_ref, p, p_ref, v_x, v_ref_x, v_y, v_ref_y = self.get_space_data(t)
            x = linspace(0, 10.0, self.num_of_sampling_points)

            plt.figure(1, figsize=(20,10))
            plt.plot(x, v, label='solved at time ' + str(t) + 's', c = color)
            plt.plot(x, v_ref, label='time_accurate at time ' + str(t) + 's', c = color, ls='--')
            plt.title(self.case_name + " Pipe Flow Mid-line velocity")
            plt.ylabel('TIME AVERAGED VELOCITY')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
            plt.grid(True)
            plt.legend()

            plt.figure(2, figsize=(20,10))
            plt.plot(x, v_x, label='solved at time ' + str(t) + 's', c = color)
            plt.plot(x, v_ref_x, label='time_accurate at time ' + str(t) + 's', c = color, ls='--')
            plt.title(self.case_name + " Pipe Flow Mid-line velocity x")
            plt.ylabel('TIME AVERAGED VELOCITY X')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
            plt.grid(True)
            plt.legend()

            plt.figure(3, figsize=(20,10))
            plt.plot(x, v_y, label='solved at time ' + str(t) + 's', c = color)
            plt.plot(x, v_ref_y, label='time_accurate at time ' + str(t) + 's', c = color, ls='--')
            plt.title(self.case_name + " Pipe Flow Mid-line velocity y")
            plt.ylabel('TIME AVERAGED VELOCITY Y')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
            plt.grid(True)
            plt.legend()

            plt.figure(4, figsize=(20,10))
            plt.plot(x, p, label='solved at time ' + str(t) + 's', c = color)
            plt.plot(x, p_ref, label='time_accurate ' + str(t) + 's',  c = color, ls='--')
            plt.title(self.case_name + " Pipe Flow Mid-line pressure")
            plt.ylabel('TIME AVERAGED PRESSURE')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(3, 3))
            plt.grid(True)
            plt.legend()

        # with WorkFolderScope(self.fig_folder_name):
        #     plt.figure(1, figsize=(20,10))
        #     savefig(self.case_name + '_velocity',dpi=600)

        #     plt.figure(2, figsize=(20,10))
        #     savefig(self.case_name + '_velocity_x',dpi=600)

        #     plt.figure(3, figsize=(20,10))
        #     savefig(self.case_name + '_velocity_y',dpi=600)

        #     plt.figure(4, figsize=(20,10))
        #     savefig(self.case_name + '_pressure',dpi=600)

        plt.show()


re100 = PostProcessing("Re100", ".","pipe_flow.json")
re100.plot()