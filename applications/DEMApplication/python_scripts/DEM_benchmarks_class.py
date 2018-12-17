from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *                                  # importing the Kratos Library
from KratosMultiphysics.DEMApplication import *
CheckForPreviousImport()                                          # check that KratosMultiphysics was imported in the main script
import shutil
from glob import glob
from math import pi, sin, cos, tan, atan, fabs

from os import system

def initialize_time_parameters(benchmark_number):

    number_of_coeffs_of_restitution = 1

    if benchmark_number==1:

        end_time                      = 0.0005
        dt                              = 6.4e-8 # Complies Rayleigh's condition
        graph_print_interval            = 0.000005
        number_of_points_in_the_graphic = 6

    elif benchmark_number==2:

        end_time                      = 0.007
        dt                              = 3e-7 # Complies Rayleigh's condition????????????????
        graph_print_interval            = 0.0001
        number_of_points_in_the_graphic = 6

    elif benchmark_number==3:

        end_time                      = 0.00031
        dt                              = 8.1e-9 #1.1e-9 # Complies Rayleigh's condition
        graph_print_interval            = 0.000001
        number_of_points_in_the_graphic = 6

    elif benchmark_number==4:

        end_time                      = 0.0002  #0.00003
        dt                              = 2e-8 #1.9e-9 # Complies Rayleigh's condition
        graph_print_interval            = 0.000001
        number_of_points_in_the_graphic = 17

    elif benchmark_number==5:

        end_time                      = 0.0000005
        dt                              = 3.6e-11  #3.6e-12 # Complies Rayleigh's condition
        graph_print_interval            = 0.00000005
        number_of_points_in_the_graphic = 17

    elif benchmark_number==6:

        end_time                      = 0.01
        dt                              = 1.0e-6  #1.0e-7 # Complies Rayleigh's condition ????????????????
        graph_print_interval            = 0.00025
        number_of_points_in_the_graphic = 17

    elif benchmark_number==7:

        end_time                      = 0.0005
        dt                              = 4.4614e-7 #4.4614e-8 # Complies Rayleigh's condition ????????????????
        graph_print_interval            = 0.000005
        number_of_points_in_the_graphic = 17

    elif benchmark_number==8:

        end_time                      = 0.02
        dt                              = 2.0e-6 #5.0e-7 # Complies Rayleigh's condition
        graph_print_interval            = 0.0001
        number_of_points_in_the_graphic = 17

    elif benchmark_number==9:

        end_time                      = 0.001 #0.0005
        dt                              = 5.0e-8 # 3.4e-8 # Complies Rayleigh's condition
        graph_print_interval            = 0.000005
        number_of_points_in_the_graphic = 6

    elif benchmark_number==10:

        end_time                      = 0.00015 #0.0005
        dt                              = 2.0e-8  #3.6e-12 # Complies Rayleigh's condition
        graph_print_interval            = 0.00001
        number_of_points_in_the_graphic = 10
        number_of_coeffs_of_restitution = 4

    elif benchmark_number==11:

        end_time                      = 0.00015 #0.0005
        dt                              = 1.0e-7 #3.6e-12 # Complies Rayleigh's condition
        graph_print_interval            = 0.00001
        number_of_points_in_the_graphic = 10
        number_of_coeffs_of_restitution = 4

    elif benchmark_number==12:

        end_time                      = 0.1
        dt                              = 5.0e-7
        graph_print_interval            = 1e-4
        number_of_points_in_the_graphic = 1

    elif benchmark_number==13:

        end_time                      = 2.0
        dt                              = 1.0e-4
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==14:

        end_time                      = 2.0
        dt                              = 1.0e-4
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==15:

        end_time                      = 2.0
        dt                              = 1.0e-4
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==16:

        end_time                      = 1.0
        dt                              = 0.50e-4
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==17:

        end_time                      = 1.0
        dt                              = 1.0e-6
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==20:          # Normal compression

        end_time                      = 0.01
        dt                              = 1e-5
        graph_print_interval            = 1e-5   # utilitzo com a output freq del grafic de punts
        number_of_points_in_the_graphic = 1

    elif benchmark_number==21:          # Normal compression with indentation

        end_time                      = 0.01
        dt                              = 1e-5
        graph_print_interval            = 1e-5
        number_of_points_in_the_graphic = 1

    elif benchmark_number==22:          # Tensile

        end_time                      = 0.05
        dt                              = 1e-5
        graph_print_interval            = 1e-5
        number_of_points_in_the_graphic = 1

    elif benchmark_number==23:          # Tensile with indentation

        end_time                      = 0.05
        dt                              = 1e-5
        graph_print_interval            = 1e-5
        number_of_points_in_the_graphic = 1

    elif benchmark_number==24:          # Shear

        end_time                      = 8e-5
        dt                              = 1e-7
        graph_print_interval            = 1e-7
        number_of_points_in_the_graphic = 1

    elif benchmark_number==25:          # Shear + radius expansion

        end_time                      = 8e-5
        dt                              = 1e-7
        graph_print_interval            = 1e-7
        number_of_points_in_the_graphic = 1

    elif benchmark_number==26:          #

        end_time                      = 0.1
        dt                              = 1e-5
        graph_print_interval            = 1e-4
        number_of_points_in_the_graphic = 1

    elif benchmark_number==27:          #UCS TEST

        end_time                      = 0.05
        dt                              = 5e-7
        graph_print_interval            = 5e-4
        number_of_points_in_the_graphic = 1

    elif benchmark_number==28:          #PENDULO3D . not ready

        end_time                      = 100
        dt                              = 1e-4
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==30:

        end_time                      = 0.5
        dt                              = 1.0e-3
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==31:

        end_time                      = 0.5
        dt                              = 1.0e-3
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==32:

        end_time                      = 0.5
        dt                              = 1.0e-6
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    elif benchmark_number==33:

        end_time                      = 0.5
        dt                              = 1.0e-6
        graph_print_interval            = 1e-2
        number_of_points_in_the_graphic = 1

    else: #benchmark_number==68:        #

        end_time                      = 1e-3
        dt                              = 1e-6
        graph_print_interval            = 1e-7
        number_of_points_in_the_graphic = 1

    return end_time, dt, graph_print_interval, number_of_points_in_the_graphic, number_of_coeffs_of_restitution

def PrintResultsMessage(test_number, it_is_success, error, elapsed_time, error_filename = 'errors.err'):
    with open(error_filename, 'a') as error_file:
        name = str(test_number)
        error_file.write('DEM Benchmark ' + name + ':')

        if it_is_success:
            error_file.write(' OK!........ Test ' + name + ' SUCCESSFUL (error: '
                             + str(round(error, 2)) + ', time: '
                             + str(round(elapsed_time, 2)) + 's.'')\n')
        else:
            error_file.write(' KO!........ Test ' + name + ' FAILED (error: ' + str(error) + ')\n')

class Benchmark1:

    def __init__(self):
        self.number = 1
        self.initial_normal_vel = 10.0

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        for node in modelpart.Nodes:
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X, -self.initial_normal_vel)
            else:
                node.SetSolutionStepValue(VELOCITY_X,  self.initial_normal_vel)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        pass

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        normal_contact_force_outfile_name = 'variables_for_node_1.txt'
        gnuplot_script_name = 'benchmark1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + normal_contact_force_outfile_name + "' every 20 u 1:8 w lp lt -1 lw 1.5 ps 1 pt 4")
        self.gnuplot_outfile.close()
        #print_gnuplot_files_on_screen(gnuplot_script_name)

        error1, error2, error3 = self.compute_errors(normal_contact_force_outfile_name)
        it_is_success = error1 < 1.0 and error2 < 1.0 and error3 < 1.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def compute_errors(self, normal_contact_force_outfile_name):

        Chung_data = []; DEM_data = []

        with open('paper_data/benchmark1_graph1.dat') as inf:
            for line in inf:
                Chung_data.append(float(line))

        with open(normal_contact_force_outfile_name) as inf:
            for line in inf:
                parts = line.split()
                if parts[0] == '#Time':
                    break
            for line in inf:
                parts = line.split()
                DEM_data.append(float(parts[7]))

        error = fabs(max(DEM_data) - float(Chung_data[0]))/float(Chung_data[0])

        print("Error in restitution numbers =", 100*error,"%")

        error1 = 100*error

        error2 = error3 = 0

        return error1, error2, error3

class Benchmark2:

    def __init__(self):
        self.number = 2
        self.initial_normal_vel = -0.2

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)


    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        pass

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        normal_contact_force_outfile_name = 'variables_for_node_2.txt'
        gnuplot_script_name = 'benchmark2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + normal_contact_force_outfile_name + "' every 10 u 1:10 w lp lt 3 lw 1.5 ps 1 pt 6")
        self.gnuplot_outfile.close()
        #print_gnuplot_files_on_screen(gnuplot_script_name)

        error1, error2, error3 = self.compute_errors(normal_contact_force_outfile_name)
        it_is_success = error1 < 1.0 and error2 < 1.0 and error3 < 1.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def compute_errors(self, normal_contact_force_outfile_name):

        Chung_data = []; DEM_data = []

        with open('paper_data/benchmark2_graph1.dat') as inf:
            for line in inf:
                Chung_data.append(float(line))

        with open(normal_contact_force_outfile_name) as inf:
            for line in inf:
                parts = line.split()
                if parts[0] == '#Time':
                    break
            for line in inf:
                parts = line.split()
                DEM_data.append(float(parts[9]))

        error = fabs(max(DEM_data) - float(Chung_data[0]))/float(Chung_data[0])

        print("Error in restitution numbers =", 100*error,"%")

        error1 = 100*error

        error2 = error3 = 0

        return error1, error2, error3

class Benchmark3:

    def __init__(self):
        self.number = 3
        self.restitution_numbers_list = []
        self.initial_normal_vel = 0
        self.generated_data = None

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        #number = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)

        if number_of_points_in_the_graphic == 1:
            number = 0
        else:
            number = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)

        for node in modelpart.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(VELOCITY_Z)
            modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = number

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        for node in modelpart.Nodes:
            final_vel = node.GetSolutionStepValue(VELOCITY_Z)

        restitution_coefficient = -final_vel / self.initial_normal_vel
        self.restitution_numbers_list.append(restitution_coefficient)

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.output_filename = "benchmark3_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat'
        self.generated_data = open(self.output_filename, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            first_col = 1/(number_of_points_in_the_graphic-1) * i
            self.generated_data.write("%6.4f %11.8f" % (first_col, self.restitution_numbers_list[i]) + '\n')
        self.generated_data.close()

        gnuplot_script_name = 'benchmark3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + self.output_filename + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4, '"\
                                                      + self.output_filename + "' u 1:3 w lp lt 2 lw 1.5 ps 2 pt 6")
        self.gnuplot_outfile.close()

        self.create_gnuplot_scripts(self.output_filename, dt)

        error1, error2, error3 = self.compute_errors(self.output_filename)
        it_is_success = error1 < 1.0 and error2 < 1.0 and error3 < 1.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, output_filename, dt):

        gnuplot_script_name_1 = 'benchmark3_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Coefficient of restitution'\nset ylabel 'Damping ratio'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:1][0:1] '" + output_filename + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark3_graph1.dat' w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark3_graph1.dat' w lp ls 2 t 'Cast iron'\n")
        self.gnuplot_outfile.close()

        #print_gnuplot_files_on_screen(gnuplot_script_name_1)

    def compute_errors(self, output_filename):

        lines_Chung = lines_DEM = list(range(0, 6))
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark3_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split()
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0

        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        generated_data_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_Chung_data

        print("Error in restitution numbers =", 100*generated_data_error,"%")

        error1 = 100*generated_data_error

        error2 = error3 = 0

        return error1, error2, error3


class Benchmark4:

    def __init__(self):
        self.number = 4
        self.initial_module_vel = 3.9
        self.initial_tangential_vel = 0
        self.radius = 0.0025
        self.degrees = 0
        self.angles_list = []
        self.tangential_restitution_coefficient_list = []
        self.final_angular_vel_list = []
        self.rebound_angle_list = []
        self.final_angular_vel_list_outfile = None
        self.rebound_angle_list_outfile = None
        self.tangential_restitution_coefficient_list_outfile = None

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  self.initial_module_vel * sin(self.degrees * pi / 180.0)
        initial_normal_vel = -self.initial_module_vel * cos(self.degrees * pi / 180.0)

        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.initial_tangential_vel)
            node.SetSolutionStepValue(VELOCITY_Z, initial_normal_vel)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        for node in modelpart.Nodes:

            final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_X)
            final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Y)
            final_normal_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
            final_tangential_contact_velocity = final_tangential_center_velocity + final_angular_vel * self.radius
            rebound_angle = 180 / pi * atan(final_tangential_contact_velocity / final_normal_center_velocity)
            tangential_restitution_coefficient = final_tangential_center_velocity / self.initial_tangential_vel

        self.final_angular_vel_list.append(final_angular_vel)
        self.rebound_angle_list.append(rebound_angle)
        self.tangential_restitution_coefficient_list.append(tangential_restitution_coefficient)
        self.angles_list.append(self.degrees)

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.tangential_restitution_coefficient_list_outfile_name = "benchmark4_dt_" + str(dt) + '_tangential_restitution_coefficient_list_data.dat'
        self.final_angular_vel_list_outfile_name = "benchmark4_dt_" + str(dt) + '_final_angular_vel_list_data.dat'
        self.rebound_angle_list_outfile_name = "benchmark4_dt_" + str(dt) + '_rebound_angle_list_data.dat'
        self.tangential_restitution_coefficient_list_outfile = open(self.tangential_restitution_coefficient_list_outfile_name, 'w')
        self.final_angular_vel_list_outfile = open(self.final_angular_vel_list_outfile_name, 'w')
        self.rebound_angle_list_outfile = open(self.rebound_angle_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.tangential_restitution_coefficient_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.tangential_restitution_coefficient_list[i]) + '\n')
            self.final_angular_vel_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.final_angular_vel_list[i]) + '\n')
            self.rebound_angle_list_outfile.write("%14.8f %14.8f" % (self.angles_list[i], self.rebound_angle_list[i]) + '\n')
        self.tangential_restitution_coefficient_list_outfile.close()
        self.final_angular_vel_list_outfile.close()
        self.rebound_angle_list_outfile.close()

        self.create_gnuplot_scripts(self.tangential_restitution_coefficient_list_outfile_name, self.final_angular_vel_list_outfile_name,\
                                    self.rebound_angle_list_outfile_name, dt)

        error1, error2, error3 = self.compute_errors(self.tangential_restitution_coefficient_list_outfile_name, self.final_angular_vel_list_outfile_name,\
                                    self.rebound_angle_list_outfile_name)
        it_is_success = error1 < 2.0 and error2 < 2.0 and error3 < 2.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, tangential_restitution_coefficient_list_outfile_name, final_angular_vel_list_outfile_name,\
                               rebound_angle_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark4_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][.4:1] '" + tangential_restitution_coefficient_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph1.dat' index 0 w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph1.dat' index 1 w lp ls 2 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph1.dat' index 2 w p pt 7 ps 2 lt -1 t 'Experimental'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark4_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Incident angle (deg)'\nset ylabel 'Final angular velocity (rad/s)'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][-750:0] '" + final_angular_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph2.dat' index 0 w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph2.dat' index 1 w lp ls 2 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph2.dat' index 2 w p pt 7 ps 2 lt -1 t 'Experimental'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_3 = 'benchmark4_comparison_3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_3, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Incident angle (deg)'\nset ylabel 'Rebound angle (deg)'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:90][-30:90] '" + rebound_angle_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph3.dat' index 0 w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph3.dat' index 1 w lp ls 2 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark4_graph3.dat' index 2 w p pt 7 ps 2 lt -1 t 'Experimental'\n")
        self.gnuplot_outfile.close()
        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        print_gnuplot_files_on_screen(gnuplot_script_name_3)'''

    def compute_errors(self, tangential_restitution_coefficient_list_outfile_name, final_angular_vel_list_outfile_name, rebound_angle_list_outfile_name):

        lines_Chung = list(range(17, 30)); lines_DEM = list(range(0, 8)) + list(range(9, 16, 2)) + [16]
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark4_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(tangential_restitution_coefficient_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_tangential_restitution_coefficient_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_tangential_restitution_coefficient_error+=fabs(i-j)
        final_tangential_restitution_coefficient_error/=summation_of_Chung_data
        print("Error in tangential restitution coefficient =", 100*final_tangential_restitution_coefficient_error,"%")

        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark4_graph2.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(final_angular_vel_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_angular_vel_total_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_angular_vel_total_error+=fabs(i-j)
        final_angular_vel_total_error/=summation_of_Chung_data
        print("Error in final angular vel =", 100*final_angular_vel_total_error,"%")

        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark4_graph3.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(rebound_angle_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_rebound_angle_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_rebound_angle_error+=fabs(i-j)
        final_rebound_angle_error/=summation_of_Chung_data
        print("Error in final rebound angle =", 100*final_rebound_angle_error,"%")

        error1 = 100*final_tangential_restitution_coefficient_error
        error2 = 100*final_angular_vel_total_error
        error3 = 100*final_rebound_angle_error

        return error1, error2, error3


class Benchmark5:

    def __init__(self):
        self.number = 5
        self.initial_normal_vel = -5.0
        self.initial_tangential_vel = 0
        self.radius = 0.00001
        self.Vst_div_mu_per_Vcn_list = []
        self.Vst_prima_div_mu_per_Vcn_prima_list = []
        self.r_w1_prima_div_mu_per_Vcn_list = []
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile = None
        self.r_w1_prima_div_mu_per_Vcn_list_outfile = None

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  -self.initial_normal_vel * tan(degrees * pi / 180.0)

        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.initial_tangential_vel)
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        mu = 0.3

        for node in modelpart.Nodes:
            final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_X)
            final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Y)
            final_normal_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
            Vst_div_mu_per_Vcn = -self.initial_tangential_vel / (mu * self.initial_normal_vel)
            Vst_prima_div_mu_per_Vcn_prima = (final_tangential_center_velocity + final_angular_vel * self.radius) / (mu * final_normal_center_velocity)
            r_w1_prima_div_mu_per_Vcn = -self.radius * final_angular_vel / (mu * self.initial_normal_vel)

        self.Vst_div_mu_per_Vcn_list.append(Vst_div_mu_per_Vcn)
        self.Vst_prima_div_mu_per_Vcn_prima_list.append(Vst_prima_div_mu_per_Vcn_prima)
        self.r_w1_prima_div_mu_per_Vcn_list.append(r_w1_prima_div_mu_per_Vcn)

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name = "benchmark5_dt_" + str(dt) + '_Vst_prima_div_mu_per_Vcn_prima_list_data.dat'
        self.r_w1_prima_div_mu_per_Vcn_list_outfile_name = "benchmark5_dt_" + str(dt) + '_r_w1_prima_div_mu_per_Vcn_list_data.dat'
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile = open(self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, 'w')
        self.r_w1_prima_div_mu_per_Vcn_list_outfile = open(self.r_w1_prima_div_mu_per_Vcn_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.Vst_prima_div_mu_per_Vcn_prima_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_mu_per_Vcn_list[i], self.Vst_prima_div_mu_per_Vcn_prima_list[i]) + '\n')
            self.r_w1_prima_div_mu_per_Vcn_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_mu_per_Vcn_list[i], self.r_w1_prima_div_mu_per_Vcn_list[i]) + '\n')
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile.close()
        self.r_w1_prima_div_mu_per_Vcn_list_outfile.close()

        self.create_gnuplot_scripts(self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, self.r_w1_prima_div_mu_per_Vcn_list_outfile_name, dt)

        error1, error2, error3 = self.compute_errors(self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, self.r_w1_prima_div_mu_per_Vcn_list_outfile_name)
        it_is_success = error1 < 2.0 and error2 < 2.0 and error3 < 2.0
        error_measure = error1 + error2 + error3
        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, r_w1_prima_div_mu_per_Vcn_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark5_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:14][-4:6] '" + Vst_prima_div_mu_per_Vcn_prima_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph1.dat' index 0 w lp ls 1 t 'Steel',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph1.dat' index 1 w lp ls 2 t 'Polyethylene',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph1.dat' index 2 w p pt 7 ps 2 lt -1 t 'FEM'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark5_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Normalized final angular velocity'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:20][-6:0] '" + r_w1_prima_div_mu_per_Vcn_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph2.dat' index 0 w lp ls 1 t 'Steel',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph2.dat' index 1 w lp ls 2 t 'Polyethylene',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph2.dat' index 2 w p pt 7 ps 2 lt -1 t 'FEM'\n")
        self.gnuplot_outfile.close()
        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)'''

    def compute_errors(self, Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, r_w1_prima_div_mu_per_Vcn_list_outfile_name):

        lines_Chung = list(range(49, 53)); lines_DEM = list(range(11, 15)) # Sliding regime for the time being
        #lines_Chung = list(range(38, 53)); lines_DEM = list(range(0, 15)) # Whole diagram
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark5_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(Vst_prima_div_mu_per_Vcn_prima_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_Vst_prima_div_mu_per_Vcn_prima_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_Vst_prima_div_mu_per_Vcn_prima_error+=fabs(i-j)

        final_Vst_prima_div_mu_per_Vcn_prima_error/=summation_of_Chung_data

        print("Error in final Vst prima div mu per Vcn prima =", 100*final_Vst_prima_div_mu_per_Vcn_prima_error,"%")

        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark5_graph2.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(r_w1_prima_div_mu_per_Vcn_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_r_w1_prima_div_mu_per_Vcn_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_r_w1_prima_div_mu_per_Vcn_error+=fabs(i-j)

        final_r_w1_prima_div_mu_per_Vcn_error/=summation_of_Chung_data
        print("Error in final r w1 prima div mu per Vcn =", 100*final_r_w1_prima_div_mu_per_Vcn_error,"%")

        error1 = 100*final_Vst_prima_div_mu_per_Vcn_prima_error
        error2 = 100*final_r_w1_prima_div_mu_per_Vcn_error
        error3 = 0

        return error1, error2, error3


class Benchmark6:

    def __init__(self):
        self.number = 6
        self.initial_normal_vel = -0.2
        self.initial_tangential_vel = 0
        self.radius = 0.1
        self.special_quantity_list = []
        self.beta_list = []
        self.Vst_div_Vcn_list = []
        self.Vst_prima_div_Vcn_prima_list = []
        self.beta_list_outfile = None
        self.Vst_prima_div_Vcn_prima_list_outfile = None

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel = -self.initial_normal_vel * tan(degrees * pi / 180.0) # Here is tangential of the contact point, only. In X axis
        initial_angular_vel = -self.initial_tangential_vel / self.radius # In Y axis

        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)
            node.SetSolutionStepValue(ANGULAR_VELOCITY_Y, initial_angular_vel)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        mu = 0.4
        restitution_coeff = 0.5

        for node in modelpart.Nodes:
            special_quantity = -3.5 * mu * (1.0 + restitution_coeff) * self.initial_normal_vel / self.initial_tangential_vel
            final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)
            final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_X)
            final_normal_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
            beta = -(final_tangential_center_velocity - final_angular_vel * self.radius)/ self.initial_tangential_vel
            Vst_div_Vcn = -self.initial_tangential_vel / self.initial_normal_vel
            Vst_prima_div_Vcn_prima = (final_tangential_center_velocity - final_angular_vel * self.radius) / final_normal_center_velocity

        self.special_quantity_list.append(special_quantity)
        self.beta_list.append(beta)
        self.Vst_div_Vcn_list.append(Vst_div_Vcn)
        self.Vst_prima_div_Vcn_prima_list.append(Vst_prima_div_Vcn_prima)

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.beta_list_outfile_name = "benchmark6_dt_" + str(dt) + '_beta_list_data.dat'
        self.Vst_prima_div_Vcn_prima_list_outfile_name = "benchmark6_dt_" + str(dt) + '_Vst_prima_div_Vcn_prima_data.dat'
        self.beta_list_outfile = open(self.beta_list_outfile_name, 'w')
        self.Vst_prima_div_Vcn_prima_list_outfile = open(self.Vst_prima_div_Vcn_prima_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.beta_list_outfile.write("%14.8f %14.8f" % (self.special_quantity_list[i], self.beta_list[i]) + '\n')
            self.Vst_prima_div_Vcn_prima_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_Vcn_list[i], self.Vst_prima_div_Vcn_prima_list[i]) + '\n')
        self.beta_list_outfile.close()
        self.Vst_prima_div_Vcn_prima_list_outfile.close()

        self.create_gnuplot_scripts(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name, dt)

        error1, error2, error3 = self.compute_errors(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name)
        it_is_success = error1 < 3.0 and error2 < 3.0 and error3 < 3.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, beta_list_outfile_name, Vst_prima_div_Vcn_prima_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark6_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:25][-1:.6] '" + beta_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark6_graph1.dat' index 0 w lp ls 1 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark6_graph1.dat' index 1 w lp ls 2 t 'Nylon'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark6_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Tangent of incident angle'\nset ylabel 'Tangent of recoil angle'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:7][-2:8] '" + Vst_prima_div_Vcn_prima_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark6_graph2.dat' index 0 w lp ls 1 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark6_graph2.dat' index 1 w lp ls 2 t 'Nylon'\n")
        self.gnuplot_outfile.close()
        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)'''

    def compute_errors(self, beta_list_outfile_name, Vst_prima_div_Vcn_prima_list_outfile_name):

        lines_Chung = list(range(1, 7)); lines_DEM = list(range(16, 10, -1)) # Sliding regime for the time being
        #lines_Chung = list(range(1, 17)); lines_DEM = list(range(0, 16)) # Whole diagram
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark6_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0

        with open(beta_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_beta_list_outfile_name_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        DEM_data.reverse()
        for i, j in zip(DEM_data, Chung_data):
            final_beta_list_outfile_name_error+=fabs(i-j)

        final_beta_list_outfile_name_error/=summation_of_Chung_data
        print("Error in final beta =", 100*final_beta_list_outfile_name_error,"%")

        lines_Chung = list(range(13, 17)); lines_DEM = list(range(12, 16)) # Sliding regime for the time being
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark6_graph2.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(Vst_prima_div_Vcn_prima_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_Vst_prima_div_Vcn_prima_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)


        for i, j in zip(DEM_data, Chung_data):
            final_Vst_prima_div_Vcn_prima_error+=fabs(i-j)

        final_Vst_prima_div_Vcn_prima_error/=summation_of_Chung_data
        print("Error in final Vst prima div Vcn =", 100*final_Vst_prima_div_Vcn_prima_error,"%")

        error1 = 100*final_beta_list_outfile_name_error
        error2 = 100*final_Vst_prima_div_Vcn_prima_error
        error3 = 0

        return error1, error2, error3

class Benchmark7:

    def __init__(self):
        self.number = 7
        self.initial_angular_vel = 0
        self.final_tangential_center_vel_list_outfile = None
        self.final_angular_vel_list_outfile = None
        self.initial_angular_vel_list = []
        self.final_tangential_center_vel_list = []
        self.final_angular_vel_list = []

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        initial_normal_vel = 0.2
        radius = 0.1
        degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_angular_vel =  initial_normal_vel / radius * tan(degrees * pi / 180.0) # Here is tangential of the contact point, only

        for node in modelpart.Nodes:
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X,  initial_normal_vel)
                node.SetSolutionStepValue(ANGULAR_VELOCITY_Y,  self.initial_angular_vel)
            else:
                node.SetSolutionStepValue(VELOCITY_X, -initial_normal_vel)
                node.SetSolutionStepValue(ANGULAR_VELOCITY_Y, -self.initial_angular_vel)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        for node in modelpart.Nodes:
            if node.Id == 1:
                final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
                final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)

        self.initial_angular_vel_list.append(self.initial_angular_vel)
        self.final_tangential_center_vel_list.append(final_tangential_center_velocity)
        self.final_angular_vel_list.append(final_angular_vel)

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.final_tangential_center_vel_list_outfile_name = "benchmark7_dt_" + str(dt) + '_final_tangential_center_vel_list_data.dat'
        self.final_angular_vel_list_outfile_name = "benchmark7_dt_" + str(dt) + '_final_angular_vel_list_data.dat'
        self.final_tangential_center_vel_list_outfile = open(self.final_tangential_center_vel_list_outfile_name, 'w')
        self.final_angular_vel_list_outfile = open(self.final_angular_vel_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.final_tangential_center_vel_list_outfile.write("%14.8f %14.8f" % (self.initial_angular_vel_list[i], self.final_tangential_center_vel_list[i]) + '\n')
            self.final_angular_vel_list_outfile.write("%14.8f %14.8f" % (self.initial_angular_vel_list[i], self.final_angular_vel_list[i]) + '\n')
        self.final_tangential_center_vel_list_outfile.close()
        self.final_angular_vel_list_outfile.close()

        gnuplot_script_name = 'benchmark7_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set multiplot layout 2, 1; set grid; set bmargin 0; set format x \"\"; set ytics -5, 5; set key bottom;\
                                    plot [0:25][-10:10] '" + self.final_tangential_center_vel_list_outfile_name + "' w lp lw 1.5 ps 2 pt 4;\
                                    set bmargin; set tmargin 0; set format x \"%g\"; set ytics 0, 5, 20; set key top;\
                                    plot [0:25][0:25] '" + self.final_angular_vel_list_outfile_name + "' w lp lw 1.5 lt 3 ps 2 pt 6; unset multiplot")
        self.gnuplot_outfile.close()

        self.create_gnuplot_scripts(self.final_tangential_center_vel_list_outfile_name, self.final_angular_vel_list_outfile_name, dt)

        error1, error2, error3 = self.compute_errors(self.final_tangential_center_vel_list_outfile_name, self.final_angular_vel_list_outfile_name)
        it_is_success = error1 < 1.0 and error2 < 1.0 and error3 < 1.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, final_tangential_center_vel_list_outfile_name, final_angular_vel_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark7_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:25][-10:10] '" + final_tangential_center_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark7_graph1.dat' w lp ls 1 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark7_graph1.dat' w lp ls 2 t 'Copper'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark7_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Initial angular velocity (rad/s)'\nset ylabel 'Final angular velocity (rad/s)'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:25][0:25] '" + final_angular_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark7_graph2.dat' w lp ls 1 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark7_graph2.dat' w lp ls 2 t 'Copper'\n")
        self.gnuplot_outfile.close()
        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)'''

    def compute_errors(self, final_tangential_center_vel_list_outfile_name, final_angular_vel_list_outfile_name):

        lines_Chung = []; lines_DEM = []; lines_Chung = list(range(0, 17)); lines_DEM = list(range(0, 17))
        Chung_data = []; DEM_data = []
        i = 0
        with open('paper_data/benchmark7_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split()
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(final_tangential_center_vel_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_tangential_center_vel_error = 0

        for i, j in zip(DEM_data, Chung_data):
            final_tangential_center_vel_error+=fabs(i-j)
        print("Error in final tangential center vel =", final_tangential_center_vel_error)

        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark7_graph2.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split()
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(final_angular_vel_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_angular_vel_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_angular_vel_error+=fabs(i-j)

        final_angular_vel_error/=summation_of_Chung_data
        print("Error in final angular vel =", 100*final_angular_vel_error,"%")

        error1 = 100*final_tangential_center_vel_error
        error2 = 100*final_angular_vel_error
        error3 = 0

        return error1, error2, error3


class Benchmark8:

    def __init__(self):
        self.number = 8
        self.initial_normal_vel = 0.2
        self.initial_tangential_vel = 0
        self.radius = 0.1
        self.special_quantity_list = []
        self.beta_list = []
        self.Vst_div_Vcn_list = []
        self.Vst_prima_div_Vcn_prima_list = []
        self.beta_list_outfile = None
        self.Vst_prima_div_Vcn_prima_list_outfile = None

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        degrees = 90 - 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  self.initial_normal_vel * tan(degrees * pi / 180.0) # Here is tangential of the contact point, only
        initial_angular_vel    =  -self.initial_tangential_vel / self.radius

        for node in modelpart.Nodes:
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X, self.initial_normal_vel)
                node.SetSolutionStepValue(ANGULAR_VELOCITY_Y, initial_angular_vel)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        mu = 0.4
        restitution_coeff = 0.5

        for node in modelpart.Nodes:
            if node.Id == 1:
                special_quantity = 3.5 * mu * (1.0 + restitution_coeff) * self.initial_normal_vel / self.initial_tangential_vel
                final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)
                final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
                final_normal_center_velocity = node.GetSolutionStepValue(VELOCITY_X)
                beta = -(final_tangential_center_velocity - final_angular_vel * self.radius)/ self.initial_tangential_vel
                Vst_div_Vcn = self.initial_tangential_vel / self.initial_normal_vel
                Vst_prima_div_Vcn_prima = -(final_tangential_center_velocity - final_angular_vel * self.radius) / final_normal_center_velocity

        self.special_quantity_list.append(special_quantity)
        self.beta_list.append(beta)
        self.Vst_div_Vcn_list.append(Vst_div_Vcn)
        self.Vst_prima_div_Vcn_prima_list.append(Vst_prima_div_Vcn_prima)

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.beta_list_outfile_name = 'benchmark8_dt_' + str(dt) + 's_beta_list_data.dat'
        self.Vst_prima_div_Vcn_prima_list_outfile_name = 'benchmark8_dt_' + str(dt) + 's_Vst_prima_div_Vcn_prima_list_data.dat'
        self.beta_list_outfile = open(self.beta_list_outfile_name, 'w')
        self.Vst_prima_div_Vcn_prima_list_outfile = open(self.Vst_prima_div_Vcn_prima_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.beta_list_outfile.write("%14.8f %14.8f" % (self.special_quantity_list[i], self.beta_list[i]) + '\n')
            self.Vst_prima_div_Vcn_prima_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_Vcn_list[i], self.Vst_prima_div_Vcn_prima_list[i]) + '\n')

        self.beta_list_outfile.close()
        self.Vst_prima_div_Vcn_prima_list_outfile.close()

        self.create_gnuplot_scripts(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name, dt)

        error1, error2, error3 = self.compute_errors(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name)
        it_is_success = error1 < 3.0 and error2 < 3.0 and error3 < 3.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, beta_list_outfile_name, Vst_prima_div_Vcn_prima_list_outfile_name, dt):

        gnuplot_script_name_1 = 'benchmark8_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:25][-1:.6] '" + beta_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark8_graph1.dat' index 0 w lp ls 1 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark8_graph1.dat' index 1 w lp ls 2 t 'Nylon'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark8_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Tangent of incident angle'\nset ylabel 'Tangent of recoil angle'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:8][-2:8] '" + Vst_prima_div_Vcn_prima_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark8_graph2.dat' index 0 w lp ls 1 t 'Al. alloy',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark8_graph2.dat' index 1 w lp ls 2 t 'Nylon'\n")
        self.gnuplot_outfile.close()
        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)'''

    def compute_errors(self, beta_list_outfile_name, Vst_prima_div_Vcn_prima_list_outfile_name):

        lines_Chung = []; lines_DEM = []; lines_Chung = list(range(1, 7)); lines_DEM = list(range(0, 6)) # Sliding regime for the time being
        #lines_Chung = []; lines_DEM = []; lines_Chung = list(range(1, 18)); lines_DEM = list(range(0, 17)) # Whole diagram
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark8_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(beta_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_beta_list_outfile_name_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            final_beta_list_outfile_name_error+=fabs(i-j)

        final_beta_list_outfile_name_error/=summation_of_Chung_data
        print("Error in final beta =", 100*final_beta_list_outfile_name_error,"%")

        lines_Chung = []; lines_DEM = []; lines_DEM = list(range(4, 0, -1)); lines_Chung = list(range(13, 17)) # Sliding regime for the time being
        #lines_Chung = list(range(1, 17)); lines_DEM = list(range(0, 16)) # Whole diagram

        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark8_graph2.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split(',')
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(Vst_prima_div_Vcn_prima_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_Vst_prima_div_Vcn_prima_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        DEM_data.reverse()
        for i, j in zip(DEM_data, Chung_data):
            final_Vst_prima_div_Vcn_prima_error+=fabs(i-j)

        final_Vst_prima_div_Vcn_prima_error/=summation_of_Chung_data
        print("Error in final Vst prima div Vcn =", 100*final_Vst_prima_div_Vcn_prima_error,"%")

        error1 = 100*final_beta_list_outfile_name_error
        error2 = 100*final_Vst_prima_div_Vcn_prima_error
        error3 = 0

        return error1, error2, error3

class Benchmark9:

    def __init__(self):
        self.number = 9
        self.initial_normal_vel = 200.0
        self.restitution_numbers_list = []
        self.generated_data = None

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        if number_of_points_in_the_graphic == 1:
            number = 0
        else:
            number = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)

        for node in modelpart.Nodes:

            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X,  self.initial_normal_vel)
                node.SetSolutionStepValue(VELOCITY_Z, 0.0)
                modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = number
            else:
                node.SetSolutionStepValue(VELOCITY_X, -self.initial_normal_vel)
                node.SetSolutionStepValue(VELOCITY_Z, 0.0)
                modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = number

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        for node in modelpart.Nodes:
            if node.Id == 1:
                final_vel = node.GetSolutionStepValue(VELOCITY_X)

        restitution_coefficient = -final_vel / self.initial_normal_vel
        self.restitution_numbers_list.append(restitution_coefficient)

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.output_filename = "benchmark9_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat'
        self.generated_data = open(self.output_filename, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            if number_of_points_in_the_graphic == 1:
                first_col = 0
            else:
                first_col = 1/(number_of_points_in_the_graphic-1) * i
            self.generated_data.write("%6.4f %11.8f" % (first_col, self.restitution_numbers_list[i]) + '\n')
        self.generated_data.close()

        gnuplot_script_name = 'benchmark9_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + self.output_filename + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4, '"\
                                                      + self.output_filename + "' u 1:3 w lp lt 2 lw 1.5 ps 2 pt 6")
        self.gnuplot_outfile.close()

        self.create_gnuplot_scripts(self.output_filename, dt)

        error1, error2, error3 = self.compute_errors(self.output_filename)
        it_is_success = error1 < 1.0 and error2 < 1.0 and error3 < 1.0
        error_measure = error1 + error2 + error3

        PrintResultsMessage(self.number, it_is_success, error_measure, elapsed_time)

    def create_gnuplot_scripts(self, output_filename, dt):

        gnuplot_script_name_1 = 'benchmark9_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Coefficient of restitution'\nset ylabel 'Damping ratio'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:1][0:1] '" + output_filename + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark9_graph1.dat' w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark9_graph1.dat' w lp ls 2 t 'Cast iron'\n")
        self.gnuplot_outfile.close()

        #print_gnuplot_files_on_screen(gnuplot_script_name_1)

    def compute_errors(self, output_filename):

        lines_Chung = lines_DEM = list(range(0, 6));
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0
        i = 0
        with open('paper_data/benchmark9_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split()
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        generated_data_error = 0

        for j in Chung_data:
            summation_of_Chung_data+=abs(j)

        for i, j in zip(DEM_data, Chung_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_Chung_data

        print("Error in restitution numbers =", 100*generated_data_error,"%")

        error1 = 100*generated_data_error

        error2 = error3 = 0

        return error1, error2, error3


class Benchmark10: ########## LINEAR THORNTON

    def __init__(self):
        self.number = 10
        self.initial_normal_vel = -5.0
        self.initial_tangential_vel = 0
        self.radius = 0.025
        self.normalized_impact_angle_list = []
        self.normalized_rebound_tangential_surface_vel_list = []
        self.normalized_rebound_angular_velocity_list = []
        self.tangential_coefficient_of_restitution_list = []
        self.normalized_rebound_tangential_surface_vel_list_outfile = None
        self.normalized_rebound_angular_velocity_list_outfile = None
        self.tangential_coefficient_of_restitution_list_outfile = None
        self.coeff_of_restitution = -1.0
        self.coeff_of_rest_string = None
        self.lines_Thornton = []
        self.lines_DEM = []
        self.degrees = 0

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):   # Change this function name from 'set_initial_data' to 'set_initial_data'

        if iteration == 1:
            self.degrees = 1
        else:
            self.degrees = 50 * (iteration - 1)/number_of_points_in_the_graphic

        if coeff_of_restitution_iteration==1:
            self.coeff_of_restitution=0.25
            self.coeff_of_rest_string='025'
            self.lines_Thornton = [12, 13, 15, 16, 18, 19]
            self.lines_DEM = [0, 1, 3, 4, 5, 6]
        elif coeff_of_restitution_iteration==2:
            self.coeff_of_restitution=0.50
            self.coeff_of_rest_string='050'
            self.lines_Thornton = [14, 15, 17, 18, 20, 22, 23]
            self.lines_DEM = [0, 1, 3, 4, 5, 6, 7]
        elif coeff_of_restitution_iteration==3:
            self.coeff_of_restitution=0.75
            self.coeff_of_rest_string='075'
            self.lines_Thornton = [14, 15, 17, 18, 19, 22, 23, 24]
            self.lines_DEM = [0, 1, 3, 4, 5, 6, 7, 8]
        else:
            self.coeff_of_restitution=0.90
            self.coeff_of_rest_string='090'
            self.lines_Thornton = [13, 14, 16, 17, 18, 21, 22, 23]
            self.lines_DEM = [0, 1, 3, 4, 5, 6, 7, 8]

        self.initial_tangential_vel = -self.initial_normal_vel * tan(self.degrees * pi / 180.0)

        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.initial_tangential_vel)
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)
            modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = self.coeff_of_restitution

        print(self.coeff_of_restitution)

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        mu = 0.1

        for node in modelpart.Nodes:
            final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_X)
            final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Y)
            normalized_impact_angle = 2.0 * tan(self.degrees * pi / 180.0) / (mu * (1 + self.coeff_of_restitution))
            normalized_rebound_tangential_surface_vel = -2.0 * (final_tangential_center_velocity + final_angular_vel * self.radius) / (self.initial_normal_vel * mu * (1 + self.coeff_of_restitution))
            normalized_rebound_angular_velocity = -2.0 * self.radius * final_angular_vel / (self.initial_normal_vel * mu * (1 + self.coeff_of_restitution))
            tangential_coefficient_of_restitution = 5.0/7.0 + 2.0 * normalized_rebound_tangential_surface_vel / (7.0 * normalized_impact_angle)

        self.normalized_impact_angle_list.append(normalized_impact_angle)
        self.normalized_rebound_tangential_surface_vel_list.append(normalized_rebound_tangential_surface_vel)
        self.normalized_rebound_angular_velocity_list.append(normalized_rebound_angular_velocity)
        self.tangential_coefficient_of_restitution_list.append(tangential_coefficient_of_restitution)

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.normalized_rebound_tangential_surface_vel_list_outfile_name = "benchmark10_dt_" + str(dt) + '_normalized_rebound_tangential_surface_vel_list_data.dat'
        self.normalized_rebound_angular_velocity_list_outfile_name = "benchmark10_dt_" + str(dt) + '_normalized_rebound_angular_velocity_list_data.dat'
        self.tangential_coefficient_of_restitution_list_outfile_name = "benchmark10_dt_" + str(dt) + '_tangential_coefficient_of_restitution_list_data.dat'

        self.normalized_rebound_tangential_surface_vel_list_outfile = open(self.normalized_rebound_tangential_surface_vel_list_outfile_name, 'w')
        self.normalized_rebound_angular_velocity_list_outfile = open(self.normalized_rebound_angular_velocity_list_outfile_name, 'w')
        self.tangential_coefficient_of_restitution_list_outfile = open(self.tangential_coefficient_of_restitution_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.normalized_rebound_tangential_surface_vel_list_outfile.write("%14.8f %14.8f" % (self.normalized_impact_angle_list[i], self.normalized_rebound_tangential_surface_vel_list[i]) + '\n')
            self.normalized_rebound_angular_velocity_list_outfile.write("%14.8f %14.8f" % (self.normalized_impact_angle_list[i], self.normalized_rebound_angular_velocity_list[i]) + '\n')
            self.tangential_coefficient_of_restitution_list_outfile.write("%14.8f %14.8f" % (self.normalized_impact_angle_list[i], self.tangential_coefficient_of_restitution_list[i]) + '\n')
        self.normalized_rebound_tangential_surface_vel_list_outfile.close()
        self.normalized_rebound_angular_velocity_list_outfile.close()
        self.tangential_coefficient_of_restitution_list_outfile.close()

        self.create_gnuplot_scripts(self.normalized_rebound_tangential_surface_vel_list_outfile_name,
                                    self.normalized_rebound_angular_velocity_list_outfile_name,
                                    self.tangential_coefficient_of_restitution_list_outfile_name,
                                    self.coeff_of_rest_string, dt)

        error1, error2, error3 = self.compute_errors(self.normalized_rebound_tangential_surface_vel_list_outfile_name,
                                                     self.normalized_rebound_angular_velocity_list_outfile_name,
                                                     self.tangential_coefficient_of_restitution_list_outfile_name)

        coeff_of_rest = '%.2f' % self.coeff_of_restitution

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')

        if (coeff_of_rest=='0.25'):
            error_file.write("\n===== THORNTON PAPER TESTS. FULL REGIME. LINEAR LAW =====\n\n")

        error_file.write("DEM Benchmark 10:")

        if (error1 < 5.0 and error2 < 5.0 and error3 < 5.0):
            error_file.write(" OK!........ Test 10 (e=" + coeff_of_rest + ") SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 10 (e=" + coeff_of_rest + ") FAILED\n")
        error_file.close()

        self.normalized_impact_angle_list = []
        self.normalized_rebound_tangential_surface_vel_list = []
        self.normalized_rebound_angular_velocity_list = []
        self.tangential_coefficient_of_restitution_list = []

    def create_gnuplot_scripts(self, normalized_rebound_tangential_surface_vel_list_outfile_name,
                                     normalized_rebound_angular_velocity_list_outfile_name,
                                     tangential_coefficient_of_restitution_list_outfile_name,
                                     coeff_of_rest_string, dt):

        gnuplot_script_name_1 = 'benchmark10_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Normalized rebound tangential surface velocity'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:10][-2:3] '" + normalized_rebound_tangential_surface_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/bench_10_norm_reb_tang_e_" + coeff_of_rest_string + ".dat' index 1 w lp ls 1 t 'Paper data'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark10_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Normalized final angular velocity'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:14][-6:0] '" + normalized_rebound_angular_velocity_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/bench_10_norm_reb_ang_vel_e_" + coeff_of_rest_string + ".dat' index 1 w lp ls 1 t 'Paper data'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_3 = 'benchmark10_comparison_3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_3, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Tangential coefficient of restitution'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:10][0.5:1.0] '" + tangential_coefficient_of_restitution_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/bench_10_tang_coeff_rest_e_" + coeff_of_rest_string + ".dat' index 1 w lp ls 1 t 'Paper data'\n")
        self.gnuplot_outfile.close()

        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        print_gnuplot_files_on_screen(gnuplot_script_name_3)

        '''

    def compute_errors(self, normalized_rebound_tangential_surface_vel_list_outfile_name,
                             normalized_rebound_angular_velocity_list_outfile_name,
                             tangential_coefficient_of_restitution_list_outfile_name):
        #
        Thornton_data = []; DEM_data = []; summation_of_Thornton_data = 0
        i = 0
        path = "paper_data/bench_10_norm_reb_tang_e_" + self.coeff_of_rest_string + ".dat"
        with open(path) as inf:
            for line in inf:
                if i in self.lines_Thornton:
                    parts = line.split(',')
                    Thornton_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(normalized_rebound_tangential_surface_vel_list_outfile_name) as inf:
            for line in inf:
                if i in self.lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_normalized_rebound_tangential_surface_vel_error = 0

        for j in Thornton_data:
            summation_of_Thornton_data+=abs(j)

        for i, j in zip(DEM_data, Thornton_data):
            final_normalized_rebound_tangential_surface_vel_error+=fabs(i-j)

        final_normalized_rebound_tangential_surface_vel_error/=summation_of_Thornton_data

        print("Error in normalized rebound tangential surface velocity =", 100*final_normalized_rebound_tangential_surface_vel_error,"%")

        #
        Thornton_data = []; DEM_data = []; summation_of_Thornton_data = 0
        i = 0
        path = "paper_data/bench_10_norm_reb_ang_vel_e_" + self.coeff_of_rest_string + ".dat"
        with open(path) as inf:
            for line in inf:
                if i in self.lines_Thornton:
                    parts = line.split(',')
                    Thornton_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(normalized_rebound_angular_velocity_list_outfile_name) as inf:
            for line in inf:
                if i in self.lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_normalized_rebound_angular_velocity_error = 0

        for j in Thornton_data:
            summation_of_Thornton_data+=abs(j)

        for i, j in zip(DEM_data, Thornton_data):
            final_normalized_rebound_angular_velocity_error+=fabs(i-j)

        final_normalized_rebound_angular_velocity_error/=summation_of_Thornton_data
        print("Error in normalized rebound angular velocity =", 100*final_normalized_rebound_angular_velocity_error,"%")

        #
        Thornton_data = []; DEM_data = []; summation_of_Thornton_data = 0
        i = 0
        path = "paper_data/bench_10_tang_coeff_rest_e_" + self.coeff_of_rest_string + ".dat"
        with open(path) as inf:
            for line in inf:
                if i in self.lines_Thornton:
                    parts = line.split(',')
                    Thornton_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(tangential_coefficient_of_restitution_list_outfile_name) as inf:
            for line in inf:
                if i in self.lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_tangential_coefficient_of_restitution_error = 0

        for j in Thornton_data:
            summation_of_Thornton_data+=abs(j)

        for i, j in zip(DEM_data, Thornton_data):
            final_tangential_coefficient_of_restitution_error+=fabs(i-j)

        final_tangential_coefficient_of_restitution_error/=summation_of_Thornton_data
        print("Error in final tangential coefficient of restitution =", 100*final_tangential_coefficient_of_restitution_error,"%")
        #
        error1 = 100*final_normalized_rebound_tangential_surface_vel_error
        error2 = 100*final_normalized_rebound_angular_velocity_error
        error3 = 100*final_tangential_coefficient_of_restitution_error

        return error1, error2, error3

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass


class Benchmark11: ########## HERTZIAN THORNTON

    def __init__(self):
        self.number = 11
        self.initial_normal_vel = -5.0
        self.initial_tangential_vel = 0
        self.radius = 0.025
        self.normalized_impact_angle_list = []
        self.normalized_rebound_tangential_surface_vel_list = []
        self.normalized_rebound_angular_velocity_list = []
        self.tangential_coefficient_of_restitution_list = []
        self.normalized_rebound_tangential_surface_vel_list_outfile = None
        self.normalized_rebound_angular_velocity_list_outfile = None
        self.tangential_coefficient_of_restitution_list_outfile = None
        self.coeff_of_restitution = -1.0
        self.coeff_of_rest_string = None
        self.lines_Thornton = []
        self.lines_DEM = []
        self.degrees = 0

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):   # Change this function name from 'set_initial_data' to 'set_initial_data'

        if iteration == 1:
            self.degrees = 1
        else:
            self.degrees = 50 * (iteration - 1)/number_of_points_in_the_graphic

        if coeff_of_restitution_iteration==1:
            self.coeff_of_restitution=0.25
            self.coeff_of_rest_string='025'
            self.lines_Thornton = [1, 2, 4, 5, 7, 8]
            self.lines_DEM = [0, 1, 3, 4, 5, 6]
        elif coeff_of_restitution_iteration==2:
            self.coeff_of_restitution=0.50
            self.coeff_of_rest_string='050'
            self.lines_Thornton = [1, 2, 4, 5, 7, 9, 10]
            self.lines_DEM = [0, 1, 3, 4, 5, 6, 7]
        elif coeff_of_restitution_iteration==3:
            self.coeff_of_restitution=0.75
            self.coeff_of_rest_string='075'
            self.lines_Thornton = [1, 2, 4, 5, 6, 8, 9, 10]
            self.lines_DEM = [0, 1, 3, 4, 5, 6, 7, 8]
        else:
            self.coeff_of_restitution=0.90
            self.coeff_of_rest_string='090'
            #self.lines_Thornton = [1, 2, 4, 5, 6, 8, 9]
            #self.lines_DEM = [0, 1, 3, 4, 5, 7, 8]
            self.lines_Thornton = [1, 2, 4, 5, 6, 7, 8, 9]
            self.lines_DEM = [0, 1, 3, 4, 5, 6, 7, 8]

        self.initial_tangential_vel = -self.initial_normal_vel * tan(self.degrees * pi / 180.0)

        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.initial_tangential_vel)
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)
            modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = self.coeff_of_restitution

        print(self.coeff_of_restitution)

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):

        mu = 0.1

        for node in modelpart.Nodes:
            final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_X)
            final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Y)
            normalized_impact_angle = 2.0 * tan(self.degrees * pi / 180.0) / (mu * (1 + self.coeff_of_restitution))
            normalized_rebound_tangential_surface_vel = -2.0 * (final_tangential_center_velocity + final_angular_vel * self.radius) / (self.initial_normal_vel * mu * (1 + self.coeff_of_restitution))
            normalized_rebound_angular_velocity = -2.0 * self.radius * final_angular_vel / (self.initial_normal_vel * mu * (1 + self.coeff_of_restitution))
            tangential_coefficient_of_restitution = 5.0/7.0 + 2.0 * normalized_rebound_tangential_surface_vel / (7.0 * normalized_impact_angle)

        self.normalized_impact_angle_list.append(normalized_impact_angle)
        self.normalized_rebound_tangential_surface_vel_list.append(normalized_rebound_tangential_surface_vel)
        self.normalized_rebound_angular_velocity_list.append(normalized_rebound_angular_velocity)
        self.tangential_coefficient_of_restitution_list.append(tangential_coefficient_of_restitution)

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        self.normalized_rebound_tangential_surface_vel_list_outfile_name = "benchmark11_dt_" + str(dt) + '_normalized_rebound_tangential_surface_vel_list_data.dat'
        self.normalized_rebound_angular_velocity_list_outfile_name = "benchmark11_dt_" + str(dt) + '_normalized_rebound_angular_velocity_list_data.dat'
        self.tangential_coefficient_of_restitution_list_outfile_name = "benchmark11_dt_" + str(dt) + '_tangential_coefficient_of_restitution_list_data.dat'

        self.normalized_rebound_tangential_surface_vel_list_outfile = open(self.normalized_rebound_tangential_surface_vel_list_outfile_name, 'w')
        self.normalized_rebound_angular_velocity_list_outfile = open(self.normalized_rebound_angular_velocity_list_outfile_name, 'w')
        self.tangential_coefficient_of_restitution_list_outfile = open(self.tangential_coefficient_of_restitution_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.normalized_rebound_tangential_surface_vel_list_outfile.write("%14.8f %14.8f" % (self.normalized_impact_angle_list[i], self.normalized_rebound_tangential_surface_vel_list[i]) + '\n')
            self.normalized_rebound_angular_velocity_list_outfile.write("%14.8f %14.8f" % (self.normalized_impact_angle_list[i], self.normalized_rebound_angular_velocity_list[i]) + '\n')
            self.tangential_coefficient_of_restitution_list_outfile.write("%14.8f %14.8f" % (self.normalized_impact_angle_list[i], self.tangential_coefficient_of_restitution_list[i]) + '\n')
        self.normalized_rebound_tangential_surface_vel_list_outfile.close()
        self.normalized_rebound_angular_velocity_list_outfile.close()
        self.tangential_coefficient_of_restitution_list_outfile.close()

        self.create_gnuplot_scripts(self.normalized_rebound_tangential_surface_vel_list_outfile_name,
                                    self.normalized_rebound_angular_velocity_list_outfile_name,
                                    self.tangential_coefficient_of_restitution_list_outfile_name,
                                    self.coeff_of_rest_string, dt)

        error1, error2, error3 = self.compute_errors(self.normalized_rebound_tangential_surface_vel_list_outfile_name,
                                                     self.normalized_rebound_angular_velocity_list_outfile_name,
                                                     self.tangential_coefficient_of_restitution_list_outfile_name)

        coeff_of_rest = '%.2f' % self.coeff_of_restitution

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')

        if (coeff_of_rest=='0.25'):
            error_file.write("\n==== THORNTON PAPER TESTS. FULL REGIME. HERTZIAN LAW ====\n\n")

        error_file.write("DEM Benchmark 11:")

        if (error1 < 6.0 and error2 < 6.0 and error3 < 6.0):
            error_file.write(" OK!........ Test 11 (e=" + coeff_of_rest + ") SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 11 (e=" + coeff_of_rest + ") FAILED\n")
        error_file.close()

        self.normalized_impact_angle_list = []
        self.normalized_rebound_tangential_surface_vel_list = []
        self.normalized_rebound_angular_velocity_list = []
        self.tangential_coefficient_of_restitution_list = []

    def create_gnuplot_scripts(self, normalized_rebound_tangential_surface_vel_list_outfile_name,
                                     normalized_rebound_angular_velocity_list_outfile_name,
                                     tangential_coefficient_of_restitution_list_outfile_name,
                                     coeff_of_rest_string, dt):

        gnuplot_script_name_1 = 'benchmark11_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Normalized rebound tangential surface velocity'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:10][-2:3] '" + normalized_rebound_tangential_surface_vel_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/bench_10_norm_reb_tang_e_" + coeff_of_rest_string + ".dat' index 0 w lp ls 1 t 'Paper data'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_2 = 'benchmark11_comparison_2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_2, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Normalized final angular velocity'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:14][-6:0] '" + normalized_rebound_angular_velocity_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/bench_10_norm_reb_ang_vel_e_" + coeff_of_rest_string + ".dat' index 0 w lp ls 1 t 'Paper data'\n")
        self.gnuplot_outfile.close()

        gnuplot_script_name_3 = 'benchmark11_comparison_3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_3, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Normalized incident angle'\nset ylabel 'Tangential coefficient of restitution'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:10][0.5:1.0] '" + tangential_coefficient_of_restitution_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/bench_10_tang_coeff_rest_e_" + coeff_of_rest_string + ".dat' index 0 w lp ls 1 t 'Paper data'\n")
        self.gnuplot_outfile.close()

        '''
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        print_gnuplot_files_on_screen(gnuplot_script_name_3)
        '''

    def compute_errors(self, normalized_rebound_tangential_surface_vel_list_outfile_name,
                             normalized_rebound_angular_velocity_list_outfile_name,
                             tangential_coefficient_of_restitution_list_outfile_name):
        #
        Thornton_data = []; DEM_data = []; summation_of_Thornton_data = 0
        i = 0
        path = "paper_data/bench_10_norm_reb_tang_e_" + self.coeff_of_rest_string + ".dat"
        with open(path) as inf:
            for line in inf:
                if i in self.lines_Thornton:
                    parts = line.split(',')
                    Thornton_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(normalized_rebound_tangential_surface_vel_list_outfile_name) as inf:
            for line in inf:
                if i in self.lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_normalized_rebound_tangential_surface_vel_error = 0

        for j in Thornton_data:
            summation_of_Thornton_data+=abs(j)

        for i, j in zip(DEM_data, Thornton_data):
            final_normalized_rebound_tangential_surface_vel_error+=fabs(i-j)

        final_normalized_rebound_tangential_surface_vel_error/=summation_of_Thornton_data

        print("Error in normalized rebound tangential surface velocity =", 100*final_normalized_rebound_tangential_surface_vel_error,"%")

        #
        Thornton_data = []; DEM_data = []; summation_of_Thornton_data = 0
        i = 0
        path = "paper_data/bench_10_norm_reb_ang_vel_e_" + self.coeff_of_rest_string + ".dat"
        with open(path) as inf:
            for line in inf:
                if i in self.lines_Thornton:
                    parts = line.split(',')
                    Thornton_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(normalized_rebound_angular_velocity_list_outfile_name) as inf:
            for line in inf:
                if i in self.lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_normalized_rebound_angular_velocity_error = 0

        for j in Thornton_data:
            summation_of_Thornton_data+=abs(j)

        for i, j in zip(DEM_data, Thornton_data):
            final_normalized_rebound_angular_velocity_error+=fabs(i-j)

        final_normalized_rebound_angular_velocity_error/=summation_of_Thornton_data
        print("Error in normalized rebound angular velocity =", 100*final_normalized_rebound_angular_velocity_error,"%")

        #
        Thornton_data = []; DEM_data = []; summation_of_Thornton_data = 0
        i = 0
        path = "paper_data/bench_10_tang_coeff_rest_e_" + self.coeff_of_rest_string + ".dat"
        with open(path) as inf:
            for line in inf:
                if i in self.lines_Thornton:
                    parts = line.split(',')
                    Thornton_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(tangential_coefficient_of_restitution_list_outfile_name) as inf:
            for line in inf:
                if i in self.lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_tangential_coefficient_of_restitution_error = 0

        for j in Thornton_data:
            summation_of_Thornton_data+=abs(j)

        for i, j in zip(DEM_data, Thornton_data):
            final_tangential_coefficient_of_restitution_error+=fabs(i-j)

        final_tangential_coefficient_of_restitution_error/=summation_of_Thornton_data
        print("Error in final tangential coefficient of restitution =", 100*final_tangential_coefficient_of_restitution_error,"%")

        #
        error1 = 100*final_normalized_rebound_tangential_surface_vel_error
        error2 = 100*final_normalized_rebound_angular_velocity_error
        error3 = 100*final_tangential_coefficient_of_restitution_error

        return error1, error2, error3

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        pass


class Benchmark12: ########## ROLLING FRICTION

    def __init__(self):
        self.number = 12

        self.balls_graph_counter = 1

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            total_angular_velocity_z = 0.0

            for node in modelpart.Nodes:
                if node.Id == 1:
                   angular_velocity_z = node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)
                   total_angular_velocity_z += angular_velocity_z

                del node

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_angular_velocity_z).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("==== WENSRICH PAPER TEST. ROLLING FRICTION ====\n\n")
        error_file.write("DEM Benchmark 12:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 12 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 12 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):

        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")

        error1 = 100*generated_data_error

        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass


class Benchmark13: ########## DEM-FEM Facet

    def __init__(self):
        self.number = 13

        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.velocity_list_outfile_name, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            total_velocity_x = 0.0
            total_velocity_z = 0.0

            for node in modelpart.Nodes:
                if node.Id == 1:
                   velocity_x = node.GetSolutionStepValue(VELOCITY_X)
                   velocity_z = node.GetSolutionStepValue(VELOCITY_Z)
                   total_velocity_x += velocity_x
                   total_velocity_z += velocity_z

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_x).rjust(13)+" "+str("%.6g"%total_velocity_z).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("======== DE/FE CONTACT BENCHMARKS ==========\n\n")
        error_file.write("DEM Benchmark 13:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 13 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 13 FAILED\n")
        error_file.close()

    def compute_errors(self, velocity_list_outfile_name):  #FINALIZATION STEP

        lines_DEM = list(range(0, 200));
        total_velocity_x = 0.0; total_velocity_z = 0.0
        i = 0
        with open(velocity_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    total_velocity_x += float(parts[1])
                    total_velocity_z += float(parts[2])
                i+=1

        if total_velocity_x > 0.0:  #VELOCITY_X should be 0 always
            error1 = 100
        else:
            error1 = 0

        if total_velocity_z > 0.0:  #VELOCITY_Z should be 0 always
            error2 = 100
        else:
            error2 = 0

        error3 = 0

        print("Error in velocity X =", error1,"%")

        print("Error in velocity Z =", error2,"%")

        return error1, error2, error3

class Benchmark14: ########## DEM-FEM Edge

    def __init__(self):
        self.number = 14

        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.velocity_list_outfile_name, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            total_velocity_x = 0.0
            total_velocity_z = 0.0

            for node in modelpart.Nodes:
                if node.Id == 1:
                   velocity_x = node.GetSolutionStepValue(VELOCITY_X)
                   velocity_z = node.GetSolutionStepValue(VELOCITY_Z)
                   total_velocity_x += velocity_x
                   total_velocity_z += velocity_z

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_x).rjust(13)+" "+str("%.6g"%total_velocity_z).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 14:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 14 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 14 FAILED\n")
        error_file.close()

    def compute_errors(self, velocity_list_outfile_name):  #FINALIZATION STEP

        lines_DEM = list(range(0, 200));
        total_velocity_x = 0.0; total_velocity_z = 0.0
        i = 0
        with open(velocity_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    total_velocity_x += float(parts[1])
                    total_velocity_z += float(parts[2])
                i+=1

        if total_velocity_x > 0.0:  #VELOCITY_X should be 0 always
            error1 = 100
        else:
            error1 = 0

        if total_velocity_z > 0.0:  #VELOCITY_Z should be 0 always
            error2 = 100
        else:
            error2 = 0

        error3 = 0

        print("Error in velocity X =", error1,"%")

        print("Error in velocity Z =", error2,"%")

        return error1, error2, error3

class Benchmark15: ########## DEM-FEM Vertex

    def __init__(self):
        self.number = 15

        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.velocity_list_outfile_name, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            total_velocity_x = 0.0
            total_velocity_z = 0.0

            for node in modelpart.Nodes:
                if node.Id == 1:
                   velocity_x = node.GetSolutionStepValue(VELOCITY_X)
                   velocity_z = node.GetSolutionStepValue(VELOCITY_Z)
                   total_velocity_x += velocity_x
                   total_velocity_z += velocity_z

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_x).rjust(13)+" "+str("%.6g"%total_velocity_z).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 15:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 15 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 15 FAILED\n")
        error_file.close()

    def compute_errors(self, velocity_list_outfile_name):  #FINALIZATION STEP

        lines_DEM = list(range(0, 200));
        total_velocity_x = 0.0; total_velocity_z = 0.0
        i = 0
        with open(velocity_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    total_velocity_x += float(parts[1])
                    total_velocity_z += float(parts[2])
                i+=1

        if total_velocity_x > 0.0:  #VELOCITY_X should be 0 always
            error1 = 100
        else:
            error1 = 0

        if total_velocity_z > 0.0:  #VELOCITY_Z should be 0 always
            error2 = 100
        else:
            error2 = 0

        error3 = 0

        print("Error in velocity X =", error1,"%")

        print("Error in velocity Z =", error2,"%")

        return error1, error2, error3


class Benchmark16: ########## DEM-FEM Grid

    def __init__(self):
        self.number = 16

        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.velocity_list_outfile_name, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            total_velocity_1 = 0.0
            total_velocity_2 = 0.0
            total_velocity_3 = 0.0

            for node in modelpart.Nodes:
                if node.Id == 1:
                   velocity_1 = node.GetSolutionStepValue(VELOCITY_Y)
                   total_velocity_1 += velocity_1
                if node.Id == 2:
                   velocity_2 = node.GetSolutionStepValue(VELOCITY_Y)
                   total_velocity_2 += velocity_2
                if node.Id == 3:
                   velocity_3 = node.GetSolutionStepValue(VELOCITY_Y)
                   total_velocity_3 += velocity_3

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_1).rjust(13)+" "+str("%.6g"%total_velocity_2).rjust(13)+" "+str("%.6g"%total_velocity_3).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 16:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 16 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 16 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):  #FINALIZATION STEP

        lines_analytics = lines_DEM = list(range(0, 250));
        ref_data1 = []; ref_data2 = []; DEM_data1 = []; ref_data3 = []; DEM_data1 = []; DEM_data2 = []; DEM_data3 = []; summation_of_ref_data1 = 0; summation_of_ref_data2 = 0; summation_of_ref_data3 = 0
        times = []
        i = 0
        with open('paper_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:  #with open('paper_data/reference_graph_benchmark12.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    times.append(float(parts[0]))
                    ref_data1.append(float(parts[1]))
                    ref_data2.append(float(parts[2]))
                    ref_data3.append(float(parts[3]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data1.append(float(parts[1]))
                    DEM_data2.append(float(parts[2]))
                    DEM_data3.append(float(parts[3]))
                i+=1
        final_velocity_1_error = 0
        final_velocity_2_error = 0
        final_velocity_3_error = 0

        for j in ref_data1:
            summation_of_ref_data1+=abs(j)
        for k in ref_data2:
            summation_of_ref_data2+=abs(k)
        for l in ref_data3:
            summation_of_ref_data3+=abs(l)

        for i, j in zip(DEM_data1, ref_data1):
            final_velocity_1_error+=fabs(i-j)
        final_velocity_1_error/=summation_of_ref_data1

        for k, l in zip(DEM_data2, ref_data2):
            final_velocity_2_error+=fabs(k-l)
        final_velocity_2_error/=summation_of_ref_data2

        for m, n in zip(DEM_data3, ref_data3):
            final_velocity_3_error+=fabs(m-n)
        final_velocity_3_error/=summation_of_ref_data3

        #for t, v1,v2,v3 in zip(times, DEM_data1, DEM_data2, DEM_data3):
        #    print(t, v1, v2, v3)

        print("Error in velocity sphere 1 =", 100*final_velocity_1_error,"%")

        print("Error in velocity sphere 2 =", 100*final_velocity_2_error,"%")

        print("Error in velocity sphere 3 =", 100*final_velocity_3_error,"%")

        error1 = 100*final_velocity_1_error

        error2 = 100*final_velocity_2_error

        error3 = 100*final_velocity_3_error

        return error1, error2, error3


class Benchmark17: ########## DEM-FEM Rolling

    def __init__(self):
        self.number = 17

        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.error_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.error_list_outfile_name, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            total_velocity_err         = 0.0
            total_angular_velocity_err = 0.0

            for node in modelpart.Nodes:
                if node.Id == 1:
                   velocity_1         = node.GetSolutionStepValue(VELOCITY_X)
                   angular_velocity_1 = node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)
                if node.Id == 2:
                   velocity_2         = node.GetSolutionStepValue(VELOCITY_X)
                   angular_velocity_2 = node.GetSolutionStepValue(ANGULAR_VELOCITY_Z)

            total_velocity_err         = (abs(velocity_1 - velocity_2))/(abs(velocity_2))
            total_angular_velocity_err = (abs(angular_velocity_1 - angular_velocity_2))/(abs(velocity_2))

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_err).rjust(13)+" "+str("%.6g"%total_angular_velocity_err).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.error_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 17:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 17 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 17 FAILED\n")
        error_file.close()

    def compute_errors(self, error_list_outfile_name):  #FINALIZATION STEP

        lines_DEM = list(range(0, 100));
        total_velocity_err = 0.0; total_angular_velocity_err = 0.0
        i = 0
        with open(error_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    total_velocity_err += float(parts[1])
                    total_angular_velocity_err += float(parts[2])
                i+=1

        if total_velocity_err > 1e-2:  #VELOCITY_X should be 0 always
            error1 = 100*total_velocity_err
        else:
            error1 = 0

        if total_angular_velocity_err > 1e-2:  #VELOCITY_Z should be 0 always
            error2 = 100*total_angular_velocity_err
        else:
            error2 = 0

        error3 = 0

        print("Error in velocity between meshes =", 100*total_velocity_err,"%")

        print("Error in angular velocity between meshes =", 100*total_angular_velocity_err,"%")

        return error1, error2, error3


class Benchmark20:

    def __init__(self):
        self.number = 20
        self.generated_data = None
        #self.graph_frequency = int(graph_print_interval/dt)  # def __init__(self, graph_print_interval, dt):
        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):  #INITIALIZATION STEP
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP
        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP
        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):

        '''
        gnuplot_script_name = 'benchmark3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + self.output_filename + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4, '"\
                                                      + self.output_filename + "' u 1:3 w lp lt 2 lw 1.5 ps 2 pt 6")
        self.gnuplot_outfile.close()
        self.create_gnuplot_scripts(self.output_filename, dt)
        '''

        error1, error2, error3 = self.compute_errors(self.output_filename)
        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("== BASIC CONTINUUM TESTS ==\n\n")
        error_file.write("DEM Benchmark 20:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 20 SUCCESSFUL\n")
            shutil.rmtree('benchmark20_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 20 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):  #FINALIZATION STEP
        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + str(sys.argv[1]) + '.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))   #segona component del vector ()
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")

        error1 = 100*generated_data_error

        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass
        '''
        gnuplot_script_name_1 = 'benchmark20_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Data'\nset ylabel 'Damping ratio'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:1][0:1] '" + output_filename + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark20_graph1.dat' w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark20_graph1.dat' w lp ls 2 t 'Cast iron'\n")
        self.gnuplot_outfile.close()
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        '''

class Benchmark21:

    def __init__(self):
        self.number = 21
        self.generated_data = None
        #self.graph_frequency = int(graph_print_interval/dt)  # def __init__(self, graph_print_interval, dt):
        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):  #INITIALIZATION STEP
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):        #FINALIZATION STEP
        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP
        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):
        error1, error2, error3 = self.compute_errors(self.output_filename)
        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 21:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 21 SUCCESSFUL\n")
            shutil.rmtree('benchmark21_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 21 FAILED\n")
        error_file.close()


    def compute_errors(self, output_filename):  #FINALIZATION STEP
        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + str(sys.argv[1]) + '.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))   #segona component del vector ()
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")
        error1 = 100*generated_data_error
        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass


class Benchmark22:

    def __init__(self):
        self.number = 22
        self.generated_data = None
        #self.graph_frequency = int(graph_print_interval/dt)  # def __init__(self, graph_print_interval, dt):
        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):  #INITIALIZATION STEP
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):        #FINALIZATION STEP
        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP
        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):
        error1, error2, error3 = self.compute_errors(self.output_filename)
        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 22:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 22 SUCCESSFUL\n")
            shutil.rmtree('benchmark22_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 22 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + str(sys.argv[1]) + '.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))   #segona component del vector ()
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")
        error1 = 100*generated_data_error
        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass

class Benchmark23:

    def __init__(self):
        self.number = 23
        self.generated_data = None
        #self.graph_frequency = int(graph_print_interval/dt)  # def __init__(self, graph_print_interval, dt):
        self.balls_graph_counter = 1   # deberia ser self.balls_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):  #INITIALIZATION STEP
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):        #FINALIZATION STEP
        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP
        #print("generate_graph_points bench23, graph_print_interval, dt - ", graph_print_interval, dt )
        self.graph_frequency        = int(graph_print_interval/dt)

        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.balls_graph_counter == self.graph_frequency):     #if(self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 23:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 23 SUCCESSFUL\n")
            shutil.rmtree('benchmark23_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 23 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):  #FINALIZATION STEP
        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + '23' + '.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))   #segona component del vector ()
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")
        error1 = 100*generated_data_error
        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass


class Benchmark24:

    def __init__(self):
        self.number = 24
        self.generated_data = None
        self.balls_graph_counter = 1

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        self.simulation_graph.close()

    def cross_product(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
        return c

    def ApplyNodalRotation(self, time, dt, modelpart):
        ang_vel = 20 * pi
        angular_velocity = [0, 0, ang_vel]
        rotation_matrix = [[cos(ang_vel * time), -1.0 * sin(ang_vel * time), 0], [sin(ang_vel * time), cos(ang_vel * time), 0], [0,0,1]]
        time_dt = time - dt
        rotation_matrix_minus_dt = [[cos(ang_vel * time_dt), -1.0 * sin(ang_vel * time_dt), 0], [sin(ang_vel * time_dt), cos(ang_vel * time_dt), 0], [0,0,1]] #
        centroid = [-1.0, 0.0, 0.0]
        relative_initial_node_coords, relative_node_coords, relative_node_coords_dt = [0]*3, [0]*3, [0]*3
        sum, sum_dt = 0, 0

        for node in modelpart.Nodes:
            if node.Id == 141:
                for j in range(3):
                    rot_mat = rotation_matrix[j]
                    rot_mat_dt = rotation_matrix_minus_dt[j]
                    relative_initial_node_coords[0] = node.X0 - centroid[0]
                    relative_initial_node_coords[1] = node.Y0 - centroid[1]
                    relative_initial_node_coords[2] = node.Z0 - centroid[2]
                    for i in range(3):
                        sum += rot_mat[i] * relative_initial_node_coords[i]
                        sum_dt += rot_mat_dt[i] * relative_initial_node_coords[i]
                    relative_node_coords[j], sum, relative_node_coords_dt[j], sum_dt = sum, 0, sum_dt, 0
                node.X = relative_node_coords[0] + centroid[0]
                node.Y = relative_node_coords[1] + centroid[1]
                node.Z = relative_node_coords[2] + centroid[2]

                displacement = [0]*3
                displacement[0] = node.X - node.X0
                displacement[1] = node.Y - node.Y0
                displacement[2] = node.Z - node.Z0
                node.SetSolutionStepValue(DISPLACEMENT, displacement)

                velocity = [0]*3
                velocity = self.cross_product(angular_velocity, relative_node_coords)
                node.SetSolutionStepValue(VELOCITY, velocity)

                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

                delta_displacement = [0]*3
                delta_displacement[0] = relative_node_coords[0] - relative_node_coords_dt[0]
                delta_displacement[1] = relative_node_coords[1] - relative_node_coords_dt[1]
                delta_displacement[2] = relative_node_coords[2] - relative_node_coords_dt[2]
                node.SetSolutionStepValue(DELTA_DISPLACEMENT, delta_displacement)

                particle_rotation_angle = [0]*3
                particle_rotation_angle[0] = angular_velocity[0] * time
                particle_rotation_angle[1] = angular_velocity[1] * time
                particle_rotation_angle[2] = angular_velocity[2] * time
                node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE, particle_rotation_angle)

                delta_rotation = [0]*3
                delta_rotation[0] = angular_velocity[0] * dt
                delta_rotation[1] = angular_velocity[1] * dt
                delta_rotation[2] = angular_velocity[2] * dt
                node.SetSolutionStepValue(DELTA_ROTATION, delta_rotation)

            if node.Id == 140:
                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        #print("generate_graph_points bench24, graph_print_interval, dt - ", graph_print_interval, dt )
        self.graph_frequency = int(graph_print_interval/dt)

        if self.graph_frequency < 1:
           self.graph_frequency = 1

        if (self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 24:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 24 SUCCESSFUL\n")
            shutil.rmtree('benchmark24_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 24 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + '24' + '.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[2]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[2]))   #segona component del vector ()
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")
        error1 = 100*generated_data_error
        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass


class Benchmark25:

    def __init__(self):
        self.number = 25
        self.generated_data = None
        self.balls_graph_counter = 1

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        self.simulation_graph.close()

    def cross_product(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
        return c

    def ApplyNodalRotation(self, time, dt, modelpart):
        ang_vel = 20 * pi
        angular_velocity = [0, 0, ang_vel]
        rotation_matrix = [[cos(ang_vel * time), -1.0 * sin(ang_vel * time), 0], [sin(ang_vel * time), cos(ang_vel * time), 0], [0,0,1]]
        time_dt = time - dt
        rotation_matrix_minus_dt = [[cos(ang_vel * time_dt), -1.0 * sin(ang_vel * time_dt), 0], [sin(ang_vel * time_dt), cos(ang_vel * time_dt), 0], [0,0,1]] #
        centroid = [-1.0, 0.0, 0.0]
        relative_initial_node_coords, relative_node_coords, relative_node_coords_dt = [0]*3, [0]*3, [0]*3
        sum, sum_dt = 0, 0

        for node in modelpart.Nodes:
            if node.Id == 141:
                for j in range(3):
                    rot_mat = rotation_matrix[j]
                    rot_mat_dt = rotation_matrix_minus_dt[j]
                    relative_initial_node_coords[0] = node.X0 - centroid[0]
                    relative_initial_node_coords[1] = node.Y0 - centroid[1]
                    relative_initial_node_coords[2] = node.Z0 - centroid[2]
                    for i in range(3):
                        sum += rot_mat[i] * relative_initial_node_coords[i]
                        sum_dt += rot_mat_dt[i] * relative_initial_node_coords[i]
                    relative_node_coords[j], sum, relative_node_coords_dt[j], sum_dt = sum, 0, sum_dt, 0
                node.X = relative_node_coords[0] + centroid[0]
                node.Y = relative_node_coords[1] + centroid[1]
                node.Z = relative_node_coords[2] + centroid[2]

                displacement = [0]*3
                displacement[0] = node.X - node.X0
                displacement[1] = node.Y - node.Y0
                displacement[2] = node.Z - node.Z0
                node.SetSolutionStepValue(DISPLACEMENT, displacement)

                velocity = [0]*3
                velocity = self.cross_product(angular_velocity, relative_node_coords)
                node.SetSolutionStepValue(VELOCITY, velocity)

                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

                delta_displacement = [0]*3
                delta_displacement[0] = relative_node_coords[0] - relative_node_coords_dt[0]
                delta_displacement[1] = relative_node_coords[1] - relative_node_coords_dt[1]
                delta_displacement[2] = relative_node_coords[2] - relative_node_coords_dt[2]
                node.SetSolutionStepValue(DELTA_DISPLACEMENT, delta_displacement)

                particle_rotation_angle = [0]*3
                particle_rotation_angle[0] = angular_velocity[0] * time
                particle_rotation_angle[1] = angular_velocity[1] * time
                particle_rotation_angle[2] = angular_velocity[2] * time
                node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE, particle_rotation_angle)

                delta_rotation = [0]*3
                delta_rotation[0] = angular_velocity[0] * dt
                delta_rotation[1] = angular_velocity[1] * dt
                delta_rotation[2] = angular_velocity[2] * dt
                node.SetSolutionStepValue(DELTA_ROTATION, delta_rotation)

                if time > 3.8e-5:
                    radius = 1.0001
                    node.SetSolutionStepValue(RADIUS, radius)
            if node.Id == 140:
                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):
        self.graph_frequency = int(graph_print_interval/dt)

        if self.graph_frequency < 1:
           self.graph_frequency = 1

        if (self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):
        error1, error2, error3 = self.compute_errors(self.output_filename)
        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 25:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 25 SUCCESSFUL\n")
            shutil.rmtree('benchmark25_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 25 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        lines_analytics = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + '25' + '.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    analytics_data.append(float(parts[2]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[2]))   #segona component del vector ()
                i+=1
        generated_data_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            generated_data_error+=fabs(i-j)
        generated_data_error/=summation_of_analytics_data

        print("Error in simulation =", 100*generated_data_error,"%")
        error1 = 100*generated_data_error
        error2 = error3 = 0

        return error1, error2, error3

    def create_gnuplot_scripts(self, output_filename, dt):
        pass


class Benchmark26:

    def __init__(self):
        self.number = 26

        self.generated_data = None
        self.balls_graph_counter = 1

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):

        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):

        self.graph_frequency = int(graph_print_interval/dt)

        if self.graph_frequency < 1:
           self.graph_frequency = 1

        if (self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0

            for node in modelpart.Nodes:
                if node.Id == 141:
                    force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                    force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                    self.total_force_x += force_node_x
                    self.total_force_y += force_node_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_x).rjust(13)+" "+str("%.6g"%self.total_force_y).rjust(13)+"\n")
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt = 0):
        pass

    def compute_errors(self, output_filename):
        pass

    def create_gnuplot_scripts(self, output_filename, dt):
        pass


class Benchmark27:

    def __init__(self):
        self.number = 27
        self.generated_data = None
        self.balls_graph_counter = 1
        self.rigid_graph_counter = 1

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):

        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.rigid_face_file = "benchmark" + str(sys.argv[1]) + '_rigid_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')
        self.rigid_graph = open(self.rigid_face_file, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        self.simulation_graph.close()
        self.rigid_graph.close()

    def cross_product(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
        return c

    def ApplyNodalRotation(self, time, dt, modelpart):

        ang_vel = 20 * pi
        angular_velocity = [0, 0, ang_vel]
        rotation_matrix = [[cos(ang_vel * time), -1.0 * sin(ang_vel * time), 0], [sin(ang_vel * time), cos(ang_vel * time), 0], [0,0,1]]
        time_dt = time - dt
        rotation_matrix_minus_dt = [[cos(ang_vel * time_dt), -1.0 * sin(ang_vel * time_dt), 0], [sin(ang_vel * time_dt), cos(ang_vel * time_dt), 0], [0,0,1]] #
        centroid = [-1.0, 0.0, 0.0]
        relative_initial_node_coords, relative_node_coords, relative_node_coords_dt = [0]*3, [0]*3, [0]*3
        sum, sum_dt = 0, 0

        for node in modelpart.Nodes:
            if node.Id == 999999:
                for j in range(3):
                    rot_mat = rotation_matrix[j]
                    rot_mat_dt = rotation_matrix_minus_dt[j]
                    relative_initial_node_coords[0] = node.X0 - centroid[0]
                    relative_initial_node_coords[1] = node.Y0 - centroid[1]
                    relative_initial_node_coords[2] = node.Z0 - centroid[2]
                    for i in range(3):
                        sum += rot_mat[i] * relative_initial_node_coords[i]
                        sum_dt += rot_mat_dt[i] * relative_initial_node_coords[i]
                    relative_node_coords[j], sum, relative_node_coords_dt[j], sum_dt = sum, 0, sum_dt, 0
                node.X = relative_node_coords[0] + centroid[0]
                node.Y = relative_node_coords[1] + centroid[1]
                node.Z = relative_node_coords[2] + centroid[2]

                displacement = [0]*3
                displacement[0] = node.X - node.X0
                displacement[1] = node.Y - node.Y0
                displacement[2] = node.Z - node.Z0
                node.SetSolutionStepValue(DISPLACEMENT, displacement)

                velocity = [0]*3
                velocity = self.cross_product(angular_velocity, relative_node_coords)
                node.SetSolutionStepValue(VELOCITY, velocity)

                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

                delta_displacement = [0]*3
                delta_displacement[0] = relative_node_coords[0] - relative_node_coords_dt[0]
                delta_displacement[1] = relative_node_coords[1] - relative_node_coords_dt[1]
                delta_displacement[2] = relative_node_coords[2] - relative_node_coords_dt[2]
                node.SetSolutionStepValue(DELTA_DISPLACEMENT, delta_displacement)

                particle_rotation_angle = [0]*3
                particle_rotation_angle[0] = angular_velocity[0] * time
                particle_rotation_angle[1] = angular_velocity[1] * time
                particle_rotation_angle[2] = angular_velocity[2] * time
                node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE, particle_rotation_angle)

                delta_rotation = [0]*3
                delta_rotation[0] = angular_velocity[0] * dt
                delta_rotation[1] = angular_velocity[1] * dt
                delta_rotation[2] = angular_velocity[2] * dt
                node.SetSolutionStepValue(DELTA_ROTATION, delta_rotation)

                if time > 3.8e-5:
                    radius = 1.0001
                    node.SetSolutionStepValue(RADIUS, radius)
            if node.Id == 99999:
                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):

        #self.graph_frequency = int(5e-7/dt)   #graph_print_interval/dt
        self.graph_frequency = int(graph_print_interval/1/dt)   #1 veces mas grf que bin
        #print (self.graph_frequency)
        #print (self.balls_graph_counter)
        if self.graph_frequency < 1:
           self.graph_frequency = 1

        if (self.balls_graph_counter == self.graph_frequency):
            #print (self.balls_graph_counter)
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0
            self.total_force_z = 0.0
            self.total_force_sum = 0.0

            self.total_angular_x = 0.0
            self.total_angular_y = 0.0
            self.total_angular_z = 0.0
            self.total_angular_sum = 0.0

            self.total_delta_x = 0.0
            self.total_delta_y = 0.0
            self.total_delta_z = 0.0
            self.total_delta_sum = 0.0

            for node in modelpart.Nodes:
                if node.Id == 29:
                   force_node_x = node.GetSolutionStepValue(ELASTIC_FORCES)[0]
                   force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                   force_node_z = node.GetSolutionStepValue(ELASTIC_FORCES)[2]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y
                   self.total_force_z += force_node_z

                   angular_node_x = node.GetSolutionStepValue(ANGULAR_VELOCITY)[0]
                   angular_node_y = node.GetSolutionStepValue(ANGULAR_VELOCITY)[1]
                   angular_node_z = node.GetSolutionStepValue(ANGULAR_VELOCITY)[2]
                   self.total_angular_x += angular_node_x
                   self.total_angular_y += angular_node_y
                   self.total_angular_z += angular_node_z

                   delta_node_x = node.GetSolutionStepValue(DELTA_DISPLACEMENT)[0]
                   delta_node_y = node.GetSolutionStepValue(DELTA_DISPLACEMENT)[1]
                   delta_node_z = node.GetSolutionStepValue(DELTA_DISPLACEMENT)[2]
                   self.total_delta_x += delta_node_x
                   self.total_delta_y += delta_node_y
                   self.total_delta_z += delta_node_z

            self.total_force_sum = self.total_force_x + self.total_force_y + self.total_force_z
            self.total_angular_sum = self.total_angular_x + self.total_angular_y + self.total_angular_z
            self.total_delta_sum = self.total_delta_x + self.total_delta_y + self.total_delta_z
            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_sum).rjust(13)+" "+str("%.6g"%self.total_angular_sum).rjust(13)+" "+str("%.6g"%self.total_delta_sum).rjust(13)+"\n")
            self.simulation_graph.flush()
        self.balls_graph_counter += 1

        for mesh_number in range(1, rigid_face_model_part.NumberOfMeshes()):
            if(rigid_face_model_part.GetMesh(mesh_number)[TOP]):

              self.top_mesh_nodes = rigid_face_model_part.GetMesh(mesh_number).Nodes

            if (self.rigid_graph_counter == self.graph_frequency):
                self.rigid_graph_counter = 0
                self.total_force_top = 0.0

                for node in self.top_mesh_nodes:

                  force_node_y = node.GetSolutionStepValue(ELASTIC_FORCES)[1]
                  self.total_force_top += force_node_y

                self.rigid_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_top).rjust(13)+"\n")
                self.rigid_graph.flush()
        self.rigid_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):
        error1, error2, error3 = self.compute_errors(self.output_filename)
        error4, error5, error6 = self.compute_rigid_errors(self.rigid_face_file)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 27:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 27 SUCCESSFUL (spheres)\n")
            shutil.rmtree('benchmark27_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 27 FAILED (spheres)\n")
        error_file.write("DEM Benchmark 27:")
        if (error4 < 10.0 and error5 < 10.0 and error6 < 10.0):
            error_file.write(" OK!........ Test 27 SUCCESSFUL (finite elements)\n")
        else:
            error_file.write(" KO!........ Test 27 FAILED (finite elements)\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + '27' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as current_data:
            for line in current_data:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))   #segona component del vector ()
                i+=1
        dem_error1 = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            dem_error1+=fabs(i-j)
        dem_error1/=summation_of_analytics_data

        print("Error in total force at the reference particle =", 100*dem_error1,"%")

        i = 0
        with open('paper_data/reference_graph_benchmark' + '27' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[2]))
                i+=1
        i = 0
        with open(output_filename) as current_data:
            for line in current_data:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[2]))   #segona component del vector ()
                i+=1
        dem_error2 = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            dem_error2+=fabs(i-j)
        dem_error2/=summation_of_analytics_data

        print("Error in angular velocity at the reference particle =", 100*dem_error2,"%")


        i = 0
        with open('paper_data/reference_graph_benchmark' + '27' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[3]))
                i+=1
        i = 0
        with open(output_filename) as current_data:
            for line in current_data:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[3]))   #segona component del vector ()
                i+=1
        dem_error3 = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            dem_error3+=fabs(i-j)
        dem_error3/=summation_of_analytics_data

        print("Error in delta displacement at the reference particle =", 100*dem_error3,"%")

        error1 = 100*dem_error1
        error2 = 100*dem_error2
        error3 = 100*dem_error3

        return error1, error2, error3

    def compute_rigid_errors(self, rigid_face_file):
        reference_data = lines_FEM = list(range(0, 1000));
        analytics_data = []; FEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark_rigid' + '27' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(rigid_face_file) as current_data:
            for line in current_data:
                if i in lines_FEM:
                    parts = line.split()
                    FEM_data.append(float(parts[1]))   #segona component del vector ()
                i+=1
        final_error = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(FEM_data, analytics_data):
            final_error+=fabs(i-j)
        final_error/=summation_of_analytics_data

        print("Error in FEM axial force =", 100*final_error,"%")

        error4 = 100*final_error

        error5 = error6 = 0

        return error4, error5, error6

    def create_gnuplot_scripts(self, output_filename, dt):
        pass



class Benchmark28:   #pendulo3D

    def __init__(self):
        self.number = 28
        self.generated_data = None
        self.balls_graph_counter = 1
        self.rigid_graph_counter = 1

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration):

        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.rigid_face_file = "benchmark" + str(sys.argv[1]) + '_rigid_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')
        self.rigid_graph = open(self.rigid_face_file, 'w')

    def get_final_data(self, modelpart, rigid_face_model_part, cluster_model_part):
        self.simulation_graph.close()
        self.rigid_graph.close()

    def cross_product(self, a, b):
        c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
        return c

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

        ang_vel = 20 * pi
        angular_velocity = [0, 0, ang_vel]
        rotation_matrix = [[cos(ang_vel * time), -1.0 * sin(ang_vel * time), 0], [sin(ang_vel * time), cos(ang_vel * time), 0], [0,0,1]]
        time_dt = time - dt
        rotation_matrix_minus_dt = [[cos(ang_vel * time_dt), -1.0 * sin(ang_vel * time_dt), 0], [sin(ang_vel * time_dt), cos(ang_vel * time_dt), 0], [0,0,1]] #
        centroid = [-1.0, 0.0, 0.0]
        relative_initial_node_coords, relative_node_coords, relative_node_coords_dt = [0]*3, [0]*3, [0]*3
        sum, sum_dt = 0, 0

        for node in modelpart.Nodes:
            if node.Id == 999999:
                for j in range(3):
                    rot_mat = rotation_matrix[j]
                    rot_mat_dt = rotation_matrix_minus_dt[j]
                    relative_initial_node_coords[0] = node.X0 - centroid[0]
                    relative_initial_node_coords[1] = node.Y0 - centroid[1]
                    relative_initial_node_coords[2] = node.Z0 - centroid[2]
                    for i in range(3):
                        sum += rot_mat[i] * relative_initial_node_coords[i]
                        sum_dt += rot_mat_dt[i] * relative_initial_node_coords[i]
                    relative_node_coords[j], sum, relative_node_coords_dt[j], sum_dt = sum, 0, sum_dt, 0
                node.X = relative_node_coords[0] + centroid[0]
                node.Y = relative_node_coords[1] + centroid[1]
                node.Z = relative_node_coords[2] + centroid[2]

                displacement = [0]*3
                displacement[0] = node.X - node.X0
                displacement[1] = node.Y - node.Y0
                displacement[2] = node.Z - node.Z0
                node.SetSolutionStepValue(DISPLACEMENT, displacement)

                velocity = [0]*3
                velocity = self.cross_product(angular_velocity, relative_node_coords)
                node.SetSolutionStepValue(VELOCITY, velocity)

                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

                delta_displacement = [0]*3
                delta_displacement[0] = relative_node_coords[0] - relative_node_coords_dt[0]
                delta_displacement[1] = relative_node_coords[1] - relative_node_coords_dt[1]
                delta_displacement[2] = relative_node_coords[2] - relative_node_coords_dt[2]
                node.SetSolutionStepValue(DELTA_DISPLACEMENT, delta_displacement)

                particle_rotation_angle = [0]*3
                particle_rotation_angle[0] = angular_velocity[0] * time
                particle_rotation_angle[1] = angular_velocity[1] * time
                particle_rotation_angle[2] = angular_velocity[2] * time
                node.SetSolutionStepValue(PARTICLE_ROTATION_ANGLE, particle_rotation_angle)

                delta_rotation = [0]*3
                delta_rotation[0] = angular_velocity[0] * dt
                delta_rotation[1] = angular_velocity[1] * dt
                delta_rotation[2] = angular_velocity[2] * dt
                node.SetSolutionStepValue(DELTA_ROTATION, delta_rotation)

                if time > 3.8e-5:
                    radius = 1.0001
                    node.SetSolutionStepValue(RADIUS, radius)
            if node.Id == 99999:
                angular_velocity = [0]*3
                node.SetSolutionStepValue(ANGULAR_VELOCITY, angular_velocity)

    def generate_graph_points(self, modelpart, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):

        #self.graph_frequency = int(5e-7/dt)   #graph_print_interval/dt
        self.graph_frequency = int(graph_print_interval/1/dt)   #1 veces mas grf que bin
        if self.graph_frequency < 1:
           self.graph_frequency = 1

        if (self.balls_graph_counter == self.graph_frequency):
            self.balls_graph_counter = 0
            self.total_force_x = 0.0
            self.total_force_y = 0.0
            self.total_force_z = 0.0
            self.total_force_sum = 0.0

            self.total_angular_x = 0.0
            self.total_angular_y = 0.0
            self.total_angular_z = 0.0
            self.total_angular_sum = 0.0

            self.total_delta_x = 0.0
            self.total_delta_y = 0.0
            self.total_delta_z = 0.0
            self.total_delta_sum = 0.0

            for node in modelpart.Nodes:
                if node.Id == 107:
                   force_node_x = node.GetSolutionStepValue(LOCAL_CONTACT_FORCE)[0]
                   force_node_y = node.GetSolutionStepValue(LOCAL_CONTACT_FORCE)[1]
                   force_node_z = node.GetSolutionStepValue(LOCAL_CONTACT_FORCE)[2]
                   self.total_force_x += force_node_x
                   self.total_force_y += force_node_y
                   self.total_force_z += force_node_z

                   angular_node_x = node.GetSolutionStepValue(ANGULAR_VELOCITY)[0]
                   angular_node_y = node.GetSolutionStepValue(ANGULAR_VELOCITY)[1]
                   angular_node_z = node.GetSolutionStepValue(ANGULAR_VELOCITY)[2]
                   self.total_angular_x += angular_node_x
                   self.total_angular_y += angular_node_y
                   self.total_angular_z += angular_node_z

                   delta_node_x = node.GetSolutionStepValue(DELTA_DISPLACEMENT)[0]
                   delta_node_y = node.GetSolutionStepValue(DELTA_DISPLACEMENT)[1]
                   delta_node_z = node.GetSolutionStepValue(DELTA_DISPLACEMENT)[2]
                   self.total_delta_x += delta_node_x
                   self.total_delta_y += delta_node_y
                   self.total_delta_z += delta_node_z

            self.total_force_sum = self.total_force_x + self.total_force_y + self.total_force_z
            self.total_angular_sum = self.total_angular_x + self.total_angular_y + self.total_angular_z
            self.total_delta_sum = self.total_delta_x + self.total_delta_y + self.total_delta_z
            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%self.total_force_sum).rjust(13)+" "+str("%.6g"%self.total_angular_sum).rjust(13)+" "+str("%.6g"%self.total_delta_sum).rjust(13)+"\n")
            self.simulation_graph.flush()
        self.balls_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 28:")

        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 28 SUCCESSFUL (spheres)\n")
            shutil.rmtree('benchmark28_Post_Files', ignore_errors = True)
        else:
            error_file.write(" KO!........ Test 28 FAILED (spheres)\n")
        error_file.write("DEM Benchmark 28:")


    def compute_errors(self, output_filename):
        reference_data = lines_DEM = list(range(0, 1000));
        analytics_data = []; DEM_data = []; summation_of_analytics_data = 0
        i = 0
        with open('paper_data/reference_graph_benchmark' + '28' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(output_filename) as current_data:
            for line in current_data:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))   #1 component del vector ()
                i+=1
        dem_error1 = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            dem_error1+=fabs(i-j)
        dem_error1/=summation_of_analytics_data

        print("Error in total force at the reference particle =", 100*dem_error1,"%")

        i = 0
        with open('paper_data/reference_graph_benchmark' + '28' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[2]))
                i+=1
        i = 0
        with open(output_filename) as current_data:
            for line in current_data:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[2]))   #segona component del vector ()
                i+=1
        dem_error2 = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            dem_error2+=fabs(i-j)
        dem_error2/=summation_of_analytics_data

        print("Error in angular velocity at the reference particle =", 100*dem_error2,"%")


        i = 0
        with open('paper_data/reference_graph_benchmark' + '28' + '.dat') as reference:
            for line in reference:
                if i in reference_data:
                    parts = line.split()
                    analytics_data.append(float(parts[3]))
                i+=1
        i = 0
        with open(output_filename) as current_data:
            for line in current_data:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[3]))   #3 component del vector ()
                i+=1
        dem_error3 = 0

        for j in analytics_data:
            summation_of_analytics_data+=abs(j)

        for i, j in zip(DEM_data, analytics_data):
            dem_error3+=fabs(i-j)
        dem_error3/=summation_of_analytics_data

        print("Error in delta displacement at the reference particle =", 100*dem_error3,"%")

        error1 = 100*dem_error1
        error2 = 100*dem_error2
        error3 = 100*dem_error3

        return error1, error2, error3

    def compute_rigid_errors(self, rigid_face_file):
        pass

    def create_gnuplot_scripts(self, output_filename, dt):
        pass



class Benchmark30: ########## Cylinder with imposed angular velocity (Velocity Verlet + Zhao)

    def __init__(self):
        self.number = 30

        self.cluster_graph_counter = 1   # deberia ser self.cluster_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.local_angular_velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.local_angular_velocity_list_outfile_name, 'w')

    def get_final_data(self, spheres_model_part, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, spheres_model_part, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.cluster_graph_counter == self.graph_frequency):     #if(self.cluster_graph_counter == self.graph_frequency):
            self.cluster_graph_counter = 0
            total_local_angular_velocity_x = 0.0
            total_local_angular_velocity_y = 0.0
            total_local_angular_velocity_z = 0.0

            for node in cluster_model_part.Nodes:
                current_local_angular_velocity_x = node.GetSolutionStepValue(LOCAL_ANGULAR_VELOCITY_X)
                total_local_angular_velocity_x += current_local_angular_velocity_x
                current_local_angular_velocity_y = node.GetSolutionStepValue(LOCAL_ANGULAR_VELOCITY_Y)
                total_local_angular_velocity_y += current_local_angular_velocity_y
                current_local_angular_velocity_z = node.GetSolutionStepValue(LOCAL_ANGULAR_VELOCITY_Z)
                total_local_angular_velocity_z += current_local_angular_velocity_z

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_local_angular_velocity_x).rjust(13)+" "+str("%.6g"%total_local_angular_velocity_y).rjust(13)+" "+str("%.6g"%total_local_angular_velocity_z).rjust(13)+"\n")
        self.cluster_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.local_angular_velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("===== DISCONTINUUM CLUSTERS TESTS =====\n\n")
        error_file.write("DEM Benchmark 30:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 30 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 30 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):  #FINALIZATION STEP

        lines_analytics = lines_DEM = list(range(0, 50));
        ref_data1 = []; ref_data2 = []; DEM_data1 = []; ref_data3 = []; DEM_data1 = []; DEM_data2 = []; DEM_data3 = []; summation_of_ref_data1 = 0; summation_of_ref_data2 = 0; summation_of_ref_data3 = 0
        i = 0
        with open('paper_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:  #with open('paper_data/reference_graph_benchmark30.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    ref_data1.append(float(parts[1]))
                    ref_data2.append(float(parts[2]))
                    ref_data3.append(float(parts[3]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data1.append(float(parts[1]))
                    DEM_data2.append(float(parts[2]))
                    DEM_data3.append(float(parts[3]))
                i+=1
        final_local_angular_velocity_x_error = 0
        final_local_angular_velocity_y_error = 0
        final_local_angular_velocity_z_error = 0

        for j in ref_data1:
            summation_of_ref_data1+=abs(j)
        for k in ref_data2:
            summation_of_ref_data2+=abs(k)
        for l in ref_data3:
            summation_of_ref_data3+=abs(l)

        for i, j in zip(DEM_data1, ref_data1):
            final_local_angular_velocity_x_error+=fabs(i-j)
        final_local_angular_velocity_x_error/=summation_of_ref_data1

        for k, l in zip(DEM_data2, ref_data2):
            final_local_angular_velocity_y_error+=fabs(k-l)
        final_local_angular_velocity_y_error/=summation_of_ref_data2

        for m, n in zip(DEM_data3, ref_data3):
            final_local_angular_velocity_z_error+=fabs(m-n)
        final_local_angular_velocity_z_error/=summation_of_ref_data3

        print("Error in local angular velocity X =", 100*final_local_angular_velocity_x_error,"%")

        print("Error in local angular velocity Y =", 100*final_local_angular_velocity_y_error,"%")

        print("Error in local angular velocity Z =", 100*final_local_angular_velocity_z_error,"%")

        error1 = 100*final_local_angular_velocity_x_error

        error2 = 100*final_local_angular_velocity_y_error

        error3 = 100*final_local_angular_velocity_z_error

        return error1, error2, error3

class Benchmark31: ########## Cylinder with imposed angular velocity (Symplectic Euler + Runge-Kutta)

    def __init__(self):
        self.number = 31

        self.cluster_graph_counter = 1   # deberia ser self.cluster_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.local_angular_velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.local_angular_velocity_list_outfile_name, 'w')

    def get_final_data(self, spheres_model_part, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, spheres_model_part, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.cluster_graph_counter == self.graph_frequency):     #if(self.cluster_graph_counter == self.graph_frequency):
            self.cluster_graph_counter = 0
            total_local_angular_velocity_x = 0.0
            total_local_angular_velocity_y = 0.0
            total_local_angular_velocity_z = 0.0

            for node in cluster_model_part.Nodes:
                current_local_angular_velocity_x = node.GetSolutionStepValue(LOCAL_ANGULAR_VELOCITY_X)
                total_local_angular_velocity_x += current_local_angular_velocity_x
                current_local_angular_velocity_y = node.GetSolutionStepValue(LOCAL_ANGULAR_VELOCITY_Y)
                total_local_angular_velocity_y += current_local_angular_velocity_y
                current_local_angular_velocity_z = node.GetSolutionStepValue(LOCAL_ANGULAR_VELOCITY_Z)
                total_local_angular_velocity_z += current_local_angular_velocity_z

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_local_angular_velocity_x).rjust(13)+" "+str("%.6g"%total_local_angular_velocity_y).rjust(13)+" "+str("%.6g"%total_local_angular_velocity_z).rjust(13)+"\n")
        self.cluster_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2, error3 = self.compute_errors(self.local_angular_velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 31:")

        if (error1 < 0.1 and error2 < 0.1 and error3 < 0.1):
            error_file.write(" OK!........ Test 31 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 31 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):  #FINALIZATION STEP

        lines_analytics = lines_DEM = list(range(0, 50));
        ref_data1 = []; ref_data2 = []; DEM_data1 = []; ref_data3 = []; DEM_data1 = []; DEM_data2 = []; DEM_data3 = []; summation_of_ref_data1 = 0; summation_of_ref_data2 = 0; summation_of_ref_data3 = 0
        i = 0
        with open('paper_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:  #with open('paper_data/reference_graph_benchmark31.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    ref_data1.append(float(parts[1]))
                    ref_data2.append(float(parts[2]))
                    ref_data3.append(float(parts[3]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data1.append(float(parts[1]))
                    DEM_data2.append(float(parts[2]))
                    DEM_data3.append(float(parts[3]))
                i+=1
        final_local_angular_velocity_x_error = 0
        final_local_angular_velocity_y_error = 0
        final_local_angular_velocity_z_error = 0

        for j in ref_data1:
            summation_of_ref_data1+=abs(j)
        for k in ref_data2:
            summation_of_ref_data2+=abs(k)
        for l in ref_data3:
            summation_of_ref_data3+=abs(l)

        for i, j in zip(DEM_data1, ref_data1):
            final_local_angular_velocity_x_error+=fabs(i-j)
        final_local_angular_velocity_x_error/=summation_of_ref_data1

        for k, l in zip(DEM_data2, ref_data2):
            final_local_angular_velocity_y_error+=fabs(k-l)
        final_local_angular_velocity_y_error/=summation_of_ref_data2

        for m, n in zip(DEM_data3, ref_data3):
            final_local_angular_velocity_z_error+=fabs(m-n)
        final_local_angular_velocity_z_error/=summation_of_ref_data3

        print("Error in local angular velocity X =", 100*final_local_angular_velocity_x_error,"%")

        print("Error in local angular velocity Y =", 100*final_local_angular_velocity_y_error,"%")

        print("Error in local angular velocity Z =", 100*final_local_angular_velocity_z_error,"%")

        error1 = 100*final_local_angular_velocity_x_error

        error2 = 100*final_local_angular_velocity_y_error

        error3 = 100*final_local_angular_velocity_z_error

        return error1, error2, error3

class Benchmark32: ########## Fiber cluster bouncing without any damping (Velocity Verlet + Zhao scheme)

    def __init__(self):
        self.number = 32

        self.cluster_graph_counter = 1   # deberia ser self.cluster_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.velocity_list_outfile_name, 'w')

    def get_final_data(self, spheres_model_part, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, spheres_model_part, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.cluster_graph_counter == self.graph_frequency):     #if(self.cluster_graph_counter == self.graph_frequency):
            self.cluster_graph_counter = 0
            total_velocity_z         = 0.0
            total_angular_velocity_y = 0.0

            for node in cluster_model_part.Nodes:
                current_velocity_z = node.GetSolutionStepValue(VELOCITY_Z)
                total_velocity_z += current_velocity_z
                current_angular_velocity_y = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)
                total_angular_velocity_y += current_angular_velocity_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_z).rjust(13)+" "+str("%.6g"%total_angular_velocity_y).rjust(13)+"\n")
        self.cluster_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2 = self.compute_errors(self.velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 32:")

        if (error1 < 0.1 and error2 < 0.1):
            error_file.write(" OK!........ Test 32 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 32 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):  #FINALIZATION STEP

        lines_analytics = lines_DEM = list(range(0, 100));
        ref_data1 = []; ref_data2 = []; DEM_data1 = []; DEM_data1 = []; DEM_data2 = []; summation_of_ref_data1 = 0; summation_of_ref_data2 = 0
        i = 0
        with open('paper_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:  #with open('paper_data/reference_graph_benchmark32.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    ref_data1.append(float(parts[1]))
                    ref_data2.append(float(parts[2]))
                i+=1
        i = 0
        with open(output_filename) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data1.append(float(parts[1]))
                    DEM_data2.append(float(parts[2]))
                i+=1
        final_velocity_z_error         = 0
        final_angular_velocity_y_error = 0

        for j in ref_data1:
            summation_of_ref_data1+=abs(j)
        for k in ref_data2:
            summation_of_ref_data2+=abs(k)

        for i, j in zip(DEM_data1, ref_data1):
            final_velocity_z_error+=fabs(i-j)
        final_velocity_z_error/=summation_of_ref_data1

        for k, l in zip(DEM_data2, ref_data2):
            final_angular_velocity_y_error+=fabs(k-l)
        final_angular_velocity_y_error/=summation_of_ref_data2

        print("Error in velocity Z =", 100*final_velocity_z_error,"%")

        print("Error in angular velocity Y =", 100*final_angular_velocity_y_error,"%")

        error1 = 100*final_velocity_z_error

        error2 = 100*final_angular_velocity_y_error

        return error1, error2

class Benchmark33: ########## Fiber cluster bouncing without any damping (Velocity Verlet + Runge-Kutta scheme)

    def __init__(self):
        self.number = 33

        self.cluster_graph_counter = 1   # deberia ser self.cluster_graph_counter = self.graph_frequency

    def set_initial_data(self, modelpart, rigid_face_model_part, iteration, number_of_points_in_the_graphic, coeff_of_restitution_iteration=0):

        self.velocity_list_outfile_name = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.velocity_list_outfile_name, 'w')

    def get_final_data(self, spheres_model_part, rigid_face_model_part, cluster_model_part):                 #FINALIZATION STEP

        self.simulation_graph.close()

    def ApplyNodalRotation(self, time, dt, modelpart):
        pass

    def generate_graph_points(self, spheres_model_part, rigid_face_model_part, cluster_model_part, time, graph_print_interval, dt):     #MAIN LOOP STEP

        self.graph_frequency        = int(graph_print_interval/dt)
        if self.graph_frequency < 1:
           self.graph_frequency = 1 #that means it is not possible to print results with a higher frequency than the computations delta time

        if(self.cluster_graph_counter == self.graph_frequency):     #if(self.cluster_graph_counter == self.graph_frequency):
            self.cluster_graph_counter = 0
            total_velocity_z         = 0.0
            total_angular_velocity_y = 0.0

            for node in cluster_model_part.Nodes:
                current_velocity_z = node.GetSolutionStepValue(VELOCITY_Z)
                total_velocity_z += current_velocity_z
                current_angular_velocity_y = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)
                total_angular_velocity_y += current_angular_velocity_y

            self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_velocity_z).rjust(13)+" "+str("%.6g"%total_angular_velocity_y).rjust(13)+"\n")
        self.cluster_graph_counter += 1

    def print_results(self, number_of_points_in_the_graphic, dt=0, elapsed_time=0.0):      #FINALIZATION STEP

        error1, error2 = self.compute_errors(self.velocity_list_outfile_name)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("DEM Benchmark 33:")

        if (error1 < 0.1 and error2 < 0.1):
            error_file.write(" OK!........ Test 33 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 33 FAILED\n")
        error_file.close()

    def compute_errors(self, velocity_list_outfile_name):  #FINALIZATION STEP

        lines_analytics = lines_DEM = list(range(0, 100));
        ref_data1 = []; ref_data2 = []; DEM_data1 = []; DEM_data1 = []; DEM_data2 = []; summation_of_ref_data1 = 0; summation_of_ref_data2 = 0
        i = 0
        with open('paper_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:  #with open('paper_data/reference_graph_benchmark33.dat') as inf:
            for line in inf:
                if i in lines_analytics:
                    parts = line.split()
                    ref_data1.append(float(parts[1]))
                    ref_data2.append(float(parts[2]))
                i+=1
        i = 0
        with open(velocity_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data1.append(float(parts[1]))
                    DEM_data2.append(float(parts[2]))
                i+=1
        final_velocity_z_error         = 0
        final_angular_velocity_y_error = 0

        for j in ref_data1:
            summation_of_ref_data1+=abs(j)
        for k in ref_data2:
            summation_of_ref_data2+=abs(k)

        for i, j in zip(DEM_data1, ref_data1):
            final_velocity_z_error+=fabs(i-j)
        final_velocity_z_error/=summation_of_ref_data1

        for k, l in zip(DEM_data2, ref_data2):
            final_angular_velocity_y_error+=fabs(k-l)
        final_angular_velocity_y_error/=summation_of_ref_data2

        print("Error in velocity Z =", 100*final_velocity_z_error,"%")

        print("Error in angular velocity Y =", 100*final_angular_velocity_y_error,"%")

        error1 = 100*final_velocity_z_error

        error2 = 100*final_angular_velocity_y_error

        return error1, error2

def delete_archives():

    #.......................Removing extra files
    files_to_delete_list = glob('*.time')
    files_to_delete_list.extend(glob('*.dat'))
    files_to_delete_list.extend(glob('*.gp'))
    files_to_delete_list.extend(glob('*.txt'))
    files_to_delete_list.extend(glob('*.lst'))

    for to_erase_file in files_to_delete_list:
        os.remove(to_erase_file)

    #............Getting rid of unuseful folders
    folders_to_delete_list      = glob('*Data')
    folders_to_delete_list.extend(glob('*ists'))
    folders_to_delete_list.extend(glob('*ults'))
    folders_to_delete_list.extend(glob('*he__'))
    folders_to_delete_list.extend(glob('*aphs'))
    folders_to_delete_list.extend(glob('*iles'))

    for to_erase_folder in folders_to_delete_list:
        shutil.rmtree(to_erase_folder)

def print_gnuplot_files_on_screen(gnuplot_script_name):
    system('gnuplot -persist ' + gnuplot_script_name)

def create_pdf_document(pdf_script_name):
    system('gnuplot -persist ' + pdf_script_name)
