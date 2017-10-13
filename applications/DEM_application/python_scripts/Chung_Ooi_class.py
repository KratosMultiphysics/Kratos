from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *                                  # importing the Kratos Library
from KratosMultiphysics.DEMApplication import *
CheckForPreviousImport()                                          # check that KratosMultiphysics was imported in the main script
import shutil
from glob import glob
from math import log, pi, sin, cos, tan, atan, fabs
from os import system

def initialize_time_parameters(benchmark_number):
    
    if benchmark_number==1:
        
        final_time                      = 0.0005
        dt                              = 6.4e-8 # Complies Rayleigh's condition
        output_time_step                = 0.000005
        number_of_points_in_the_graphic = 6
            
    elif benchmark_number==2:

        final_time                      = 0.007
        dt                              = 3e-6 #3e-7 # Complies Rayleigh's condition????????????????
        output_time_step                = 0.0001
        number_of_points_in_the_graphic = 6
            
    elif benchmark_number==3:
        
        final_time                      = 0.00031
        dt                              = 5.0e-8  #1.1e-9 # Complies Rayleigh's condition
        output_time_step                = 0.000001
        number_of_points_in_the_graphic = 6
        
    elif benchmark_number==4:
        
        final_time                      = 0.0002  #0.00003
        dt                              = 1.9e-8  #1.9e-9 # Complies Rayleigh's condition
        output_time_step                = 0.000001
        number_of_points_in_the_graphic = 17
                
    elif benchmark_number==5:
                
        final_time                      = 0.0000005
        dt                              = 3.3e-10  #3.6e-12 # Complies Rayleigh's condition
        output_time_step                = 0.00000005
        number_of_points_in_the_graphic = 17
                
    elif benchmark_number==6:
                
        final_time                      = 0.01
        dt                              = 5.0e-7  #1.0e-7 # Complies Rayleigh's condition ????????????????
        output_time_step                = 0.00025
        number_of_points_in_the_graphic = 17
        
    elif benchmark_number==7:
                
        final_time                      = 0.0005
        dt                              = 2e-7 #4.4614e-8 # Complies Rayleigh's condition ????????????????
        output_time_step                = 0.000005
        number_of_points_in_the_graphic = 17
                
    elif benchmark_number==8:
        
        final_time                      = 0.02
        dt                              = 1.0e-6 #5.0e-7 # Complies Rayleigh's condition
        output_time_step                = 0.0001
        number_of_points_in_the_graphic = 17
                
    elif benchmark_number==9:
                
        final_time                      = 0.001 #0.0005
        dt                              = 1e-7 #6.4e-8 # Complies Rayleigh's condition
        output_time_step                = 0.000005
        number_of_points_in_the_graphic = 6
            
    else: #benchmark_number=10
                
        final_time                      = 0.001 #0.0005
        dt                              = 6.4e-8 # Complies Rayleigh's condition
        output_time_step                = 0.000005
        number_of_points_in_the_graphic = 10
                
    return final_time, dt, output_time_step, number_of_points_in_the_graphic


class Benchmark1:
    
    def __init__(self):
        self.initial_normal_vel = 10.0
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
                   
        for node in modelpart.Nodes:
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X, -self.initial_normal_vel)
            else:
                node.SetSolutionStepValue(VELOCITY_X,  self.initial_normal_vel)

    def get_final_data(self, modelpart):
        pass
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        normal_contact_force_outfile_name = 'variables_for_node_1.txt'
        gnuplot_script_name = 'benchmark1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + normal_contact_force_outfile_name + "' every 20 u 1:8 w lp lt -1 lw 1.5 ps 1 pt 4")
        self.gnuplot_outfile.close()
        print_gnuplot_files_on_screen(gnuplot_script_name)
            

class Benchmark2:
    
    def __init__(self):
        self.initial_normal_vel = -0.2
    
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
                   
        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)
            
    def get_final_data(self, modelpart):
        pass
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        normal_contact_force_outfile_name = 'variables_for_node_1.txt'
        gnuplot_script_name = 'benchmark2_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + normal_contact_force_outfile_name + "' every 10 u 1:10 w lp lt 3 lw 1.5 ps 1 pt 6")
        self.gnuplot_outfile.close()
        print_gnuplot_files_on_screen(gnuplot_script_name)
        

class Benchmark3:
    
    def __init__(self):
        self.restitution_numbers_list = []
        self.initial_normal_vel = 0
        self.restitution_numbers_vector_list_outfile = None
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
                   
        number = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)
        
        if number_of_points_in_the_graphic == 1:
            number = 0
        else:
            number = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)
        
        for node in modelpart.Nodes:
            self.initial_normal_vel = node.GetSolutionStepValue(VELOCITY_Z)
            modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = number

    def get_final_data(self, modelpart):
     
        for node in modelpart.Nodes:
            final_vel = node.GetSolutionStepValue(VELOCITY_Z)
        
        restitution_coefficient = -final_vel / self.initial_normal_vel
        self.restitution_numbers_list.append(restitution_coefficient)
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        self.restitution_numbers_vector_list_outfile_name = "benchmark3_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat'
        self.restitution_numbers_vector_list_outfile = open(self.restitution_numbers_vector_list_outfile_name, 'w')
    
        for i in range(0, number_of_points_in_the_graphic):
            first_col = 1/(number_of_points_in_the_graphic-1) * i
            self.restitution_numbers_vector_list_outfile.write("%6.4f %11.8f" % (first_col, self.restitution_numbers_list[i]) + '\n')
        self.restitution_numbers_vector_list_outfile.close()
        
        gnuplot_script_name = 'benchmark3_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + self.restitution_numbers_vector_list_outfile_name + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4, '"\
                                                      + self.restitution_numbers_vector_list_outfile_name + "' u 1:3 w lp lt 2 lw 1.5 ps 2 pt 6")
        self.gnuplot_outfile.close()
        
        #self.create_gnuplot_scripts(self.restitution_numbers_vector_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.restitution_numbers_vector_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 3:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 3 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3 FAILED\n")
        error_file.close()
        
    def create_gnuplot_scripts(self, restitution_numbers_vector_list_outfile_name, dt):    
        
        gnuplot_script_name_1 = 'benchmark3_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Coefficient of restitution'\nset ylabel 'Damping ratio'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:1][0:1] '" + restitution_numbers_vector_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark3_graph1.dat' w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark3_graph1.dat' w lp ls 2 t 'Cast iron'\n")
        self.gnuplot_outfile.close()
                
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
                
    def compute_errors(self, restitution_numbers_vector_list_outfile_name):
        
        lines_Chung = lines_DEM = list(range(0, 6));
        Chung_data = []; DEM_data = []; summation_of_Chung_data = 0 
        i = 0
        with open('paper_data/benchmark3_graph1.dat') as inf:
            for line in inf:
                if i in lines_Chung:
                    parts = line.split()
                    Chung_data.append(float(parts[1]))
                i+=1
        i = 0
        with open(restitution_numbers_vector_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_restitution_numbers_error = 0
        
        for j in Chung_data:
            summation_of_Chung_data+=abs(j)
        
        for i, j in zip(DEM_data, Chung_data):
            final_restitution_numbers_error+=fabs(i-j)
        final_restitution_numbers_error/=summation_of_Chung_data
        
        print("Error in restitution numbers =", 100*final_restitution_numbers_error,"%")
        
        error1 = 100*final_restitution_numbers_error
        
        error2 = error3 = 0
        
        return error1, error2, error3
                

class Benchmark4:
    
    def __init__(self):
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
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
        
        self.degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  self.initial_module_vel * sin(self.degrees * pi / 180.0)
        initial_normal_vel = -self.initial_module_vel * cos(self.degrees * pi / 180.0)
                
        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.initial_tangential_vel)
            node.SetSolutionStepValue(VELOCITY_Z, initial_normal_vel)
            
    def get_final_data(self, modelpart):
     
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
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
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
        
        #self.create_gnuplot_scripts(self.tangential_restitution_coefficient_list_outfile_name, self.final_angular_vel_list_outfile_name,\
        #                            self.rebound_angle_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.tangential_restitution_coefficient_list_outfile_name, self.final_angular_vel_list_outfile_name,\
                                    self.rebound_angle_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 4:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 4 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 4 FAILED\n")                           
        error_file.close()
        
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
        
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        print_gnuplot_files_on_screen(gnuplot_script_name_3)
        
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
        self.initial_normal_vel = -5.0
        self.initial_tangential_vel = 0
        self.radius = 0.00001
        self.Vst_div_mu_per_Vcn_list = []
        self.Vst_prima_div_mu_per_Vcn_prima_list = []
        self.r_w1_prima_div_mu_per_Vcn_list = []
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile = None
        self.r_w1_prima_div_mu_per_Vcn_list_outfile = None
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
        
        degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  -self.initial_normal_vel * tan(degrees * pi / 180.0)
                
        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Y, self.initial_tangential_vel)
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)
            
    def get_final_data(self, modelpart):
     
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
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name = "benchmark5_dt_" + str(dt) + '_Vst_prima_div_mu_per_Vcn_prima_list_data.dat'
        self.r_w1_prima_div_mu_per_Vcn_list_outfile_name = "benchmark5_dt_" + str(dt) + '_r_w1_prima_div_mu_per_Vcn_list_data.dat'
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile = open(self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, 'w')
        self.r_w1_prima_div_mu_per_Vcn_list_outfile = open(self.r_w1_prima_div_mu_per_Vcn_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.Vst_prima_div_mu_per_Vcn_prima_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_mu_per_Vcn_list[i], self.Vst_prima_div_mu_per_Vcn_prima_list[i]) + '\n')
            self.r_w1_prima_div_mu_per_Vcn_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_mu_per_Vcn_list[i], self.r_w1_prima_div_mu_per_Vcn_list[i]) + '\n')
        self.Vst_prima_div_mu_per_Vcn_prima_list_outfile.close()
        self.r_w1_prima_div_mu_per_Vcn_list_outfile.close()
        
        #self.create_gnuplot_scripts(self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, self.r_w1_prima_div_mu_per_Vcn_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.Vst_prima_div_mu_per_Vcn_prima_list_outfile_name, self.r_w1_prima_div_mu_per_Vcn_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 5:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 5 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 5 FAILED\n")
        error_file.close()    
                
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
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt 3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:20][-6:0] '" + r_w1_prima_div_mu_per_Vcn_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph2.dat' index 0 w lp ls 1 t 'Steel',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph2.dat' index 1 w lp ls 2 t 'Polyethylene',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark5_graph2.dat' index 2 w p pt 7 ps 2 lt -1 t 'FEM'\n")
        self.gnuplot_outfile.close()
                          
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        
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
        self.initial_normal_vel = -0.2
        self.initial_tangential_vel = 0
        self.radius = 0.1
        self.special_quantity_list = []
        self.beta_list = []
        self.Vst_div_Vcn_list = []
        self.Vst_prima_div_Vcn_prima_list = []
        self.beta_list_outfile = None
        self.Vst_prima_div_Vcn_prima_list_outfile = None
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
        
        degrees = 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel = -self.initial_normal_vel * tan(degrees * pi / 180.0) # Here is tangential of the contact point, only. In X axis
        initial_angular_vel = -self.initial_tangential_vel / self.radius # In Y axis
                        
        for node in modelpart.Nodes:
            node.SetSolutionStepValue(VELOCITY_Z, self.initial_normal_vel)
            node.SetSolutionStepValue(ANGULAR_VELOCITY_Y, initial_angular_vel)        
            
    def get_final_data(self, modelpart):
     
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
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        self.beta_list_outfile_name = "benchmark6_dt_" + str(dt) + '_beta_list_data.dat'
        self.Vst_prima_div_Vcn_prima_list_outfile_name = "benchmark6_dt_" + str(dt) + '_Vst_prima_div_Vcn_prima_data.dat'
        self.beta_list_outfile = open(self.beta_list_outfile_name, 'w')
        self.Vst_prima_div_Vcn_prima_list_outfile = open(self.Vst_prima_div_Vcn_prima_list_outfile_name, 'w')

        for i in range(0, number_of_points_in_the_graphic):
            self.beta_list_outfile.write("%14.8f %14.8f" % (self.special_quantity_list[i], self.beta_list[i]) + '\n')
            self.Vst_prima_div_Vcn_prima_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_Vcn_list[i], self.Vst_prima_div_Vcn_prima_list[i]) + '\n')
        self.beta_list_outfile.close()
        self.Vst_prima_div_Vcn_prima_list_outfile.close()
        
        #self.create_gnuplot_scripts(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 6:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 6 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 6 FAILED\n")
        error_file.close()
        
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
                          
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        
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
        self.initial_angular_vel = 0
        self.final_tangential_center_vel_list_outfile = None
        self.final_angular_vel_list_outfile = None
        self.initial_angular_vel_list = []
        self.final_tangential_center_vel_list = []
        self.final_angular_vel_list = []
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
        
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

    def get_final_data(self, modelpart):
        
        mu = 0.4
        
        for node in modelpart.Nodes:
            if node.Id == 1:
                final_tangential_center_velocity = node.GetSolutionStepValue(VELOCITY_Z)
                final_angular_vel = node.GetSolutionStepValue(ANGULAR_VELOCITY_Y)
                        
        self.initial_angular_vel_list.append(self.initial_angular_vel)
        self.final_tangential_center_vel_list.append(final_tangential_center_velocity)
        self.final_angular_vel_list.append(final_angular_vel)
        
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
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
        
        #self.create_gnuplot_scripts(self.final_tangential_center_vel_list_outfile_name, self.final_angular_vel_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.final_tangential_center_vel_list_outfile_name, self.final_angular_vel_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 7:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 7 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 7 FAILED\n")
        error_file.close()    
        
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
                          
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
        
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
        self.initial_normal_vel = 0.2
        self.initial_tangential_vel = 0
        self.radius = 0.1
        self.special_quantity_list = []
        self.beta_list = []
        self.Vst_div_Vcn_list = []
        self.Vst_prima_div_Vcn_prima_list = []
        self.beta_list_outfile = None
        self.Vst_prima_div_Vcn_prima_list_outfile = None
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
        
        degrees = 90 - 90 / (number_of_points_in_the_graphic + 1) * iteration
        self.initial_tangential_vel =  self.initial_normal_vel * tan(degrees * pi / 180.0) # Here is tangential of the contact point, only
        initial_angular_vel    =  -self.initial_tangential_vel / self.radius
                        
        for node in modelpart.Nodes:
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X, self.initial_normal_vel)
                node.SetSolutionStepValue(ANGULAR_VELOCITY_Y, initial_angular_vel)        
                    
    def get_final_data(self, modelpart):
     
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
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        self.beta_list_outfile_name = 'benchmark8_dt_' + str(dt) + 's_beta_list_data.dat'
        self.Vst_prima_div_Vcn_prima_list_outfile_name = 'benchmark8_dt_' + str(dt) + 's_Vst_prima_div_Vcn_prima_list_data.dat'
        self.beta_list_outfile = open(self.beta_list_outfile_name, 'w')
        self.Vst_prima_div_Vcn_prima_list_outfile = open(self.Vst_prima_div_Vcn_prima_list_outfile_name, 'w')
        
        for i in range(0, number_of_points_in_the_graphic):
            self.beta_list_outfile.write("%14.8f %14.8f" % (self.special_quantity_list[i], self.beta_list[i]) + '\n')
            self.Vst_prima_div_Vcn_prima_list_outfile.write("%14.8f %14.8f" % (self.Vst_div_Vcn_list[i], self.Vst_prima_div_Vcn_prima_list[i]) + '\n')
        
        self.beta_list_outfile.close()
        self.Vst_prima_div_Vcn_prima_list_outfile.close()
        
        #self.create_gnuplot_scripts(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.beta_list_outfile_name, self.Vst_prima_div_Vcn_prima_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 8:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 8 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 8 FAILED\n")
        error_file.close()
        
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
                          
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
        print_gnuplot_files_on_screen(gnuplot_script_name_2)
                
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
        self.initial_normal_vel = 200.0
        self.restitution_numbers_list = []
        self.restitution_numbers_vector_list_outfile = None
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
                   
        if number_of_points_in_the_graphic == 1:
            number = 0
        else:
            number = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)
                
        for node in modelpart.Nodes:
            
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X,  self.initial_normal_vel)
                node.SetSolutionStepValue(VELOCITY_Z,  self.initial_normal_vel)
                modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = number
            else:
                node.SetSolutionStepValue(VELOCITY_X, -self.initial_normal_vel)
                node.SetSolutionStepValue(VELOCITY_Z,  self.initial_normal_vel)
                modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = number

    def get_final_data(self, modelpart):
        
        for node in modelpart.Nodes:
            if node.Id == 1:
                final_vel = node.GetSolutionStepValue(VELOCITY_X)
                           
        restitution_coefficient = -final_vel / self.initial_normal_vel
        self.restitution_numbers_list.append(restitution_coefficient)
    
    def print_results(self, number_of_points_in_the_graphic, dt=0):
        
        self.restitution_numbers_vector_list_outfile_name = "benchmark9_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat'
        self.restitution_numbers_vector_list_outfile = open(self.restitution_numbers_vector_list_outfile_name, 'w')
        
        for i in range(0, number_of_points_in_the_graphic):
            if number_of_points_in_the_graphic == 1:
                first_col = 0
            else:
                first_col = 1/(number_of_points_in_the_graphic-1) * i
            self.restitution_numbers_vector_list_outfile.write("%6.4f %6.4f %11.8f" % (first_col, first_col, self.restitution_numbers_list[i]) + '\n')
        self.restitution_numbers_vector_list_outfile.close()
                
        gnuplot_script_name = 'benchmark9_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name, 'w')
        self.gnuplot_outfile.write("set grid; plot '" + self.restitution_numbers_vector_list_outfile_name + "' u 1:2 w lp lt 3 lw 1.5 ps 2 pt 4, '"\
                                                      + self.restitution_numbers_vector_list_outfile_name + "' u 1:3 w lp lt 2 lw 1.5 ps 2 pt 6")
        self.gnuplot_outfile.close()
        
        #self.create_gnuplot_scripts(self.restitution_numbers_vector_list_outfile_name, dt)
        
        error1, error2, error3 = self.compute_errors(self.restitution_numbers_vector_list_outfile_name)
        
        error_filename = 'errors.txt'
        error_file = open(error_filename, 'a')
        error_file.write("Test 9:")
        
        if (error1 < 10.0 and error2 < 10.0 and error3 < 10.0):
            error_file.write(" OK!........ Test 9 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 9 FAILED\n")
        error_file.close()
        
    def create_gnuplot_scripts(self, restitution_numbers_vector_list_outfile_name, dt):    
        
        gnuplot_script_name_1 = 'benchmark9_comparison_1_dt_' + str(dt) + 's.gp'
        self.gnuplot_outfile = open(gnuplot_script_name_1, 'w')
        self.gnuplot_outfile.write("set grid\nset key left bottom\nset xlabel 'Coefficient of restitution'\nset ylabel 'Damping ratio'\nset style line 1 pt 8 lt -1 ps 3\nset style line 2 pt 9 lt  3 ps 3\n")
        self.gnuplot_outfile.write("plot [0:1][0:1] '" + restitution_numbers_vector_list_outfile_name + "' w lp lt 1 lw 1.5 ps 2 pt 5,\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark9_graph1.dat' w lp ls 1 t 'Al. oxide',\\\n")
        self.gnuplot_outfile.write("'paper_data/benchmark9_graph1.dat' w lp ls 2 t 'Cast iron'\n")
        self.gnuplot_outfile.close()
                
        print_gnuplot_files_on_screen(gnuplot_script_name_1)
                
    def compute_errors(self, restitution_numbers_vector_list_outfile_name):
        
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
        with open(restitution_numbers_vector_list_outfile_name) as inf:
            for line in inf:
                if i in lines_DEM:
                    parts = line.split()
                    DEM_data.append(float(parts[1]))
                i+=1
        final_restitution_numbers_error = 0
        
        for j in Chung_data:
            summation_of_Chung_data+=abs(j)
        
        for i, j in zip(DEM_data, Chung_data):
            final_restitution_numbers_error+=fabs(i-j)
        final_restitution_numbers_error/=summation_of_Chung_data
        
        print("Error in restitution numbers =", 100*final_restitution_numbers_error,"%")
        
        error1 = 100*final_restitution_numbers_error
        
        error2 = error3 = 0
        
        return error1, error2, error3                    

class Benchmark10:
    
    def __init__(self):
        self.initial_normal_vel = 10.0
        self.restitution_numbers_list = []
        self.restitution_numbers_vector_list_outfile = None
        
    def get_initial_data(self, modelpart, iteration, number_of_points_in_the_graphic):
                   
        if number_of_points_in_the_graphic == 1:
            coefficient_of_restitution = 0
        else:
            coefficient_of_restitution = 1.0/(number_of_points_in_the_graphic-1) * (iteration - 1)
                        
        for node in modelpart.Nodes:
            
            if node.Id == 1:
                node.SetSolutionStepValue(VELOCITY_X,  self.initial_normal_vel)
                node.SetSolutionStepValue(VELOCITY_Z,  self.initial_normal_vel)
                modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = coefficient_of_restitution
            else:
                node.SetSolutionStepValue(VELOCITY_X, -1.2 * self.initial_normal_vel)
                node.SetSolutionStepValue(VELOCITY_Z,  1.6 * self.initial_normal_vel)
                modelpart.GetProperties()[1][COEFFICIENT_OF_RESTITUTION] = coefficient_of_restitution

    def get_final_data(self, modelpart):
        
        for node in modelpart.Nodes:
            if node.Id == 1:
                final_vel = node.GetSolutionStepValue(VELOCITY_X)
                           
        restitution_coefficient = -final_vel / self.initial_normal_vel
        self.restitution_numbers_list.append(restitution_coefficient)
    
    def print_results(self, number_of_points_in_the_graphic, dt):
        
        self.restitution_numbers_vector_list_outfile = open("benchmark10_dt_" + str(dt) + '_restitution_numbers_vector_list_data.dat', 'w')
        
        for i in range(0, number_of_points_in_the_graphic):
            if number_of_points_in_the_graphic == 1:
                first_col = 0
            else:
                first_col = 1/(number_of_points_in_the_graphic-1) * i
            self.restitution_numbers_vector_list_outfile.write("%6.4f %6.4f %11.8f" % (first_col, first_col, self.restitution_numbers_list[i]) + '\n')
        self.restitution_numbers_vector_list_outfile.close()        
        

def delete_archives(nodeplotter):
    
    #.......................Removing extra files
    files_to_delete_list = glob('*.time')
    for to_erase_file in files_to_delete_list:
        os.remove(to_erase_file)
    
    #............Getting rid of unuseful folders    
    folders_to_delete_list      = glob('*Data')
    folders_to_delete_list.extend(glob('*ists'))
    folders_to_delete_list.extend(glob('*ults'))
    folders_to_delete_list.extend(glob('*he__'))
    folders_to_delete_list.extend(glob('*aphs'))
    
    for to_erase_folder in folders_to_delete_list:
        shutil.rmtree(to_erase_folder)

def print_gnuplot_files_on_screen(gnuplot_script_name):
    system('gnuplot -persist ' + gnuplot_script_name)
    
def create_pdf_document(pdf_script_name):
    system('gnuplot -persist ' + pdf_script_name)

