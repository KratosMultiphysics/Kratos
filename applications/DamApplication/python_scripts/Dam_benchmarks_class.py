from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam

KratosMultiphysics.CheckForPreviousImport() # check that KratosMultiphysics was imported in the main script

import shutil
from glob import glob
from math import pi, sin, cos, tan, atan, fabs
import sys
import os

class Benchmark201:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y

            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+"\n")
        self.simulation_graph.flush()        

    def print_results(self):
        error1, error2 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("Mechanical-PlaneStrain with hydrostatic pressure:\n\n")
        error_file.write("Dam Benchmark 2D 01:")

        if (error1 < 0.001 and error2 < 0.001):
            error_file.write(" OK!........ Test 2D 01 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 01 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX
        
        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY

        return error1, error2

    def finalize_data(self):
        self.simulation_graph.close()

class Benchmark202:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y

            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+"\n")
        self.simulation_graph.flush()        

    def print_results(self):
        error1, error2 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("Mechanical-PlaneStress with hydrostatic pressure:\n\n")
        error_file.write("Dam Benchmark 2D 02:")

        if (error1 < 0.001 and error2 < 0.001):
            error_file.write(" OK!........ Test 2D 02 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 02 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX
        
        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY

        return error1, error2

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark203:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStrain with hydrostatic pressure and heat flux:\n\n")
        error_file.write("Dam Benchmark 2D 03:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 03 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 03 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark204:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStress with hydrostatic pressure and heat flux:\n\n")
        error_file.write("Dam Benchmark 2D 04:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 04 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 04 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark205:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStrain with hydrostatic pressure and imposed temperature:\n\n")
        error_file.write("Dam Benchmark 2D 05:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 05 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 05 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark206:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStress with hydrostatic pressure and imposed temperature:\n\n")
        error_file.write("Dam Benchmark 2D 06:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 06 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 06 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark207:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStrain with hydrostatic pressure and Bofang:\n\n")
        error_file.write("Dam Benchmark 2D 07:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 07 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 07 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark208:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStress with hydrostatic pressure and Bofang:\n\n")
        error_file.write("Dam Benchmark 2D 08:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 08 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 08 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark209:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStrain with heat flux:\n\n")
        error_file.write("Dam Benchmark 2D 09:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 09 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 09 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark210:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStress with heat flux:\n\n")
        error_file.write("Dam Benchmark 2D 10:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 10 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 10 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()

class Benchmark211:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStrain with imposed temperature:\n\n")
        error_file.write("Dam Benchmark 2D 11:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 11 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 11 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()

class Benchmark212:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 1:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
            if node.Id == 3340:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        self.simulation_graph.flush()

    def print_results(self):
        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical-PlaneStress with imposed temperature:\n\n")
        error_file.write("Dam Benchmark 2D 12:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 2D 12 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 2D 12 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_temp.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_temp.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_temp:
            summation_of_reference_data_temp += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(m-n)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_temp

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()

class Benchmark301:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0

        for node in modelpart.Nodes:
            if node.Id == 186:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
                total_displacement_z += displacement_z
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_displacement_z).rjust(13)+"\n")
        
        self.simulation_graph.flush()

    def print_results(self):

        error1, error2, error3 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("Mechanical-PlaneStrain with hydrostatic pressure:\n\n")
        error_file.write("Dam Benchmark 3D 01:")

        if (error1 < 0.001 and error2 < 0.001 and error3 < 0.001):
            error_file.write(" OK!........ Test 3D 01 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3D 01 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_displZ = []; calculation_data_displZ = []; summation_of_reference_data_displZ = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_displZ.append(float(parts[3]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_displZ.append(float(parts[3]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_displZ = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_displZ:
            summation_of_reference_data_displZ += fabs(k)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_displZ, reference_data_displZ):
            generated_data_error_displZ += fabs(m-n)
        generated_data_error_displZ /= summation_of_reference_data_displZ

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in Z displacement =", 100 * generated_data_error_displZ, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_displZ

        return error1, error2, error3

    def finalize_data(self):
        self.simulation_graph.close()

class Benchmark302:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 186:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
                total_displacement_z += displacement_z
            if node.Id == 4050:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_displacement_z).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        
        self.simulation_graph.flush()

    def print_results(self):

        error1, error2, error3, error4 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical with hydrostatic pressure and heat flux:\n\n")
        error_file.write("Dam Benchmark 3D 02:")

        if (error1 < 1.0 and error2 < 1.0 and error3 < 1.0 and error4 < 1.0):
            error_file.write(" OK!........ Test 3D 02 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3D 02 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_displZ = []; calculation_data_displZ = []; summation_of_reference_data_displZ = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_displZ.append(float(parts[3]))
                reference_data_temp.append(float(parts[4]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_displZ.append(float(parts[3]))
                calculation_data_temp.append(float(parts[4]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_displZ = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_displZ:
            summation_of_reference_data_displZ += fabs(k)
        for l in reference_data_temp:
            summation_of_reference_data_temp += fabs(l)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_displZ, reference_data_displZ):
            generated_data_error_displZ += fabs(m-n)
        generated_data_error_displZ /= summation_of_reference_data_displZ
        
        for o, p in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(o-p)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in Z displacement =", 100 * generated_data_error_displZ, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_displZ
        error4 = 100 * generated_data_error_temp

        return error1, error2, error3, error4

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark303:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 186:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
                total_displacement_z += displacement_z
            if node.Id == 4050:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_displacement_z).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        
        self.simulation_graph.flush()

    def print_results(self):

        error1, error2, error3, error4 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical with hydrostatic pressure and imposed temperature:\n\n")
        error_file.write("Dam Benchmark 3D 03:")

        if (error1 < 1.0 and error2 < 1.0 and error3 < 1.0 and error4 < 1.0):
            error_file.write(" OK!........ Test 3D 03 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3D 03 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_displZ = []; calculation_data_displZ = []; summation_of_reference_data_displZ = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_displZ.append(float(parts[3]))
                reference_data_temp.append(float(parts[4]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_displZ.append(float(parts[3]))
                calculation_data_temp.append(float(parts[4]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_displZ = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_displZ:
            summation_of_reference_data_displZ += fabs(k)
        for l in reference_data_temp:
            summation_of_reference_data_temp += fabs(l)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_displZ, reference_data_displZ):
            generated_data_error_displZ += fabs(m-n)
        generated_data_error_displZ /= summation_of_reference_data_displZ
        
        for o, p in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(o-p)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in Z displacement =", 100 * generated_data_error_displZ, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_displZ
        error4 = 100 * generated_data_error_temp

        return error1, error2, error3, error4

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark304:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 186:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
                total_displacement_z += displacement_z
            if node.Id == 4050:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_displacement_z).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        
        self.simulation_graph.flush()

    def print_results(self):

        error1, error2, error3, error4 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical with hydrostatic pressure and Bofang:\n\n")
        error_file.write("Dam Benchmark 3D 04:")

        if (error1 < 1.0 and error2 < 1.0 and error3 < 1.0 and error4 < 1.0):
            error_file.write(" OK!........ Test 3D 04 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3D 04 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_displZ = []; calculation_data_displZ = []; summation_of_reference_data_displZ = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_displZ.append(float(parts[3]))
                reference_data_temp.append(float(parts[4]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_displZ.append(float(parts[3]))
                calculation_data_temp.append(float(parts[4]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_displZ = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_displZ:
            summation_of_reference_data_displZ += fabs(k)
        for l in reference_data_temp:
            summation_of_reference_data_temp += fabs(l)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_displZ, reference_data_displZ):
            generated_data_error_displZ += fabs(m-n)
        generated_data_error_displZ /= summation_of_reference_data_displZ
        
        for o, p in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(o-p)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in Z displacement =", 100 * generated_data_error_displZ, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_displZ
        error4 = 100 * generated_data_error_temp

        return error1, error2, error3, error4

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark305:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 186:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
                total_displacement_z += displacement_z
            if node.Id == 4050:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_displacement_z).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        
        self.simulation_graph.flush()

    def print_results(self):

        error1, error2, error3, error4 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical with heat flux:\n\n")
        error_file.write("Dam Benchmark 3D 05:")

        if (error1 < 1.0 and error2 < 1.0 and error3 < 1.0 and error4 < 1.0):
            error_file.write(" OK!........ Test 3D 05 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3D 05 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_displZ = []; calculation_data_displZ = []; summation_of_reference_data_displZ = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_displZ.append(float(parts[3]))
                reference_data_temp.append(float(parts[4]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_displZ.append(float(parts[3]))
                calculation_data_temp.append(float(parts[4]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_displZ = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_displZ:
            summation_of_reference_data_displZ += fabs(k)
        for l in reference_data_temp:
            summation_of_reference_data_temp += fabs(l)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_displZ, reference_data_displZ):
            generated_data_error_displZ += fabs(m-n)
        generated_data_error_displZ /= summation_of_reference_data_displZ
        
        for o, p in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(o-p)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in Z displacement =", 100 * generated_data_error_displZ, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_displZ
        error4 = 100 * generated_data_error_temp

        return error1, error2, error3, error4

    def finalize_data(self):
        self.simulation_graph.close()
        
class Benchmark306:

    def __init__(self):
        pass

    def set_initial_data(self):
        self.output_filename = "benchmark" + str(sys.argv[1]) + '_graph.dat'
        self.simulation_graph = open(self.output_filename, 'w')

    def generate_graph_points(self, modelpart, time):
        total_displacement_x = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_temperature = 0.0

        for node in modelpart.Nodes:
            if node.Id == 186:
                displacement_x = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                displacement_y = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                displacement_z = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                total_displacement_x += displacement_x
                total_displacement_y += displacement_y
                total_displacement_z += displacement_z
            if node.Id == 4050:
                temperature = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                total_temperature += temperature
            del node

        self.simulation_graph.write(str("%.8g"%time).rjust(12)+" "+str("%.6g"%total_displacement_x).rjust(13)+" "+str("%.6g"%total_displacement_y).rjust(13)+" "+str("%.6g"%total_displacement_z).rjust(13)+" "+str("%.6g"%total_temperature).rjust(13)+"\n")
        
        self.simulation_graph.flush()

    def print_results(self):

        error1, error2, error3, error4 = self.compute_errors(self.output_filename)

        error_filename = 'errors.err'
        error_file = open(error_filename, 'a')
        error_file.write("\n\n")
        error_file.write("ThermoMechanical with imposed temperature:\n\n")
        error_file.write("Dam Benchmark 3D 06:")

        if (error1 < 1.0 and error2 < 1.0 and error3 < 1.0 and error4 < 1.0):
            error_file.write(" OK!........ Test 3D 06 SUCCESSFUL\n")
        else:
            error_file.write(" KO!........ Test 3D 06 FAILED\n")
        error_file.close()

    def compute_errors(self, output_filename):
        reference_data_displX = []; calculation_data_displX = []; summation_of_reference_data_displX = 0.0
        reference_data_displY = []; calculation_data_displY = []; summation_of_reference_data_displY = 0.0
        reference_data_displZ = []; calculation_data_displZ = []; summation_of_reference_data_displZ = 0.0
        reference_data_temp = []; calculation_data_temp = []; summation_of_reference_data_temp = 0.0

        with open('reference_data/benchmark' + str(sys.argv[1]) + '_graph.dat') as inf:
            for line in inf:
                parts = line.split()
                reference_data_displX.append(float(parts[1]))
                reference_data_displY.append(float(parts[2]))
                reference_data_displZ.append(float(parts[3]))
                reference_data_temp.append(float(parts[4]))

        with open(self.output_filename) as inf:
            for line in inf:
                parts = line.split()
                calculation_data_displX.append(float(parts[1]))
                calculation_data_displY.append(float(parts[2]))
                calculation_data_displZ.append(float(parts[3]))
                calculation_data_temp.append(float(parts[4]))
                
        generated_data_error_displX = 0.0
        generated_data_error_displY = 0.0
        generated_data_error_displZ = 0.0
        generated_data_error_temp = 0.0

        for i in reference_data_displX:
            summation_of_reference_data_displX += fabs(i)
        for j in reference_data_displY:
            summation_of_reference_data_displY += fabs(j)
        for k in reference_data_displZ:
            summation_of_reference_data_displZ += fabs(k)
        for l in reference_data_temp:
            summation_of_reference_data_temp += fabs(l)

        for i, j in zip(calculation_data_displX, reference_data_displX):
            generated_data_error_displX += fabs(i-j)
        generated_data_error_displX /= summation_of_reference_data_displX

        for k, l in zip(calculation_data_displY, reference_data_displY):
            generated_data_error_displY += fabs(k-l)
        generated_data_error_displY /= summation_of_reference_data_displY
        
        for m, n in zip(calculation_data_displZ, reference_data_displZ):
            generated_data_error_displZ += fabs(m-n)
        generated_data_error_displZ /= summation_of_reference_data_displZ
        
        for o, p in zip(calculation_data_temp, reference_data_temp):
            generated_data_error_temp += fabs(o-p)
        generated_data_error_temp /= summation_of_reference_data_temp

        print("Error in X displacement =", 100 * generated_data_error_displX, "%")
        print("Error in Y displacement =", 100 * generated_data_error_displY, "%")
        print("Error in Z displacement =", 100 * generated_data_error_displZ, "%")
        print("Error in temperature =", 100 * generated_data_error_temp, "%")

        error1 = 100 * generated_data_error_displX
        error2 = 100 * generated_data_error_displY
        error3 = 100 * generated_data_error_displZ
        error4 = 100 * generated_data_error_temp

        return error1, error2, error3, error4

    def finalize_data(self):
        self.simulation_graph.close()

def delete_archives():

    #.......................Removing extra files
    files_to_delete_list = glob('*.time')
    files_to_delete_list.extend(glob('*.dat'))
    files_to_delete_list.extend(glob('*.gp'))
    files_to_delete_list.extend(glob('*.lst'))
    files_to_delete_list.extend(glob('*.post.res'))
    files_to_delete_list.extend(glob('*.post.msh'))
    files_to_delete_list.extend(glob('*.info'))

    for to_erase_file in files_to_delete_list:
        os.remove(to_erase_file)
