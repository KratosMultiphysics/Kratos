# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.OptimizationApplication import optimizer
from KratosMultiphysics.KratosUnittest import TestCase
import csv, os, shutil

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

# Defining the model_part
model = KM.Model()

# =======================================================================================================
# Create optimizer & Perform optimization
# =======================================================================================================

optimizer = optimizer.CreateOptimizer(parameters["optimization_settings"],model)
optimizer.model_parts_controller.Initialize()
optimizer.analyses_controller.Initialize()
optimizer.responses_controller.Initialize()
optimizer.controls_controller.Initialize()
optimizer.optimizations_controller.Initialize()
optimizer.optimizations_controller.OptimizeAll()

# =======================================================================================================
# Test results and clean directory
# =======================================================================================================
optimization_log_filename = "optimization_log.csv"

# Testing
original_directory = os.getcwd()

with open(optimization_log_filename, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    last_line = None
    for line in reader:
        if not line:
            continue
        else:
            last_line = line

    resulting_mass = float(last_line[1].strip())
    resulting_strain_energy = float(last_line[5].strip())
    resulting_sina_alpha = float(last_line[12].strip())

    # # Check against specifications
    TestCase().assertAlmostEqual(resulting_mass, 3736.87,4)
    TestCase().assertAlmostEqual(resulting_strain_energy, 0.00529714, 6)
    TestCase().assertAlmostEqual(resulting_sina_alpha, 0.87266, 5)

shutil.rmtree(original_directory+"/Optimization_Results")
shutil.rmtree(original_directory+"/Primal_Results")
os.remove(optimization_log_filename)

os.chdir(original_directory)
