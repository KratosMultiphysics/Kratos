# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from KratosMultiphysics.KratosUnittest import TestCase
from algorithm_trust_region import Projector
from custom_math import PerformBisectioning
import csv

# =======================================================================================================
# Set and read input data
# =======================================================================================================

algorithm_settings = Parameters("""
{
    "name"                          : "trust_region",
    "max_step_length"               : 0.1,
    "max_iterations"                : 10,
    "far_away_length"               : 2.0,
    "subopt_max_itr"                : 50,
    "subopt_tolerance"              : 1e-10,
    "bisectioning_max_itr"          : 30,
    "bisectioning_tolerance"        : 1e-8,
    "obj_share_during_correction"   : 1
}""")

# lens and dirs are obtained from the initial iteration of a constrained optimization of a shell
len_obj = 3.1113059839426676
len_eqs = [9.116181801652377]
len_ineqs = [0.0, 12.354599999999998]

dir_obj = []
dir_eqs = [[]]
dir_ineqs = [[],[]]

with open('dirs.csv', 'r') as input_file:
    reader = csv.reader(input_file)
    for line in reader:
        dir_obj.append(float(line[0].strip()))
        dir_eqs[0].append(float(line[1].strip()))
        dir_ineqs[0].append(float(line[2].strip()))
        dir_ineqs[1].append(float(line[3].strip()))

# =======================================================================================================
# Perform test
# =======================================================================================================

projector = Projector(len_obj, dir_obj, len_eqs, dir_eqs, len_ineqs, dir_ineqs, algorithm_settings)

# Test projection
len_obj_test = 0.01
inactive_threshold = 100
test_norm_dX, is_projection_sucessfull = projector.RunProjection(len_obj_test, inactive_threshold)

# Test bisectioning
len_obj = algorithm_settings["obj_share_during_correction"].GetDouble()
func = lambda threshold: projector.RunProjection(len_obj, threshold)

threshold_min = 0
threshold_max = 1.3
bi_target = 1
bi_tolerance = algorithm_settings["bisectioning_tolerance"].GetDouble()
bi_max_itr = algorithm_settings["bisectioning_max_itr"].GetInt()
l_threshold_result, bi_itrs, bi_err = PerformBisectioning(func, threshold_min, threshold_max, bi_target, bi_tolerance, bi_max_itr)

# =======================================================================================================
# Check results
# =======================================================================================================

TestCase().assertAlmostEqual(test_norm_dX, 19.961965107754867, 6)
TestCase().assertTrue(is_projection_sucessfull)

TestCase().assertAlmostEqual(l_threshold_result, 0.5578154608607293, 6)
TestCase().assertAlmostEqual(bi_err, 3.9182936895088005e-09, 8)
TestCase().assertEqual(bi_itrs, 25)

# =======================================================================================================