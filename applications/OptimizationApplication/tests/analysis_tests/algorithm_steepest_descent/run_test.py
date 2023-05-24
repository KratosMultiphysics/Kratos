import KratosMultiphysics as Kratos
from KratosMultiphysics.KratosUnittest import TestCase
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
import csv, os, shutil

if __name__ == "__main__":
    # os.chdir("/mnt/suneth/data/PostDoc/6_Simulation_Data/14_OptApp/2_reza_app_comparison/1_mass_opt/2_new_proposal")
    with open("optimization_parameters.json", "r") as file_input:
        parameters = Kratos.Parameters(file_input.read())

    model = Kratos.Model()
    analysis = OptimizationAnalysis(model, parameters)
    analysis.Run()

    # # Check against specifications
    mass_value = []
    start_found = False
    count = 0
    with open("summary.csv", 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for line in reader:
            if "#  STEP" in line:
                start_found = True
                continue
            if start_found and count < 11:
                mass_value.append(float(line[1].strip()))
                count += 1

        TestCase().assertEqual(mass_value[0], 1.386000000e+04)
        TestCase().assertEqual(mass_value[2], 1.266000000e+04)
        TestCase().assertEqual(mass_value[5], 1.086000000e+04)
        TestCase().assertEqual(mass_value[7], 9.660000000e+03)
        TestCase().assertEqual(mass_value[10], 7.860000000e+03)

