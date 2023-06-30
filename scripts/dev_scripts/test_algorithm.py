import os
import shutil
from glob import glob

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis


def pre_operations():
    try:
        shutil.rmtree("./Optimization_Results")
    except:
        pass
    try:
        shutil.rmtree("./other_results")
    except:
        pass


def post_operations():
    os.makedirs("./other_results")

    # remove primal analysis files
    for file in glob('./Structure*.h5'):
        os.remove(file)

    for file in glob('./*.csv')+glob('./*.time')+glob('./*.html'):
        os.rename(str(file), "./other_results"+str(file)[1:])


working_folder = "measurement_residual_test"
os.chdir(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), working_folder)))

pre_operations()

with open("optimization_parameters.json", "r") as file_input:
    parameters = Kratos.Parameters(file_input.read())
model = Kratos.Model()
analysis = OptimizationAnalysis(model, parameters)
analysis.Run()

post_operations()
