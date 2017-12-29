from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import *
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
import os

parameter_file = open("ProjectParameters.json", 'r')
ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

newStrategyModule = __import__(ProjectParameters[
    "co_simulation_solver_settings"]["type"].GetString())
solver = newStrategyModule.CreateSolver(
    ProjectParameters["co_simulation_solver_settings"])
solver.Initialize()

delta_time = ProjectParameters["problem_data"]["time_step"].GetDouble()
# Start time
time = ProjectParameters["problem_data"]["start_time"].GetDouble()
# End time
end_time = ProjectParameters["problem_data"]["end_time"].GetDouble()

step = 0
time = 0
# Solving the problem (time integration)
while (time <= end_time):
    time = time + delta_time
    step = step + 1

    solver.InitializeSolutionStep()
    print('############## ')
    print('Step :: ', step)
    print('Time :: ', time)

    solver.Solve()

    solver.FinalizeSolutionStep()