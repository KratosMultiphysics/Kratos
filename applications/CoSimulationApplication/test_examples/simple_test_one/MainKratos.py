from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics import *
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
import tools

parameter_file = open("ProjectParameters.json", 'r')
ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

participating_solvers = tools.GetSolvers(ProjectParameters["solvers"])

new_strategy_module = __import__(ProjectParameters["co_simulation_solver_settings"]["type"].GetString())
solver = new_strategy_module.CreateSolver(ProjectParameters)
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

    solver.InitializeTimeStep()
    print('############## ')
    print('Step :: ', step)
    print('Time :: ', time)

    solver.SolveTimeStep()
    
    solver.FinalizeTimeStep()

