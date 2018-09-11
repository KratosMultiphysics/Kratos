from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import json
import KratosMultiphysics
from KratosMultiphysics import *
import KratosMultiphysics.CoSimulationApplication as CoSimulationApplication
import custom_convergence_criteria as ccc
from custom_co_simulation_solver_interfaces.co_simulation_base_solver import CoSimulationBaseSolver

from custom_co_simulation_coupled_solvers import *
import custom_co_simulation_coupled_solvers as cc
import custom_co_simulation_coupled_solvers.co_simulation_coupled_solver_factory as fact

filename = '/home/aditya/work/Kratos/applications/CoSimulationApplication/test_examples/test_one/test_json.json'

#Read JSON data into the datastore variable
if filename:
    with open(filename, 'r') as f:
        json_data = json.load(f)


gs = fact.CreateCoupledSolver(json_data)