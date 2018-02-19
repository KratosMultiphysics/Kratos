from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def GetSolvers(SolversDataList):
    SolversList = {}
    num_solvers = SolversDataList.size()

    for i in range(0,num_solvers):
        solver_module = __import__(SolversDataList[i]['type'].GetString())
        solver_name = SolversDataList[i]['name'].GetString()
        solver = solver_module.Create(SolversDataList[i])
        SolversList[solver_name] = solver
        
    return SolversList


def Extract(d, keys):
    return dict((k, d[k]) for k in keys if k in d)


def GetSolverCoSimulationDetails(participants):
    num_participants = participants.size()
    solver_cosim_details = {}
    for i in range(num_participants):
        participant = participants[i]
        solver_cosim_details[participant['name'].GetString()] = participant

    return solver_cosim_details

    