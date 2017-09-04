from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI
try:
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    have_trilinos = True
except ImportError:
    have_trilinos = False

def CreateConvergenceAccelerator(configuration):

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object.")

    convergence_accelerator_type = configuration["solver_type"].GetString()

    if(convergence_accelerator_type == "Relaxation"):
        return KratosFSI.AitkenConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "MVQN"):
        return KratosFSI.MVQNFullJacobianConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "MVQN_recursive"):
        return KratosFSI.MVQNRecursiveJacobianConvergenceAccelerator(configuration)

    else:
        raise Exception("Convergence accelerator not found. Asking for : " + convergence_accelerator_type)

    return convergence_accelerator

def CreateTrilinosConvergenceAccelerator(configuration):

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object.")

    if not have_trilinos:
        raise Exception("Trying to create a Trilinos convergence accelerator, but TrilinosApplication could not be found.")

    convergence_accelerator_type = configuration["solver_type"].GetString()
    
    if(convergence_accelerator_type == "Relaxation"):
        return KratosTrilinos.TrilinosAitkenConvergenceAccelerator(configuration)

    else:
        raise Exception("Trilinos convergence accelerator not found. Asking for : " + convergence_accelerator_type)

    return convergence_accelerator
