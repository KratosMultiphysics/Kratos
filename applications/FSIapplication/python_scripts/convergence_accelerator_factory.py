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

def CreateTrilinosConvergenceAccelerator(interface_model_part, epetra_communicator, configuration):

    if not have_trilinos:
        raise Exception("Trying to create a Trilinos convergence accelerator, but TrilinosApplication could not be found.")

    if (type(interface_model_part) != KratosMultiphysics.ModelPart):
        raise Exception("First input in Trilinos convergence accelerator factory is expceted to be provided as a Kratos ModelPart object.")

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("Third input in Trilinos convergence accelerator factory is expected to be provided as a Kratos Parameters object.")

    convergence_accelerator_type = configuration["solver_type"].GetString()
    
    if(convergence_accelerator_type == "Relaxation"):
        return KratosTrilinos.TrilinosAitkenConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "MVQN_recursive"):
        return KratosTrilinos.TrilinosMVQNRecursiveJacobianConvergenceAccelerator(interface_model_part, epetra_communicator, configuration)

    else:
        raise Exception("Trilinos convergence accelerator not found. Asking for : " + convergence_accelerator_type)

    return convergence_accelerator
