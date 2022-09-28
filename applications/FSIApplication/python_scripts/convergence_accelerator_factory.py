from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FSIApplication as KratosFSI

import KratosMultiphysics.kratos_utilities as KratosUtilities
have_trilinos = KratosUtilities.CheckIfApplicationsAvailable("TrilinosApplication")
is_distributed = KratosMultiphysics.DataCommunicator.GetDefault().IsDistributed()
if have_trilinos and is_distributed:
    try:
        import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    except Exception as e:
        raise Exception("Trying to create a Trilinos convergence accelerator, but TrilinosApplication could not be found.")

have_linear_solvers = KratosUtilities.CheckIfApplicationsAvailable("LinearSolversApplication")
if have_linear_solvers:
    import KratosMultiphysics.LinearSolversApplication as KratosLinearSolvers

def CreateConvergenceAccelerator(configuration):

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object.")

    convergence_accelerator_type = configuration["solver_type"].GetString()

    if(convergence_accelerator_type == "constant_relaxation"):
        return KratosFSI.ConstantRelaxationConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "Relaxation"):
        return KratosFSI.AitkenConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "MVQN"):
        return KratosFSI.MVQNFullJacobianConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "MVQN_randomized_SVD"):
        if not have_linear_solvers:
            err_msg = "MVQN with randomized SVD Jacobian requires the \'LinearSolversApplication\'."
            raise Exception(err_msg)
        bdc_svd = KratosLinearSolvers.EigenDenseBDCSVD()
        col_piv_qr = KratosLinearSolvers.EigenDenseColumnPivotingHouseholderQRDecomposition()
        return KratosFSI.MVQNRandomizedSVDConvergenceAccelerator(col_piv_qr, bdc_svd, configuration)

    elif(convergence_accelerator_type == "MVQN_recursive"):
        return KratosFSI.MVQNRecursiveJacobianConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "IBQN_MVQN"):
        return KratosFSI.IBQNMVQNConvergenceAccelerator(configuration)

    elif(convergence_accelerator_type == "IBQN_MVQN_randomized_SVD"):
        if not have_linear_solvers:
            err_msg = "MVQN with randomized SVD Jacobian requires the \'LinearSolversApplication\'."
            raise Exception(err_msg)
        bdc_svd = KratosLinearSolvers.EigenDenseBDCSVD()
        col_piv_qr = KratosLinearSolvers.EigenDenseColumnPivotingHouseholderQRDecomposition()
        return KratosFSI.IBQNMVQNRandomizedSVDConvergenceAccelerator(col_piv_qr, bdc_svd, configuration)

    else:
        raise Exception("Convergence accelerator not found. Asking for : " + convergence_accelerator_type)

    return convergence_accelerator

def CreateTrilinosConvergenceAccelerator(interface_model_part, epetra_communicator, configuration):

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
