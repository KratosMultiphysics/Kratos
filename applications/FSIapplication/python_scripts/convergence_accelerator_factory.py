from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI
try:
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    have_trilinos = True
except ImportError:
    have_trilinos = False

_aitken_defaults = KratosMultiphysics.Parameters("""{   "solver_type"       : "Relaxation",
                                                        "acceleration_type" : "Aitken",
                                                        "w_0"               : 0.825  }""")

_mvqn_defaults = KratosMultiphysics.Parameters("""{ "solver_type" : "MVQN",
                                                    "w_0"         : 0.825   }""")

_mvqn_defaults_recursive = KratosMultiphysics.Parameters("""{ "solver_type" : "MVQN_recursive",
                                                              "w_0"         : 0.825,
                                                              "buffer_size" : 7           }""")

def CreateConvergenceAccelerator(configuration):

    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object.")

    convergence_accelerator_type = configuration["solver_type"].GetString()

    if(convergence_accelerator_type == "Relaxation"):

        configuration.ValidateAndAssignDefaults(_aitken_defaults)

        if (configuration["acceleration_type"].GetString() == "Aitken"):
            convergence_accelerator = KratosFSI.AitkenConvergenceAccelerator(configuration["w_0"].GetDouble())

    elif(convergence_accelerator_type == "MVQN"):

        configuration.ValidateAndAssignDefaults(_mvqn_defaults)

        convergence_accelerator = KratosFSI.MVQNFullJacobianConvergenceAccelerator(configuration["w_0"].GetDouble())

    elif(convergence_accelerator_type == "MVQN_recursive"):

        configuration.ValidateAndAssignDefaults(_mvqn_defaults_recursive)

        convergence_accelerator = KratosFSI.MVQNRecursiveJacobianConvergenceAccelerator(configuration["w_0"].GetDouble(),
                                                                                        configuration["buffer_size"].GetInt())

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

        configuration.ValidateAndAssignDefaults(_aitken_defaults)

        if (configuration["acceleration_type"].GetString() == "Aitken"):
            convergence_accelerator = KratosTrilinos.TrilinosAitkenConvergenceAccelerator(configuration["w_0"].GetDouble())

    elif(convergence_accelerator_type == "MVQN_recursive"):

        configuration.ValidateAndAssignDefaults(_mvqn_defaults_recursive)

        convergence_accelerator = KratosTrilinos.TrilinosMVQNRecursiveJacobianConvergenceAccelerator(interface_model_part,
                                                                                                     epetra_communicator,
                                                                                                     configuration["w_0"].GetDouble(),
                                                                                                     configuration["buffer_size"].GetInt())

    else:
        raise Exception("Trilinos convergence accelerator not found. Asking for : " + convergence_accelerator_type)

    return convergence_accelerator
