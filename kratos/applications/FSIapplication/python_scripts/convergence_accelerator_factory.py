from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI

def CreateConvergenceAccelerator(configuration):
    
    if(type(configuration) != KratosMultiphysics.Parameters):
        raise Exception("Input is expected to be provided as a Kratos Parameters object.")
    
    convergence_accelerator_type = configuration["solver_type"].GetString()
    
    if(convergence_accelerator_type == "Relaxation"):
        
        defaults_relaxation = KratosMultiphysics.Parameters("""{
                                                                "solver_type"       : "Relaxation",
                                                                "acceleration_type" : "Aitken",
                                                                "w_0"               : 0.825
                                                               }""")
                                                               
        configuration.ValidateAndAssignDefaults(defaults_relaxation)        
        
        if (configuration["acceleration_type"].GetString() == "Aitken"):
            convergence_accelerator = KratosFSI.AitkenConvergenceAccelerator(configuration["w_0"].GetDouble())       
            
    elif(convergence_accelerator_type == "MVQN"):
        
        defaults_MVQN = KratosMultiphysics.Parameters("""{
                                                          "solver_type" : "MVQN",
                                                          "w_0"         : 0.825
                                                         }""")
                                                               
        configuration.ValidateAndAssignDefaults(defaults_MVQN)
        
        convergence_accelerator = KratosFSI.MVQNFullJacobianConvergenceAccelerator(configuration["w_0"].GetDouble())
            
    elif(convergence_accelerator_type == "MVQN_recursive"):
        
        defaults_MVQN_recursive = KratosMultiphysics.Parameters("""{
                                                                    "solver_type" : "MVQN_recursive",
                                                                    "w_0"         : 0.825,
                                                                    "buffer_size" : 7
                                                                   }""")
                                                               
        configuration.ValidateAndAssignDefaults(defaults_MVQN_recursive)
        
        convergence_accelerator = KratosFSI.MVQNRecursiveJacobianConvergenceAccelerator(configuration["w_0"].GetDouble(),
                                                                                        configuration["buffer_size"].GetInt())
        
    else:
        raise Exception("Convergence accelerator not found. Asking for : " + convergence_accelerator_type)

    return convergence_accelerator
