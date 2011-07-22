from KratosTrilinosApplication import *

def MultilevelLinearSolver(tolerance, max_iterations):
    #settings for the iterative solver
    aztec_parameters = ParameterList()
    aztec_parameters.set("AZ_solver","AZ_bicgstab");
##    aztec_parameters.set("AZ_output","AZ_none");
    aztec_parameters.set("AZ_output",10);

    #settings of the ML solver
    MLList = ParameterList()
    
    default_settings = EpetraDefaultSetter()
    default_settings.SetDefaults(MLList,"SA");
    
    MLList.set("ML output", 10);
    MLList.set("max levels",3);
    MLList.set("increasing or decreasing","increasing");
    MLList.set("aggregation: type", "MIS");
##    MLList.set("coarse: type","Amesos-Superludist");
    MLList.set("smoother: type","Chebyshev");
    MLList.set("smoother: sweeps",3);
    MLList.set("smoother: pre or post", "both");
    MLList.set("ML output",0);

    #tolerance = 1e-4
    nit_max = 1000

    linear_solver =  MultiLevelSolver(aztec_parameters,MLList,tolerance,nit_max);

    return linear_solver

    
