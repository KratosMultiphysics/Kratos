##with refinement
def ExecuteMultilevelAnalisys(lev_fine, data_coarse=[ModelCoarse, ParametersCoarse] ):
    #generate sample
    sample_parameters = data_fine["sample_parameters"]
    sample_generator.Execute(sample_parameters)

    if(lev_fine == 0):
        ...
    else:
        for i in range(lev_fine):
            #solve coarse domain

            #estimate error and refine mesh

            #solve fine domain

            #solve fine domain

            coarse_data = fine_data ... go on with the loop

        return(    (lev_fine, lev_coarse), QoIfine, QoIcoarse, total_time  )


##with precomputed meshes
def ExecuteMultilevelAnalisys(lev_fine, data_coarse=[ModelCoarse, ParametersCoarse], data_fine ):
    #generate sample
    sample_parameters = data_fine["sample_parameters"]
    sample_generator.Execute(sample_parameters)

    if(lev_fine == 0):
        ...
    else:
        results_coarse = ExecuteAnalysis(lev_fine-1, data[lev_fine-1])
        results_fine = ExecuteAnalysis(lev_fine, data[lev_fine])
 
        return(    (lev_fine, lev_coarse), QoIfine, QoIcoarse, total_time  )




MLMC = MLMC_Driver(5)
MLMC["number_of_samples"][0] = 1000

#screening phase
for level in range(MLMC["number_of_levels"]+1):
    for instance in range (MLMC["number_of_samples"]):
        MLMC.AddResult( ExecuteMultilevelAnalysis(level, ) )

MLMC.FinalizeScreening() ##compute MLMC coefficients, variance and number of iterations

MLMC.PrintScreeningResults()

###MLMC
while(not converged):
    for level in range(MLMC["number_of_levels"]+1):
        for instance in range (MLMC["number_of_samples"]):
            MLMC.AddResult( ExecuteMultilevelAnalysis(level, MLMC["simulation_parameters"][level]) )
    MLMC.ComputeContinuationParameters()
    converged = MLMC.IsConverged()



class MLMC_Driver
    def __init__(max_levels):
        ...

    

    def Screening(...)

    def EstimateNumber

    
    