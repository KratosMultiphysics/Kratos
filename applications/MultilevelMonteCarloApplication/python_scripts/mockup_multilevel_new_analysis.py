class PotentialFlowMonteCarlo(PotentialFlow):
    '''Main analysis stage for Monte Carlo simulations'''
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(PotentialFlowMonteCarlo,self).__init__(input_model,input_parameters)

    '''Introduce here the stochasticity in the right hand side defining the forcing function and apply the stochastic contribute'''
    def ModifyInitialProperties(self):
        do things using self.sample, e.g.
        for node in self.model.GetModelPart("MLMCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y) # this forcing presents an analytical solution
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)



def GenerateSample():
    alpha = 2.0
    beta = 6.0
    number_samples = 1
    sample = np.random.beta(alpha,beta,number_samples)
    return sample

def EvaluateQuantityOfInterest(simulation):
    return QoI

def ExecuteMultilevelMonteCarloAnalisys(current_MLMC_level,
                                        pickled_coarse_model,
                                        pickled_coarse_parameters,
                                        size_meshes,
                                        pickled_settings_metric_refinement,
                                        pickled_settings_remesh_refinement,
                                        current_analysis_stage):
    sample = GenerateSample() # <--------------- SAMPLE GENERATION
    mlmc_results_class = mlmc.MultilevelMonteCarloResults()
    if (current_MLMC_level == 0):
        mlmc_results_class,pickled_current_model,pickled_current_parameters = ExecuteMultilevelMonteCarloAnalisys_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement,sample,current_MLMC_level,mlmc_results_class,current_analysis_stage)
    else:
        for level in range(current_MLMC_level+1):
            mlmc_results_class,pickled_current_model,pickled_current_parameters = ExecuteMultilevelMonteCarloAnalisys_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement,sample,level,mlmc_results_class,current_analysis_stage)
    return mlmc_results_class,current_MLMC_level

@ExaquteTask(returns=3)
def ExecuteMultilevelMonteCarloAnalisys_Task(current_MLMC_level,pickled_coarse_model,pickled_coarse_parameters,size_meshes,pickled_settings_metric_refinement,pickled_settings_remesh_refinement,sample,current_level,mlmc_results_class,current_analysis_stage):
    UNSERIALIZE
    '''refine if current current_level > 0, adaptive refinement based on the solution of previous level'''
    if (current_level > 0):
        UNSERIALIZE
        REFINEMENT
    simulation = current_analysis_stage(current_model,current_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    SERIALIZE
    UPDATE MLMC RESULT CLASS
    return mlmc_results_class,pickled_finer_model,pickled_finer_parameters

@ExaquteTask(parameter_file_name=FILE_IN,returns=2)
def SerializeModelParameters_Task(parameter_file_name):
    SERIALIZE
    return pickled_model,pickled_parameters

#TODO: to execute this function in a compss task, metric_refinement_parameters and remeshing_refinement_parameters should be read from a file
def SerializeRefinementParameters(metric_refinement_parameters,remeshing_refinement_parameters):
    SERIALIZE
    return pickled_metric_refinement_parameters,pickled_remesh_refinement_parameters


if __name__ == '__main__':

    SERIALIZE PARAMETERS, MODEL, ETC
    '''contruct MultilevelMonteCarlo class'''
    mlmc_class = mlmc.MultilevelMonteCarlo(settings_user_defined)
    mlmc.analysis = PotentialFlowMonteCarlo
    for level in range(mlmc_class.current_number_levels+1):
        for instance in range (mlmc_class.number_samples[level]):
            mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys(level,pickled_model,pickled_parameters,mlmc_class.sizes_mesh,pickled_settings_metric_refinement,pickled_settings_remesh_refinement))

    # here it is equal
    # '''finalize screening phase'''
    # mlmc_class.FinalizeScreeningPhase()
    # mlmc_class.ScreeningInfoScreeningPhase()
    # '''start MLMC phase'''
    # while mlmc_class.convergence is not True:
    #     '''initialize MLMC phase'''
    #     mlmc_class.InitializeMLMCPhase()
    #     mlmc_class.ScreeningInfoInitializeMLMCPhase()
    #     '''MLMC execution phase'''
    #     for level in range (mlmc_class.current_number_levels+1):
    #         for instance in range (mlmc_class.difference_number_samples[level]):
    #             mlmc_class.AddResults(ExecuteMultilevelMonteCarloAnalisys(level,pickled_model,pickled_parameters,mlmc_class.sizes_mesh,pickled_settings_metric_refinement,pickled_settings_remesh_refinement))
    #     '''finalize MLMC phase'''
    #     mlmc_class.FinalizeMLMCPhase()
    #     mlmc_class.ScreeningInfoFinalizeMLMCPhase()
