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
    return Q


@ExaquteTask(returns=1)
def ExecuteMonteCarloAnalysis_Task(pickled_model,pickled_parameters,current_analysis_stage):
    '''overwrite the old model serializer with the unpickled one'''
    model_serializer = pickle.loads(pickled_model)
    current_model = KratosMultiphysics.Model()
    model_serializer.Load("ModelSerialization",current_model)
    del(model_serializer)
    '''overwrite the old parameters serializer with the unpickled one'''
    serialized_parameters = pickle.loads(pickled_parameters)
    current_parameters = KratosMultiphysics.Parameters()
    serialized_parameters.Load("ParametersSerialization",current_parameters)
    del(serialized_parameters)
    sample = GenerateSample() # <------------------ in MLMC will be in the loop!!!!
    simulation = current_analysis_stage(current_model,current_parameters,sample)
    simulation.Run()
    QoI = EvaluateQuantityOfInterest(simulation)
    return QoI


if __name__ == '__main__':

    '''customize setting parameters of the ML simulation'''
    settings_MC_simulation = KratosMultiphysics.Parameters("""
    {
    }
    """)
    '''contruct MonteCarlo class'''
    mc_class = mc.MonteCarlo(settings_MC_simulation)
    mc_class.SetAnalysis(MonteCarloAnalysis)
    SERIALIZE
    while mc_class.convergence is not True:
        mc_class.InitializeMCPhase()
        mc_class.ScreeningInfoInitializeMCPhase()
        for instance in range (mc_class.difference_number_samples[0]):
            mc_class.AddResults(ExecuteMonteCarloAnalysis_Task(pickled_model,pickled_parameters,mc_class.GetAnalysis()))
        mc_class.FinalizeMCPhase()
        mc_class.ScreeningInfoFinalizeMCPhase()
