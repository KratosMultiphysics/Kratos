from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.MultilevelMonteCarloApplication as KratosMLMC

# Avoid printing of Kratos informations
KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) # avoid printing of Kratos things

# Importing the problem analysis stage class
from analysis_stage import AnalysisStage

# Import Monte Carlo library
import mc_utilities_new_analysis as mc
import cmlmc_utilities_new_analysis as mlmc


'''
This Analysis Stage implementation solves the elliptic PDE in (0,1)^2 with zero Dirichlet boundary conditions
-lapl(u) = xi*f,    f= -432*x*(x-1)*y*(y-1)
                    f= -432*(x**2+y**2-x-y)
where xi is a Beta(2,6) random variable, and computes statistic of the QoI
Q = int_(0,1)^2 u(x,y)dxdy
'''
class SimulationScenario(AnalysisStage):
    '''Main analysis stage for Monte Carlo simulations'''
    def __init__(self,input_model,input_parameters,sample):
        self.sample = sample
        super(SimulationScenario,self).__init__(input_model,input_parameters)
        self._GetSolver().main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

    def _CreateSolver(self):
        import convection_diffusion_stationary_solver
        return convection_diffusion_stationary_solver.CreateSolver(self.model,self.project_parameters["solver_settings"])

    '''Introduce here the stochasticity in the right hand side defining the forcing function and apply the stochastic contribute'''
    def ModifyInitialProperties(self):
        for node in self.model.GetModelPart("MLMCLaplacianModelPart").Nodes:
            coord_x = node.X
            coord_y = node.Y
            # forcing = -432.0 * coord_x * (coord_x - 1) * coord_y * (coord_y - 1)
            forcing = -432.0 * (coord_x**2 + coord_y**2 - coord_x - coord_y) # this forcing presents the below commented analytical solution
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX,forcing*self.sample)

    '''
    function evaluating the QoI of the problem: int_{domain} TEMPERATURE(x,y) dx dy
    midpoint rule used to compute the integral
    '''
    def EvaluateQuantityOfInterest(self):
        KratosMultiphysics.CalculateNodalAreaProcess(self._GetSolver().main_model_part,2).Execute()
        Q = 0.0
        for node in self._GetSolver().main_model_part.Nodes:
            Q = Q + (node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)*node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        return Q

####################################################################################################


if __name__ == '__main__':

    # '''set the ProjectParameters.json path'''
    # project_parameters_path = "../tests/PoissonSquareTest/parameters_poisson_finer.json"
    # '''customize setting parameters of the ML simulation'''
    # settings_MC_simulation = KratosMultiphysics.Parameters("""
    # {
    #     "tolerance" : 0.1,
    #     "cphi" : 5e-1,
    #     "batch_size" : 20,
    #     "convergence_criteria" : "MC_higher_moments_sequential_stopping_rule"
    # }
    # """)
    # '''contruct MonteCarlo or MultilevelMonteCarlo class'''
    # mc_manager = mc.MonteCarlo(settings_MC_simulation,project_parameters_path,SimulationScenario)
    # '''execute algorithm'''
    # mc_manager.Run()



    '''set the ProjectParameters.json path'''
    project_parameters_path = "../tests/PoissonSquareTest/parameters_poisson_coarse.json"
    '''customize setting parameters of the ML simulation'''
    settings_ML_simulation = KratosMultiphysics.Parameters("""
    {
        "tol0"                            : 0.25,
        "tolF"                            : 0.1,
        "cphi"                            : 1.0,
        "number_samples_screening"        : 25,
        "Lscreening"                      : 2,
        "Lmax"                            : 4,
        "initial_mesh_size"               : 0.5
    }
    """)
    '''customize setting parameters of the metric of the adaptive refinement utility'''
    custom_metric_refinement_parameters = KratosMultiphysics.Parameters("""
        {
            "hessian_strategy_parameters"           :{
                    "metric_variable"               : ["TEMPERATURE"],
                    "estimate_interpolation_error"  : false,
                    "interpolation_error"           : 0.004
            },
            "anisotropy_remeshing"                  : true,
            "anisotropy_parameters":{
                "reference_variable_name"           : "TEMPERATURE",
                "hmin_over_hmax_anisotropic_ratio"  : 0.15,
                "boundary_layer_max_distance"       : 1.0,
                "interpolation"                     : "Linear"
            },
            "local_gradient_variable"               : "TEMPERATURE"
        }
    """)
    '''customize setting parameters of the remesh of the adaptive refinement utility'''
    custom_remesh_refinement_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"                            : 0
        }
    """)
    '''contruct MultilevelMonteCarlo class'''
    mlmc_manager = mlmc.MultilevelMonteCarlo(settings_ML_simulation,project_parameters_path,custom_metric_refinement_parameters,custom_remesh_refinement_settings,SimulationScenario)
    mlmc_manager.RunScreening()
    mlmc_manager.Run()
