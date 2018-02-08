from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_implicit_dynamic_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return PfemDynamicMechanicalSolver(main_model_part, custom_settings)

class PfemDynamicMechanicalSolver(BaseSolver.ImplicitMechanicalSolver):
    """The pfem solid mechanics dynamic solver to add Pfem solid variables

    This class creates the mechanical solvers for dynamic analysis.

    Public member variables:

    See solid_mechanics_solver.py for more information.
    """    
    def __init__(self, main_model_part, custom_settings): 
        
        super(PfemDynamicMechanicalSolver, self).__init__(main_model_part, custom_settings)



    def _create_solution_scheme(self):
        
        #scheme_type = self.settings["solver_type"].GetString()
	#integration_method = self.time_integration_settings["integration_method"].GetString()
        
#        beta = self.settings["beta"].GetDouble()
        damp_factor_m = 0.0
#        alphaM =-0.1;
        alphaM = 0.0
        dynamic = 1
        alphaF = 0
        beta = 0.3025
        gamma = 0.6
#        beta = 0.25*(1-alphaM+alphaF)*(1-alphaM+alphaF);
#        gamma = 0.5-alphaM+alphaF;
        mechanical_scheme = KratosPfemSolid.ResidualBasedUWBossakScheme(alphaM, dynamic, alphaF, beta, gamma)

        self.process_info[KratosMultiphysics.NEWMARK_BETA] = beta
        self.process_info[KratosMultiphysics.NEWMARK_GAMMA] = gamma
        self.process_info[KratosMultiphysics.BOSSAK_ALPHA] = alphaM

        time_integration_method = KratosSolid.BossakMethod()
        time_integration_method.AddToProcessInfo(KratosSolid.TIME_INTEGRATION_METHOD, time_integration_method, self.process_info)
        time_integration_method.SetParameters(self.process_info)

        mechanical_scheme = KratosPfemSolid.NewResidualBasedUWBossakScheme()
                    
        return mechanical_scheme
 
