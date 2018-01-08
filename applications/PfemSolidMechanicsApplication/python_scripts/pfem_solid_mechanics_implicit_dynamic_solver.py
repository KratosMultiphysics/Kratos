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
        
        integration_method = self.time_integration_settings["integration_method"].GetString()

        #if(integration_method == "Newmark"):           
        damp_factor_m = 0.0
        mechanical_scheme = KratosPfemSolid.ResidualBasedUWBossakScheme(damp_factor_m, 1.0)
                    
        return mechanical_scheme
 
