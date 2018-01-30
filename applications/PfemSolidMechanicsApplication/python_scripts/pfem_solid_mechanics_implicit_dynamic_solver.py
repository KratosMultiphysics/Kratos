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

    def SetVariables(self):

        super(PfemDynamicMechanicalSolver, self).SetVariables()
        self.nodal_variables = self.nodal_variables + ['CONTACT_FORCE']

        if self._check_input_dof("WATER_DISPLACEMENT"):
            self.dof_variables = self.dof_variables + ['WATER_DISPLACEMENT','WATER_VELOCITY','WATER_ACCELERATION']
            self.dof_reactions = self.dof_reactions + ['WATER_DISPLACEMENT_REACTION','WATER_VELOCITY_REACTION','WATER_ACCELERATION_REACTION']


    def _create_solution_scheme(self):
        
        scheme_type = self.settings["scheme_type"].GetString()
        
        damp_factor_m = 0.0
        mechanical_scheme = KratosPfemSolid.ResidualBasedUWBossakScheme(damp_factor_m, 1.0)
                    
        return mechanical_scheme
 
