from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_static_solver as BaseSolver

def CreateSolver(main_model_part, custom_settings):
    return PfemStaticMechanicalSolver(main_model_part, custom_settings)

class PfemStaticMechanicalSolver(BaseSolver.StaticMechanicalSolver):
    """The pfem solid mechanics static solver to add Pfem solid variables

    This class creates the mechanical solvers for static analysis.

    Public member variables:

    See solid_mechanics_solver.py for more information.
    """    
    def __init__(self, main_model_part, custom_settings): 
        
        super(PfemStaticMechanicalSolver, self).__init__(main_model_part, custom_settings)

    def SetVariables(self):

        super(PfemStaticMechanicalSolver, self).SetVariables()
        self.nodal_variables = self.nodal_variables + ['CONTACT_FORCE']

        if self._check_input_dof("WATER_PRESSURE"):
            # Add specific variables for the problem (pressure dofs)
            self.dof_variables = self.dof_variables + ['WATER_PRESSURE']
            self.dof_reactions = self.dof_reactions + ['REACTION_WATER_PRESSURE']
            
        if self._check_input_dof("JACOBIAN"):
            # Add specific variables for the problem (jacobian dofs)
            self.dof_variables = self.dof_variables + ['JACOBIAN']
            self.dof_reactions = self.dof_reactions + ['REACTION_JACOBIAN']


 
