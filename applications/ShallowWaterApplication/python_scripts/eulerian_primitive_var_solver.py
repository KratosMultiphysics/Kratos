from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

## Import base class file
import shallow_water_base_solver

def CreateSolver(model_part, custom_settings):
    return EulerianPrimitiveVarSolver(model_part, custom_settings)

class EulerianPrimitiveVarSolver(shallow_water_base_solver.ShallowWaterBaseSolver):
    def Solve(self):
        (self.solver).Solve()
