from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.PFEM2Application as PFEM2

from KratosMultiphysics.PFEM2Application.pfem_2_base_solver import PFEM2BaseSolver

def CreateSolver(model, custom_settings):
    return PFEM2MonolithicSolver(model, custom_settings)

class PFEM2MonolithicSolver(PFEM2BaseSolver):
    def Initialize(self):
        super(PFEM2MonolithicSolver, self).Initialize()
        # Pressure projection calculation
        model_part = self.GetComputingModelPart()
        self.explicit_strategy = PFEM2.PFEM2_Explicit_Strategy(model_part,self.domain_size, self.settings["move_mesh_flag"].GetBool())

    def FinalizeSolutionStep(self):
        self.get_mesh_strategy().FinalizeSolutionStep()
        self._CalculatePressureProjection()
        self.get_particles_stage().ExecuteFinalizeSolutionStep()

    def _CalculatePressureProjection(self):
        self.GetComputingModelPart().ProcessInfo.SetValue(KM.FRACTIONAL_STEP, 10)
        self.explicit_strategy.InitializeSolutionStep()
        self.explicit_strategy.AssembleLoop()
        self.explicit_strategy.FinalizeSolutionStep()

    def _create_time_scheme(self):
        time_scheme = KM.ResidualBasedIncrementalUpdateStaticScheme()
        return time_scheme
