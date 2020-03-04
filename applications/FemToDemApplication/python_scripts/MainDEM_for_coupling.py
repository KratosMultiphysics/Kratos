from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics.DEMApplication
import KratosMultiphysics.DEMApplication.main_script as MainDEM
import KratosMultiphysics.FemToDemApplication as FEMDEM

class DEM_for_coupling_Solution(MainDEM.Solution):

    def Info(self):
        print("DEM part of the FEM-DEM application")

    def SetAnalyticParticleWatcher(self):
        pass

    def AddVariables(self):
        super(DEM_for_coupling_Solution).AddVariables()
        # For averaging forces when substepping
        spheres_model_part.AddNodalSolutionStepVariable(CONTACT_IMPULSE)
