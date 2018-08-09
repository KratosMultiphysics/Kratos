# This script couples DEM and Nvidia Flex...

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import sys

# Kratos
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.NvidiaFlexApplication import *

import main_script

class Solution(main_script.Solution):

    def Run(self):

        self.nvidia_flex_wrapper = FlexWrapper(self.spheres_model_part, self.creator_destructor)
        super(Solution, self).Run()

    def SolverSolve(self):
        if self.step < 2:
        #    for node in self.spheres_model_part.Nodes:
        #        node.SetSolutionStepValue(VELOCITY_X, 0.0)
        #        node.SetSolutionStepValue(VELOCITY_Y, 0.0)
        #        node.SetSolutionStepValue(VELOCITY_Z, 0.0)
            self.nvidia_flex_wrapper.UpdateFlex()

        self.nvidia_flex_wrapper.SolveTimeSteps(100*self.dt, 100)
        self.nvidia_flex_wrapper.TransferDataFromFlexToKratos()


if __name__=="__main__":
    Solution().Run()