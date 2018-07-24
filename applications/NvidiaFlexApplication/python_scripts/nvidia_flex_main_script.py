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
        
        self.nvidia_wrapper = FlexWrapper()
        self.nvidia_wrapper.FlexWrapperFunction()
        super(Solution, self).Run()

if __name__=="__main__":
    Solution().Run()