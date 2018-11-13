import KratosSwimmingDEM as script
import os

import pre_calculated_fluid_analysis
test = script.Solution(pre_calculated_fluid_analysis)
test.Run()
