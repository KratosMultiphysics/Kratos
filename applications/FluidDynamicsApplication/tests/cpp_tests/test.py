
from KratosMultiphysics import *
# If you want to run tests defined in an application, import it here.

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication

from fluid_dynamics_analysis import FluidDynamicsAnalysis

import sys
import time

Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level
# Tester.ListOfAllTestCases()

print(" Element in hydrstatic case ---- below ----- ")
Tester.RunTestCases("TestElementTwoFluidNavierStokes2D3NHydrostatic") # for 2D

print(" Element with BEHR2004 in hydrstatic case ---- below ----- ")
Tester.RunTestCases("TestElementTwoFluidNavierStokes2D3NHydrostaticBehr")

print(" Element with BEHR2004 in dynamic case ---- below ----- ")
Tester.RunTestCases("TestElementTwoFluidNavierStokes2D3NMotionBehr") 

print(" Element with BEHR2004 in dynamic case ---- below ----- ")
Tester.RunTestCases("TestElementTwoFluidNavierStokes2D3NMotionBehr2") 

print(" Element with BEHR2004 in 3D dynamic case ---- below ----- ")
Tester.RunTestCases("TestElementTwoFluidNavierStokes2D3NMotionBehr3") 
