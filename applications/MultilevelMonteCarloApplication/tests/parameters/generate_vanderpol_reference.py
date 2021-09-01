"""Generate reference data to test VanDerPol class.

Simple script to be run without argument. It creates or overwrite a JSON file containing reference data for the test case solverWrapperTest.VanDerPolTest. It is merely a call to solverWrapperTest.VanDerPolTest.generate_reference_solution with preset parameters."""

from solverWrapperTest import VanDerPol, VanDerPolTest

vdp = VanDerPol(
    duration=20,
    timestep=0.1,
    damping=1.3,
    initial=(1.5, 0.8),
    forcingAmplitude=2.1,
    index=(0,),
    outputDimension=1,
)
file_name = "parameters/vanderpol_reference_solution.json"

VanDerPolTest.generate_reference_solution(file_name, vdp)
