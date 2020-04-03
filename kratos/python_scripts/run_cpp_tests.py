import KratosMultiphysics as KM

try:
    import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
except ImportError as e:
    KM.Logger.PrintWarning("StructuralMechanicsApplication is not available")
try:
    import KratosMultiphysics.FluidDynamicsApplication as FluidDynamicsApplication
except ImportError as e:
    KM.Logger.PrintWarning("FluidDynamicsApplication is not available")
try:
    import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
except ImportError as e:
    KM.Logger.PrintWarning("ConvectionDiffusionApplication is not available")

KM.Tester.SetVerbosity(KM.Tester.Verbosity.TESTS_OUTPUTS)
KM.Tester.RunAllTestCases()
