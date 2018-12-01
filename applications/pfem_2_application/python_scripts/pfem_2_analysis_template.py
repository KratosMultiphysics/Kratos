from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosPFEM2Application as Pfem2
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

def InstantiatePFEM2izedAnalysis(base_class,model,parameters):
    class PFEM2Analysis(base_class):
        def __init__(self,model,parameters):
            super(PFEM2Analysis, self).__init__(model,parameters)

        def _GetSimulationName(self):
            return "PFEM2 Analysis"

    return PFEM2Analysis(model,parameters)

if __name__ == '__main__':
    from sys import argv

    if len(argv) == 2: # Base class file name is being passed from outside
        base_class_name = argv[1]
    else: # using default name
        base_class_name = "fluid_dynamics"

    if base_class_name == "fluid_dynamics":
        import KratosMultiphysics.FluidDynamicsApplication
        from fluid_dynamics_analysis import FluidDynamicsAnalysis
        base_class = FluidDynamicsAnalysis
    elif base_class_name == "shallow_water":
        import KratosMultiphysics.ShallowWaterApplication
        from shallow_water_analysis import ShallowWaterAnalysis
        base_class = ShallowWaterAnalysis
    elif base_class_name == "fluid_transport":
        import KratosMultiphysics.ConvectionDiffusionApplication
        import KratosMultiphysics.FluidDynamicsApplication
        import KratosMultiphysics.FluidTransportApplication
        from fluid_transport_analysis import FluidTransportAnalysis
        base_class = FluidTransportAnalysis
    else:
        raise Exception("Trying to derive PFEM2Analysis from an unknown class name!")

    model = Kratos.Model()

    parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    simulation = InstantiatePFEM2izedAnalysis(base_class,model,parameters)
    simulation.Run()

