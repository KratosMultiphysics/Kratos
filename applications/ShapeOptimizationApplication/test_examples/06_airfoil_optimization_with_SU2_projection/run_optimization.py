
# Import Kratos core and apps
import KratosMultiphysics as KM

# Additional imports
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.interface_su2 import InterfaceSU2
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass

# =======================================================================================================
# Define external analyzer
# =======================================================================================================

with open("parameters.json",'r') as parameter_file:
    parameters = KM.Parameters(parameter_file.read())

interface_su2 = InterfaceSU2(parameters["su2_interface_settings"])
# interface_su2.WriteSU2MeshAsMDPA()

# Definition of external analyzer
class CustomSU2Analyzer(AnalyzerBaseClass):

    # --------------------------------------------------------------------------
    def __init__( self ):
        interface_su2.InitializeNewSU2Project()

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        if optimization_iteration == 1:
            interface_su2.WriteNodesAsSU2MeshMotionFile(current_design.GetNodes())
        else:
            previos_iteration = optimization_iteration-1
            interface_su2.WriteNodesAsSU2MeshMotionFile(current_design.GetNodes(),"DESIGNS/DSN_"+str(previos_iteration).zfill(3))

        if communicator.isRequestingValueOf("drag") or communicator.isRequestingValueOf("lift"):
            update_mesh = True
            [drag,lift] = interface_su2.ComputeValues(["DRAG","LIFT"], update_mesh, optimization_iteration)

            if communicator.isRequestingValueOf("drag"):
                communicator.reportValue("drag", drag)

            if communicator.isRequestingValueOf("lift"):
                communicator.reportValue("lift", lift)

        if communicator.isRequestingGradientOf("drag"):
            update_mesh = False
            [drag_gradient] = interface_su2.ComputeGradient(["DRAG"], update_mesh, optimization_iteration)
            communicator.reportGradient("drag", drag_gradient)

        if communicator.isRequestingGradientOf("lift"):
            update_mesh = False
            [lift_gradient] = interface_su2.ComputeGradient(["LIFT"], update_mesh, optimization_iteration)
            communicator.reportGradient("lift", lift_gradient)

# =======================================================================================================
# Perform optimization
# =======================================================================================================
model = KM.Model()

# Create optimizer and perform optimization
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model, CustomSU2Analyzer())
optimizer.Optimize()