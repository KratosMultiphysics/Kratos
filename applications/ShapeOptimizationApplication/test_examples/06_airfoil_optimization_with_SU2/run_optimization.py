# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from interface_su2 import InterfaceSU2
from analyzer_base import AnalyzerBaseClass

# =======================================================================================================
# Define external analyzer
# =======================================================================================================

with open("parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

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

        if communicator.isRequestingValueOf("drag"):
            update_mesh = True
            [value] = interface_su2.ComputeValues(["DRAG"], update_mesh, optimization_iteration)
            communicator.reportValue("drag", value)

        if communicator.isRequestingGradientOf("drag"):
            update_mesh = False
            gradient = interface_su2.ComputeGradient(["DRAG"], update_mesh, optimization_iteration)
            communicator.reportGradient("drag", gradient)


# =======================================================================================================
# Perform optimization
# =======================================================================================================

optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters, optimization_model_part, CustomSU2Analyzer())
optimizer.Optimize()