# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.MappingApplication import *
from KratosMultiphysics.MeshMovingApplication import *

# Additional imports
from structural_response import StrainEnergyResponseFunction
from interface_su2 import InterfaceSU2
from analyzer_base import AnalyzerBaseClass
from custom_variable_utilities import ReadNodalVariableToDictionary
from mesh_moving_analysis import MeshMovingAnalysis
import os

# =======================================================================================================
# Preprocessign
# =======================================================================================================

with open("GENERAL_parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

interface_su2 = InterfaceSU2(parameters["su2_interface_settings"])
# interface_su2.WriteSU2MeshAsMDPA()

csm_model = Model()
csm_response = StrainEnergyResponseFunction("csm_response", parameters["csm_response_settings"], csm_model)

mesh_moving_analysis = MeshMovingAnalysis(csm_model, parameters["mesh_motion_settings"])
csm_model.GetModelPart("CSM_model_thin").AddNodalSolutionStepVariable(MESH_CHANGE)

# Change threading layer, otherwise the suprocess module used in the interface.py in SU2 will hang, (only)
# when Intel solvers are used as linear solvers for calculating the structure. This is a known issue of
# the Intel compiler 2018 and will be fixed in future releases. Compare:
# 1) https://github.com/numpy/numpy/issues/10060
# 2) https://software.intel.com/en-us/forums/intel-c-compiler/topic/758961
os.environ["MKL_THREADING_LAYER"] = "TBB"

# =======================================================================================================
# Define external analyzer
# =======================================================================================================
class CustomSU2Analyzer(AnalyzerBaseClass):

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        # interface_su2.InitializeNewSU2Project()
        csm_response.Initialize()
        mesh_moving_analysis.Initialize()

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        # # Compute CFD
        # if communicator.isRequestingValueOf("pressure_loss") or \
        #    communicator.isRequestingGradientOf("pressure_loss") or \
        #    communicator.isRequestingValueOf("strain_energy"):

        #     # Mesh update
        #     if optimization_iteration == 1:
        #         interface_su2.WriteNodesAsSU2MeshMotionFile(current_design.GetNodes())
        #     else:
        #         previos_iteration = int(optimization_iteration-1)
        #         interface_su2.WriteNodesAsSU2MeshMotionFile(current_design.GetNodes(),"DESIGNS/DSN_"+str(previos_iteration).zfill(3))

        #     # Caluclate value (primal field)
        #     update_mesh = True
        #     [value] = interface_su2.ComputeValues(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)

        #     if communicator.isRequestingValueOf("pressure_loss"):
        #         communicator.reportValue("pressure_loss", value)

        #     # Calculate gradient
        #     if communicator.isRequestingGradientOf("pressure_loss"):
        #         update_mesh = False
        #         [pressure_gradient] = interface_su2.ComputeGradient(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)
        #         communicator.reportGradient("pressure_loss", pressure_gradient)

        # Compute CSM
        if communicator.isRequestingValueOf("strain_energy") or \
           communicator.isRequestingGradientOf("strain_energy"):

            # Some variables
            csm_mdpa = csm_model.GetModelPart("CSM_model_thin")
            csm_support_mdpa = csm_mdpa.GetSubModelPart("supports")
            csm_outer_surface_mdpa = csm_mdpa.GetSubModelPart("outer_surface")
            csm_interface_mdpa = csm_mdpa.GetSubModelPart("wet_surface")
            cfd_initerface_mdpa = current_design

            # Mesh motion / initialize new shape
            VariableUtils().SetToZero_VectorVar(MESH_DISPLACEMENT,csm_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_X, True, csm_support_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Y, True, csm_support_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Z, True, csm_support_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_X, True, csm_outer_surface_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Y, True, csm_outer_surface_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Z, True, csm_outer_surface_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_X, True, csm_interface_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Y, True, csm_interface_mdpa.Nodes)
            VariableUtils().ApplyFixity(MESH_DISPLACEMENT_Z, True, csm_interface_mdpa.Nodes)

            shape_update_mapper = MapperFactory.CreateMapper(cfd_initerface_mdpa, csm_interface_mdpa, parameters["mapper_settings"].Clone())
            shape_update_mapper.Map(SHAPE_UPDATE, MESH_DISPLACEMENT)

            time_before_mesh_update = csm_mdpa.ProcessInfo.GetValue(TIME)
            step_before_mesh_update = csm_mdpa.ProcessInfo.GetValue(STEP)

            if not mesh_moving_analysis.time < mesh_moving_analysis.end_time:
                mesh_moving_analysis.end_time += 1
            mesh_moving_analysis.RunSolutionLoop()

            csm_mdpa.ProcessInfo.SetValue(TIME, time_before_mesh_update)
            csm_mdpa.ProcessInfo.SetValue(STEP, step_before_mesh_update)

            MeshControllerUtilities(csm_mdpa).SetReferenceMeshToMesh()
            MeshControllerUtilities(csm_mdpa).LogMeshChangeAccordingInputVariable(MESH_DISPLACEMENT)

            # Initialize new structural solution
            csm_response.InitializeSolutionStep()

            # Read / calculate fluid forces
            file_with_pressures = "DESIGNS/DSN_"+str(1).zfill(3)+"/DIRECT/surface_flow.csv"
            pressure_values = interface_su2.ReadNodalValueFromCSVFile(file_with_pressures,0,4,True)

            GeometryUtilities(cfd_initerface_mdpa).ComputeUnitSurfaceNormals()

            for node in cfd_initerface_mdpa.Nodes:
                area_normal = node.GetSolutionStepValue(NORMAL)
                force = -1.0 * pressure_values[node.Id] * area_normal
                node.SetSolutionStepValue(COUPLING_VARIABLE_1,force)

            # Mapping of forces
            force_mapper = MapperFactory.CreateMapper(cfd_initerface_mdpa, csm_interface_mdpa, parameters["mapper_settings"].Clone())
            force_mapper.Map(COUPLING_VARIABLE_1, POINT_LOAD)

            # Calculate value
            csm_response.CalculateValue()
            communicator.reportValue("strain_energy", csm_response.GetValue())

            # Calculate gradient
            if communicator.isRequestingGradientOf("strain_energy"):
                csm_response.CalculateGradient()

                sensitivity_mapper = MapperFactory.CreateMapper(csm_interface_mdpa, cfd_initerface_mdpa, parameters["mapper_settings"].Clone())
                sensitivity_mapper.Map(SHAPE_SENSITIVITY, DF1DX)

                gradient = ReadNodalVariableToDictionary(cfd_initerface_mdpa, DF1DX)
                communicator.reportGradient("strain_energy", gradient)

            csm_response.FinalizeSolutionStep()

            # Clean results on model parts
            MeshControllerUtilities(csm_mdpa).SetMeshToReferenceMesh()
            MeshControllerUtilities(csm_mdpa).SetDeformationVariablesToZero()

# =======================================================================================================
# Perform optimization
# =======================================================================================================

optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model_part, CustomSU2Analyzer())
optimizer.Optimize()