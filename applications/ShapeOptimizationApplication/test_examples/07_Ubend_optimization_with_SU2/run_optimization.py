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
from gid_output_process import GiDOutputProcess
from custom_variable_utilities import WriteDictionaryDataOnNodalVariable, ReadNodalVariableToDictionary
import optimizer_factory
import os

# =======================================================================================================
# Preprocessign
# =======================================================================================================

with open("GENERAL_parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

csm_model = Model()
csm_response = StrainEnergyResponseFunction("csm_response", parameters["csm_response_settings"], csm_model)

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
        # Read interface mdpa from fluid side
        self.cfd_mdpa = ModelPart("CFD_ubend")
        self.cfd_mdpa.ProcessInfo.SetValue(DOMAIN_SIZE, 3)
        self.cfd_mdpa.AddNodalSolutionStepVariable(DF1DX)
        self.cfd_mdpa.AddNodalSolutionStepVariable(POINT_LOAD)
        self.cfd_mdpa.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)
        self.cfd_mdpa.AddNodalSolutionStepVariable(MESH_CHANGE)
        self.cfd_mdpa.AddNodalSolutionStepVariable(NORMAL)
        self.cfd_mdpa.AddNodalSolutionStepVariable(NORMALIZED_SURFACE_NORMAL)
        model_part_io = ModelPartIO("CFD_ubend")
        model_part_io.ReadModelPart(self.cfd_mdpa)

        # Define some working variables
        self.cfd_interface_mdpa = self.cfd_mdpa.GetSubModelPart("BOGEN")
        self.csm_mdpa = csm_model.GetModelPart("CSM_model_thin")
        self.csm_support_mdpa = self.csm_mdpa.GetSubModelPart("supports")
        self.csm_outer_surface_mdpa = self.csm_mdpa.GetSubModelPart("outer_surface_with_opt_conditions")
        self.csm_interface_mdpa = self.csm_mdpa.GetSubModelPart("wet_surface_with_opt_conditions")

        # Initialize SU2 interface
        self.interface_su2 = InterfaceSU2(parameters["su2_interface_settings"])
        # self.interface_su2.WriteSU2MeshAsMDPA()
        self.interface_su2.InitializeNewSU2Project()

        # Initialize output of CFD interface
        self.cfd_interface_gid_output = GiDOutputProcess( self.cfd_interface_mdpa, "CFD_interface", parameters["cfd_interface_output_settings"] )
        self.cfd_interface_gid_output.ExecuteInitialize()
        self.cfd_interface_gid_output.ExecuteBeforeSolutionLoop()

        # Initialize damping of mesh displacement on CFD interface
        cfd_interface_damping_region_name = "fix_nodes"
        damping_regions = {}
        damping_regions[cfd_interface_damping_region_name] = self.cfd_mdpa.GetSubModelPart(cfd_interface_damping_region_name)
        modified_settings_for_damping = parameters["optimization_settings"].Clone()
        modified_settings_for_damping["design_variables"]["damping"]["damping_regions"][0]["sub_model_part_name"].SetString(cfd_interface_damping_region_name)
        self.cfd_interface_damping_utils = DampingUtilities(self.cfd_interface_mdpa, damping_regions, modified_settings_for_damping)

        # Initialize CSM solver
        csm_response.Initialize()

    # --------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        # Update fluid mesh (structure is controlled by the optimization algorithm)
        self.cfd_interface_mdpa.ProcessInfo[TIME] = optimization_iteration

        generalized_vm_mapper = MapperGeneralizedVertexMorphing(self.csm_interface_mdpa, self.cfd_interface_mdpa, parameters["optimization_settings"]["design_variables"]["filter"])
        generalized_vm_mapper.Map(CONTROL_POINT_UPDATE, MESH_DISPLACEMENT)

        if parameters["optimization_settings"]["design_variables"]["damping"]["perform_damping"].GetBool():
            self.cfd_interface_damping_utils.DampNodalVariable(MESH_DISPLACEMENT)

        MeshControllerUtilities(self.cfd_interface_mdpa).UpdateMeshAccordingInputVariable(MESH_DISPLACEMENT)
        MeshControllerUtilities(self.cfd_interface_mdpa).SetReferenceMeshToMesh()
        MeshControllerUtilities(self.cfd_interface_mdpa).LogMeshChangeAccordingInputVariable(MESH_DISPLACEMENT)

        # Compute CFD
        if communicator.isRequestingValueOf("pressure_loss") or \
           communicator.isRequestingGradientOf("pressure_loss") or \
           communicator.isRequestingValueOf("strain_energy"):

            # Write mesh update to SU2 case folder
            if optimization_iteration == 1:
                self.interface_su2.WriteNodesAsSU2MeshMotionFile(self.cfd_interface_mdpa.GetNodes())
            else:
                previos_iteration = int(optimization_iteration-1)
                self.interface_su2.WriteNodesAsSU2MeshMotionFile(self.cfd_interface_mdpa.GetNodes(),"DESIGNS/DSN_"+str(previos_iteration).zfill(3))

            # Caluclate value (primal field)
            update_mesh = True
            [value] = self.interface_su2.ComputeValues(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)

            if communicator.isRequestingValueOf("pressure_loss"):
                communicator.reportValue("pressure_loss", value)

            # Calculate gradient
            if communicator.isRequestingGradientOf("pressure_loss"):
                update_mesh = False
                [pressure_gradient] = self.interface_su2.ComputeGradient(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)

                WriteDictionaryDataOnNodalVariable(pressure_gradient, self.cfd_mdpa, DF1DX)
                if parameters["optimization_settings"]["design_variables"]["damping"]["perform_damping"].GetBool():
                    self.cfd_interface_damping_utils.DampNodalVariable(DF1DX)

                generalized_vm_mapper.InverseMap(DF1DX,DF1DX_MAPPED)

                mapped_pressure_gradient = ReadNodalVariableToDictionary(self.csm_interface_mdpa, DF1DX)
                communicator.reportGradient("pressure_loss", mapped_pressure_gradient)

        # Output results on CFD interface
        self.cfd_interface_gid_output.ExecuteInitializeSolutionStep()
        self.cfd_interface_gid_output.PrintOutput()
        self.cfd_interface_gid_output.ExecuteFinalizeSolutionStep()

        # Compute CSM
        if communicator.isRequestingValueOf("strain_energy") or \
           communicator.isRequestingGradientOf("strain_energy"):

            # Initialize new structural solution
            csm_response.InitializeSolutionStep()

            # Read / calculate fluid forces
            file_with_pressures = "DESIGNS/DSN_"+str(1).zfill(3)+"/DIRECT/surface_flow.csv"
            pressure_values = self.interface_su2.ReadNodalValueFromCSVFile(file_with_pressures,0,4,True)

            GeometryUtilities(self.cfd_interface_mdpa).ComputeUnitSurfaceNormals()

            for node in self.cfd_interface_mdpa.Nodes:
                area_normal = node.GetSolutionStepValue(NORMAL)
                force = -1.0 * pressure_values[node.Id] * area_normal
                node.SetSolutionStepValue(POINT_LOAD,force)

            # Mapping of forces
            force_mapper = MapperFactory.CreateMapper(self.cfd_interface_mdpa, self.csm_interface_mdpa, parameters["mapper_settings"].Clone())
            force_mapper.Map(POINT_LOAD, POINT_LOAD)

            # Calculate value
            csm_response.CalculateValue()
            communicator.reportValue("strain_energy", csm_response.GetValue())

            # Calculate gradient
            if communicator.isRequestingGradientOf("strain_energy"):
                csm_response.CalculateGradient()
                communicator.reportGradient("strain_energy", csm_response.GetShapeGradient())

            csm_response.FinalizeSolutionStep()

            # Clean results on model parts
            MeshControllerUtilities(self.csm_mdpa).SetMeshToReferenceMesh()
            MeshControllerUtilities(self.csm_mdpa).SetDeformationVariablesToZero()

 # =======================================================================================================
# Perform optimization
# =======================================================================================================
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], csm_model.GetModelPart("CSM_model_thin"), CustomSU2Analyzer())
optimizer.Optimize()