# Import Kratos core and apps
import KratosMultiphysics as km
import KratosMultiphysics.ShapeOptimizationApplication as kso
import KratosMultiphysics.StructuralMechanicsApplication as kcsm
import KratosMultiphysics.MappingApplication as kma

# Additional imports
from KratosMultiphysics.StructuralMechanicsApplication import structural_response_function_factory
from KratosMultiphysics.ShapeOptimizationApplication import mapper_factory
from KratosMultiphysics.ShapeOptimizationApplication import optimizer_factory
from KratosMultiphysics.ShapeOptimizationApplication.analyzer_base import AnalyzerBaseClass
from KratosMultiphysics.ShapeOptimizationApplication import custom_variable_utilities as cvu
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import time, os, shutil, csv, math
from decimal import Decimal

# =================================================================================================================================================================
# Preprocessing
# =================================================================================================================================================================

# Read parameters
with open("parameters_optimization.json",'r') as parameter_file:
    parameters = km.Parameters(parameter_file.read())

# Create Model
optimization_model = km.Model()

# Change threading layer, otherwise the suprocess module used in the interface.py in SU2 will hang, (only)
# when Intel solvers are used as linear solvers for calculating the structure. This is a known issue of
# the Intel compiler 2018 and will be fixed in future releases. Compare:
# 1) https://github.com/numpy/numpy/issues/10060
# 2) https://software.intel.com/en-us/forums/intel-c-compiler/topic/758961
os.environ["MKL_THREADING_LAYER"] = "TBB"

# =================================================================================================================================================================
# Define external analyzer
# =================================================================================================================================================================

class CustomAnalyzer(AnalyzerBaseClass):
    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __init__(self):
        self.cfd_part_name = "ubend_2d"
        self.cfd_part = None
        self.cfd_interface_part_name = "design_surface"
        self.cfd_interface_part = None

        self.cfd_wall_part_name = "WALL"
        self.cfd_wall_part = None

        self.csm_interface_part_name = "wet_interface"

        self.results_folder = "Optimization_Results"
        self.reference_pressure = 1#74571.25

        self.cfd_interface_gid_output = None
        self.cfd_interface_output_parameters = km.Parameters("""
        {
            "result_file_configuration" : {
                "gidpost_flags"       : {
                    "GiDPostMode"           : "GiD_PostBinary",
                    "WriteDeformedMeshFlag" : "WriteDeformed",
                    "WriteConditionsFlag"   : "WriteConditions",
                    "MultiFileFlag"         : "SingleFile"
                },
                "file_label"          : "step",
                "output_control_type" : "step",
                "output_frequency"    : 1,
                "body_output"         : true,
                "node_output"         : false,
                "skin_output"         : false,
                "nodal_results"       : ["TRACTION", "POINT_LOAD", "NORMAL", "DF1DX", "ADJOINT_DISPLACEMENT", "MESH_CHANGE", "MESH_DISPLACEMENT", "CFD_GRADIENT"],
                "gauss_point_results" : []
            },
            "point_data_configuration"  : []
        }""")
        self.adjoint_output_parameters = km.Parameters("""
        {
            "result_file_configuration" : {
                "gidpost_flags"       : {
                    "GiDPostMode"           : "GiD_PostBinary",
                    "WriteDeformedMeshFlag" : "WriteDeformed",
                    "WriteConditionsFlag"   : "WriteConditions",
                    "MultiFileFlag"         : "SingleFile"
                },
                "file_label"          : "step",
                "output_control_type" : "step",
                "output_frequency"    : 1,
                "body_output"         : true,
                "node_output"         : false,
                "skin_output"         : false,
                "nodal_results"       : ["ADJOINT_DISPLACEMENT","SHAPE_SENSITIVITY"],
                "gauss_point_results" : []
            },
            "point_data_configuration"  : []
        }""")
        self.traction_mapper_settings = km.Parameters("""
        {
            "mapper_type" : "nearest_element",
            "echo_level"  : 0
        }""")
        self.adjoint_displacement_mapper_settings = km.Parameters("""
        {
            "mapper_type" : "nearest_element",
            "echo_level"  : 0
        }""")
        self.stress_response_settings = km.Parameters("""
        {
                "response_type"        : "adjoint_max_stress",
                "gradient_mode"        : "semi_analytic",
                "step_size"            : 1e-11,
                "critical_part_name"   : "stress_partition",
                "stress_type"          : "VON_MISES_STRESS",
                "stress_treatment"     : "mean",
                "echo_level"           : 1,
                "primal_settings"      : "parameters_analysis.json",
                "adjoint_settings"     : "auto",
                "sensitivity_settings" : {
                    "sensitivity_model_part_name"     : "Parts_structure",
                    "nodal_sensitivity_variables"     : ["SHAPE_SENSITIVITY"],
                    "element_sensitivity_variables"   : [],
                    "condition_sensitivity_variables" : [],
                    "build_mode": "static"
                }
        }""")

        # Initialize SU2 interface
        # For running simulation on cluster (it is important to define the environmental variable before the import of the interface because as soon as SU2 is intialized it reads this variable in the file interface.py!
        # os.environ['SU2_MPI_COMMAND'] = 'mpirun -machinefile /home/danielb/my_simulations/05_ubend_coupled_volume_morphing/hosts.txt -n %i %s'
        from KratosMultiphysics.ShapeOptimizationApplication.interface_su2 import InterfaceSU2
        self.interface_su2 = InterfaceSU2(parameters["su2_interface_settings"])
        # self.interface_su2.WriteSU2MeshAsMDPA()
        self.interface_su2.InitializeNewSU2Project()

        # Settings from SU2 config
        self.RESTART_FLOW_FILENAME = "restart_flow.dat"
        self.SURFACE_FLOW_FILENAME = "surface_flow"

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        # Initialize cfd data
        cfd_model = km.Model()
        self.cfd_part = cfd_model.CreateModelPart(self.cfd_part_name)
        self.cfd_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,self.interface_su2.su2_mesh_data["NDIME"])
        self.cfd_part.AddNodalSolutionStepVariable(kso.TRACTION)
        self.cfd_part.AddNodalSolutionStepVariable(kcsm.POINT_LOAD)
        self.cfd_part.AddNodalSolutionStepVariable(kcsm.ADJOINT_DISPLACEMENT)
        self.cfd_part.AddNodalSolutionStepVariable(kso.DF1DX)
        self.cfd_part.AddNodalSolutionStepVariable(km.MESH_DISPLACEMENT)
        self.cfd_part.AddNodalSolutionStepVariable(kso.MESH_CHANGE)
        self.cfd_part.AddNodalSolutionStepVariable(km.NORMAL)
        self.cfd_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)
        self.cfd_part.AddNodalSolutionStepVariable(kso.CFD_GRADIENT)

        model_part_io = km.ModelPartIO(self.cfd_part_name)
        model_part_io.ReadModelPart(self.cfd_part)

        self.cfd_interface_part = self.cfd_part.GetSubModelPart(self.cfd_interface_part_name)
        self.cfd_wall_part = self.cfd_part.GetSubModelPart(self.cfd_wall_part_name)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        # obtain new csm mesh
        new_csm_mesh = {node.Id: [node.X, node.Y, node.Z] for node in current_design.Nodes}

        # Update fluid mesh (structure is controlled by the optimization algorithm)
        # Note that the mapper should be based on the previos design to match the forward map in the structure
        # Afterwards it hwas to be updated to current design to correctly map gradients of THIS optimization iteration
        self.__RemoveVariableFromCoordinates(kso.SHAPE_UPDATE, current_design)

        vm_csm2cfd_mapper = kso.MapperVertexMorphingMatrixFree(current_design, self.cfd_interface_part, parameters["optimization_settings"]["design_variables"]["filter"].Clone())
        vm_csm2cfd_mapper.Map(kso.CONTROL_POINT_UPDATE, km.MESH_DISPLACEMENT)

        self.__AddVariableToCoordinates(kso.SHAPE_UPDATE, current_design)

        vm_csm2cfd_mapper.Update()

        kso.MeshControllerUtilities(self.cfd_interface_part).UpdateMeshAccordingInputVariable(km.MESH_DISPLACEMENT)
        kso.MeshControllerUtilities(self.cfd_interface_part).SetReferenceMeshToMesh()
        kso.MeshControllerUtilities(self.cfd_interface_part).LogMeshChangeAccordingInputVariable(km.MESH_DISPLACEMENT)
        kso.GeometryUtilities(self.cfd_interface_part).ComputeUnitSurfaceNormals()

        if optimization_iteration == 1:
            self.interface_su2.WriteNodesAsSU2MeshMotionFile(self.cfd_wall_part.GetNodes())
        else:
            previos_iteration = int(optimization_iteration-1)
            self.interface_su2.WriteNodesAsSU2MeshMotionFile(self.cfd_wall_part.GetNodes(),"DESIGNS/DSN_"+str(previos_iteration).zfill(3))

        # Evaluate fluid
        if communicator.isRequestingValueOf("pressure_loss"):
            update_mesh = True
            [value] = self.interface_su2.ComputeValues(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)
            communicator.reportValue("pressure_loss", value)

        if communicator.isRequestingGradientOf("pressure_loss"):
            update_mesh = False
            [gradient] = self.interface_su2.ComputeGradient(["SURFACE_TOTAL_PRESSURE"], update_mesh, optimization_iteration)
            cvu.WriteDictionaryDataOnNodalVariable(gradient, self.cfd_part, kso.DF1DX)

            vm_csm2cfd_mapper.InverseMap(kso.DF1DX,kso.DF1DX_MAPPED)

            dummy_gradient = {node.Id: [0,0,0] for node in current_design.Nodes}
            communicator.reportGradient("pressure_loss", dummy_gradient)

        # Evaluate structure
        if communicator.isRequestingValueOf("mean_stress") or communicator.isRequestingGradientOf("mean_stress"):
            # Process fluid results
            self.__ReadAndDimensionalizeFluidResults(opt_itr=optimization_iteration)

            # Run CSM
            displacements, adjoint_displacements, value, csm_gradients = self.__RunCSM(opt_itr=optimization_iteration, new_mesh=new_csm_mesh, displacements={})
            communicator.reportValue("mean_stress", value)

            if communicator.isRequestingGradientOf("mean_stress"):
                # Map csm gradients
                cvu.WriteDictionaryDataOnNodalVariable(csm_gradients, current_design, kso.CSM_GRADIENT)
                vm_csm2csm_mapper = kso.MapperVertexMorphingMatrixFree(current_design, current_design, parameters["optimization_settings"]["design_variables"]["filter"].Clone())
                vm_csm2csm_mapper.InverseMap(kso.CSM_GRADIENT,kso.CSM_GRADIENT_MAPPED)

                # Run another adjoint fluid force analysis using the adjoint displacement as input
                update_mesh = False
                self.__ReplaceAdjointDisplacementInSU2RestartFile(adjoint_displacements, optimization_iteration)
                [drag_gradients] = self.interface_su2.ComputeGradient(["DRAG"], update_mesh, optimization_iteration)

                # Also dimensionalize cfd gradient as CFD is non-dimensional)
                drag_gradients = {key: [self.reference_pressure*value[0],self.reference_pressure*value[1],self.reference_pressure*value[2]] for key, value in drag_gradients.items()}
                cvu.WriteDictionaryDataOnNodalVariable(drag_gradients, self.cfd_part, kso.CFD_GRADIENT)

                vm_csm2cfd_mapper.InverseMap(kso.CFD_GRADIENT,kso.CFD_GRADIENT_MAPPED)

                self.__OutputCFDInterface(optimization_iteration)

                # Combine stress sensitivity parts
                for node in current_design.Nodes:
                    cfd_sens = node.GetSolutionStepValue(kso.CFD_GRADIENT_MAPPED)
                    csm_sens = node.GetSolutionStepValue(kso.CSM_GRADIENT_MAPPED)
                    combined = [cfd_sens[0]+csm_sens[0], cfd_sens[1]+csm_sens[1], cfd_sens[2]+csm_sens[2]]
                    node.SetSolutionStepValue(kso.DC1DX_MAPPED, combined)

                # Report dummy value (actual gradient is computed in analyzer here)
                dummy_gradient = {node.Id: [0,0,0] for node in current_design.Nodes}
                communicator.reportGradient("mean_stress", dummy_gradient)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def __RemoveVariableFromCoordinates(variable, model_part):
        for node in model_part.Nodes:
            shape_update = node.GetSolutionStepValue(variable)
            node.X -= shape_update[0]
            node.Y -= shape_update[1]
            node.Z -= shape_update[2]
            node.X0 -= shape_update[0]
            node.Y0 -= shape_update[1]
            node.Z0 -= shape_update[2]

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def __AddVariableToCoordinates(variable, model_part):
        for node in model_part.Nodes:
            shape_update = node.GetSolutionStepValue(variable)
            node.X += shape_update[0]
            node.Y += shape_update[1]
            node.Z += shape_update[2]
            node.X0 += shape_update[0]
            node.Y0 += shape_update[1]
            node.Z0 += shape_update[2]

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __ReadAndDimensionalizeFluidResults(self, opt_itr):
        file_with_pressures = "DESIGNS/DSN_"+str(opt_itr).zfill(3)+"/DIRECT/"+self.SURFACE_FLOW_FILENAME+".csv"
        if self.interface_su2.su2_mesh_data["NDIME"] == 2:
            nodal_forces = self.interface_su2.ReadNodalValuesFromCSVFile(file_with_pressures,0,[5,6],1)
            nodal_forces = {key: [value[0],value[1],0.0] for key, value in nodal_forces.items()}
            nodal_force_densities = self.interface_su2.ReadNodalValuesFromCSVFile(file_with_pressures,0,[7,8],1)
            nodal_force_densities = {key: [value[0],value[1],0.0] for key, value in nodal_force_densities.items()}
        else:
            nodal_forces = self.interface_su2.ReadNodalValuesFromCSVFile(file_with_pressures,0,[6,7,8],1)
            nodal_force_densities = self.interface_su2.ReadNodalValuesFromCSVFile(file_with_pressures,0,[9,10,11],1)

        for node in self.cfd_interface_part.Nodes:
            force = nodal_forces[node.Id]
            dimensional_force = [value*self.reference_pressure for value in force]
            node.SetSolutionStepValue(kcsm.POINT_LOAD, dimensional_force)

            traction = nodal_force_densities[node.Id]
            dimensional_traction = [value*self.reference_pressure for value in traction]
            node.SetSolutionStepValue(kso.TRACTION, dimensional_traction)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __RunCSM(self, opt_itr, new_mesh, displacements):
        print("\n> Starting __RunCSM in analyzer...")
        start_time = time.time()

        analysis_model = km.Model()

        csm = structural_response_function_factory.CreateResponseFunction(self.stress_response_settings["response_type"].GetString(), self.stress_response_settings, analysis_model)

        csm_primal_part = csm.primal_model_part
        csm_primal_part.AddNodalSolutionStepVariable(kso.MESH_CHANGE)
        csm_primal_part.AddNodalSolutionStepVariable(kso.NORMALIZED_SURFACE_NORMAL)
        csm_primal_part.AddNodalSolutionStepVariable(km.NORMAL)
        csm_primal_part.AddNodalSolutionStepVariable(kso.TRACTION)
        csm_adjoint_part = csm.adjoint_model_part

        # Initialize and set time to current optimization iteration
        csm.primal_analysis.project_parameters["problem_data"]["start_time"].SetDouble(opt_itr-1)
        csm.primal_analysis.project_parameters["problem_data"]["end_time"].SetDouble(opt_itr)

        csm.adjoint_analysis.project_parameters["problem_data"]["start_time"].SetDouble(opt_itr-1)
        csm.adjoint_analysis.project_parameters["problem_data"]["end_time"].SetDouble(opt_itr)

        csm.Initialize()

        csm_primal_part.ProcessInfo.SetValue(km.STEP, opt_itr-1)
        csm_adjoint_part.ProcessInfo.SetValue(km.STEP, opt_itr-1)

        # apply mesh motion
        if new_mesh != {}:
            for node in csm_primal_part.Nodes:
                X_new = new_mesh[node.Id]
                mesh_change = [X_new[0]-node.X, X_new[1]-node.Y, X_new[2]-node.Z]
                node.SetSolutionStepValue(kso.MESH_CHANGE, mesh_change)
                node.X = X_new[0]
                node.Y = X_new[1]
                node.Z = X_new[2]
                node.X0 = X_new[0]
                node.Y0 = X_new[1]
                node.Z0 = X_new[2]

        # compute primal field if necessary
        if displacements == {}:

            # Map forces
            csm_interface_part = csm_primal_part.GetSubModelPart(self.csm_interface_part_name)
            csm_interface_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,2)
            kso.GeometryUtilities(csm_interface_part).ComputeUnitSurfaceNormals()
            csm_interface_part.ProcessInfo.SetValue(km.DOMAIN_SIZE,3)

            cfd_to_csm_mapper = kma.MapperFactory.CreateMapper(csm_interface_part, self.cfd_interface_part, self.traction_mapper_settings.Clone())
            cfd_to_csm_mapper.InverseMap(kcsm.POINT_LOAD, kcsm.POINT_LOAD, kma.Mapper.USE_TRANSPOSE)

            # Apply forces
            # for node_i in csm_interface_part.Nodes:
            #     normal = node_i.GetSolutionStepValue(km.NORMAL)
            #     area = math.sqrt( normal[0]**2+normal[1]**2+normal[2]**2 )
            #     traction = node_i.GetSolutionStepValue(kso.TRACTION)
            #     force = [value*area for value in traction]
            #     node_i.SetSolutionStepValue(kso.TRACTION, traction)
            #     node_i.SetSolutionStepValue(kcsm.POINT_LOAD, force)

            # Compute primals
            csm.InitializeSolutionStep()

            # Store primal result
            displacements = {node.Id: node.GetSolutionStepValue(km.DISPLACEMENT) for node in csm_primal_part.Nodes}
        else:
            csm_analys = csm.primal_analysis
            csm_analys._GetSolver().AdvanceInTime(current_time=opt_itr-1)
            csm_analys.InitializeSolutionStep()

            for node in csm_primal_part.Nodes:
                node.SetSolutionStepValue(km.DISPLACEMENT, displacements[node.Id])

            csm_analys.FinalizeSolutionStep()
            csm_analys.OutputSolutionStep()

        csm.CalculateValue()
        csm.CalculateGradient()
        csm.FinalizeSolutionStep()
        csm.Finalize()

        # Map adjoint displacements to cfd side
        csm_adjoint_interface_part = csm_adjoint_part.GetSubModelPart(self.csm_interface_part_name)
        cfd_to_csm_mapper = kma.MapperFactory.CreateMapper(csm_adjoint_interface_part, self.cfd_interface_part, self.adjoint_displacement_mapper_settings.Clone())
        # cfd_to_csm_mapper.InverseMap(kcsm.ADJOINT_DISPLACEMENT, kcsm.ADJOINT_DISPLACEMENT, kma.Mapper.USE_TRANSPOSE) #Actually this is correct, but it leads to unexpected osciallations on CFD side --> approximation: mapping using consistent mapping matrix
        cfd_to_csm_mapper.Map(kcsm.ADJOINT_DISPLACEMENT, kcsm.ADJOINT_DISPLACEMENT)
        adjoint_displacements = {node.Id: node.GetSolutionStepValue(kcsm.ADJOINT_DISPLACEMENT) for node in self.cfd_interface_part.Nodes}

        # Some postprocessing
        primal_gid_results_filename = "primal_results.post.bin"
        adjoint_results_filename = "adjoint_results.post.bin"
        self.__OutputMdpaAsGid(csm.adjoint_model_part, adjoint_results_filename.replace(".post.bin",""), self.adjoint_output_parameters)
        self.__CopyBinFileToResultsFolder(primal_gid_results_filename, opt_itr)
        self.__CopyBinFileToResultsFolder(adjoint_results_filename, opt_itr)

        print("> Finished __RunCSM in" ,round( time.time()-start_time, 3 ), " s.")

        return displacements, adjoint_displacements, csm.GetValue(), csm.GetShapeGradient()

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __ReplaceAdjointDisplacementInSU2RestartFile(self, adjoint_displacements, optimization_iteration):
        su2_primal_restart_file = "DESIGNS/DSN_"+str(optimization_iteration).zfill(3)+"/"+self.RESTART_FLOW_FILENAME

        # Read restart file
        with open(su2_primal_restart_file, 'r') as infile:
            csv_reader = csv.reader(infile, delimiter='\t')
            lines = list(csv_reader)

        # Overwrite restart file
        with open(su2_primal_restart_file, 'w') as outfile:
            node_id_su2_to_kratos = self.interface_su2.node_id_su2_to_kratos
            num_nodes_in_su2_mesh = len(node_id_su2_to_kratos.keys())

            # Write header (already includes adjoint structure)
            for entry in lines[0][:-1]:
                outfile.write(entry+"\t")
            outfile.write(lines[0][-1]+"\n")

            # Write nodes
            for line in lines[1:num_nodes_in_su2_mesh+1]:
                su_node_id = int(line[0])
                kratos_node_id = node_id_su2_to_kratos[su_node_id]

                if kratos_node_id in adjoint_displacements.keys():
                    adjoint_disp = adjoint_displacements[kratos_node_id]
                else:
                    adjoint_disp = [0,0,0]

                line_to_write = line[:-1]

                if self.interface_su2.su2_mesh_data["NDIME"] == 2:
                    line_to_write[-2] = '%.15E' % Decimal(str(adjoint_disp[0]))
                    line_to_write[-1] = '%.15E' % Decimal(str(adjoint_disp[1]))
                else:
                    line_to_write[-3] = '%.15E' % Decimal(str(adjoint_disp[0]))
                    line_to_write[-2] = '%.15E' % Decimal(str(adjoint_disp[1]))
                    line_to_write[-1] = '%.15E' % Decimal(str(adjoint_disp[2]))

                # Write given quanitites
                for entry in line_to_write[:-1]:
                    outfile.write(entry+"\t")
                outfile.write(line_to_write[-1]+"\n")

            # Write appendix
            for line in lines[num_nodes_in_su2_mesh+1:]:
                outfile.write(line[0]+"\n")

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __OutputCFDInterface( self, optimization_iteration ):
        self.cfd_interface_part.ProcessInfo[km.TIME] = optimization_iteration
        if self.cfd_interface_gid_output is None:
            self.cfd_interface_gid_output = GiDOutputProcess( self.cfd_interface_part, os.path.join(self.results_folder,"CFD_interface"), self.cfd_interface_output_parameters)
            self.cfd_interface_gid_output.ExecuteInitialize()
            self.cfd_interface_gid_output.ExecuteBeforeSolutionLoop()
        self.cfd_interface_gid_output.ExecuteInitializeSolutionStep()
        self.cfd_interface_gid_output.PrintOutput()
        self.cfd_interface_gid_output.ExecuteFinalizeSolutionStep()

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __CopyBinFileToResultsFolder(self, filename, optimization_iteration):
        old_name = filename
        new_name = old_name.replace(".post.bin","_DESIGN_"+str(optimization_iteration)+".post.bin")
        os.rename(old_name,new_name)
        shutil.move(new_name, self.results_folder)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    def __CopyFileUnderNewNameToResultsFolder(self, filename, new_name):
        os.rename(filename,new_name)
        shutil.move(new_name, self.results_folder)

    # -------------------------------------------------------------------------------------------------------------------------------------------------------------
    @staticmethod
    def __OutputMdpaAsGid( output_mdpa, output_filename, parameters ):
        gid_output_original = GiDOutputProcess( output_mdpa, output_filename, parameters )
        gid_output_original.ExecuteInitialize()
        gid_output_original.ExecuteBeforeSolutionLoop()
        gid_output_original.ExecuteInitializeSolutionStep()
        gid_output_original.PrintOutput()
        gid_output_original.ExecuteFinalizeSolutionStep()
        gid_output_original.ExecuteFinalize()

# =================================================================================================================================================================
# Perform optimization
# =================================================================================================================================================================

optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model, CustomAnalyzer())
optimizer.Optimize()