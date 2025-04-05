import time as timer
import os
import sys
import pathlib
import math
import numpy as np
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.DEMApplication.DEM_restart_utility import DEMRestartUtility
import KratosMultiphysics.DEMApplication.dem_default_input_parameters
from KratosMultiphysics.DEMApplication.analytic_tools import analytic_data_procedures
from KratosMultiphysics.DEMApplication.materials_assignation_utility import MaterialsAssignationUtility
from KratosMultiphysics.DEMApplication.DEM_properties_measure_utility import DEMPropertiesMeasureUtility

from importlib import import_module

if IsDistributedRun():
    if "DO_NOT_PARTITION_DOMAIN" in os.environ:
        Logger.PrintInfo("DEM", "Running under MPI........")
        from KratosMultiphysics.mpi import *
        import KratosMultiphysics.DEMApplication.DEM_procedures_mpi_no_partitions as DEM_procedures
    else:
        Logger.PrintInfo("DEM", "Running under OpenMP........")
        from KratosMultiphysics.MetisApplication import *
        from KratosMultiphysics.MPISearchApplication import *
        from KratosMultiphysics.mpi import *
        import KratosMultiphysics.DEMApplication.DEM_procedures_mpi as DEM_procedures
else:
    Logger.PrintInfo("DEM", "Running under OpenMP........")
    import KratosMultiphysics.DEMApplication.DEM_procedures as DEM_procedures

class DEMAnalysisStage(AnalysisStage):

    def GetParametersFileName(self):
        return "ProjectParametersDEM.json"

    def GetInputParameters(self):
        self.KratosPrintWarning('Warning: Calls to this method (GetInputParameters) will become deprecated in the near future.')
        parameters_file_name = self.GetParametersFileName()
        parameters_file = open(parameters_file_name, 'r')
        return Parameters(parameters_file.read())

    def LoadParametersFile(self):
        self.KratosPrintWarning('Warning: Calls to this method (LoadParametersFile) will become deprecated in the near future.')
        self.DEM_parameters = self.GetInputParameters()
        self.project_parameters = self.DEM_parameters
        default_input_parameters = self.GetDefaultInputParameters()
        self.DEM_parameters.ValidateAndAssignDefaults(default_input_parameters)
        self.FixParametersInconsistencies()

    def FixParametersInconsistencies(self): # TODO: This is here to avoid inconsistencies until the jsons become standard
        final_time = self.DEM_parameters["FinalTime"].GetDouble()
        problem_name = self.DEM_parameters["problem_name"].GetString()
        self.project_parameters["problem_data"]["end_time"].SetDouble(final_time)
        self.project_parameters["problem_data"]["problem_name"].SetString(problem_name)

    def GetDefaultInputParameters(self):
        return KratosMultiphysics.DEMApplication.dem_default_input_parameters.GetDefaultInputParameters()

    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid, IO.SKIP_TIMER)

    def GetMainPath(self):
        return os.getcwd()

    def __init__(self, model, DEM_parameters):
        self.model = model
        self.main_path = self.GetMainPath()
        self.mdpas_folder_path = self.main_path

        self.DEM_parameters = DEM_parameters    # TODO, can be improved
        self.project_parameters = DEM_parameters
        default_input_parameters = self.GetDefaultInputParameters()
        self.DEM_parameters.ValidateAndAssignDefaults(default_input_parameters)
        self.FixParametersInconsistencies()

        self.do_print_results_option = self.DEM_parameters["do_print_results_option"].GetBool()
        if not "WriteMdpaFromResults" in self.DEM_parameters.keys():
            self.write_mdpa_from_results = False
        else:
            self.write_mdpa_from_results = self.DEM_parameters["WriteMdpaFromResults"].GetBool()
        self.creator_destructor = self.SetParticleCreatorDestructor(DEM_parameters["creator_destructor_settings"])
        self.dem_fem_search = self.SetDemFemSearch()
        self.procedures = self.SetProcedures()
        self.PreUtilities = PreUtilities()

        # Set the print function TO_DO: do this better...
        self.KratosPrintInfo = self.procedures.KratosPrintInfo

        # Creating necessary directories:
        self.problem_name = self.GetProblemTypeFileName()

        [self.post_path, self.graphs_path] = self.procedures.CreateDirectories(str(self.main_path), str(self.problem_name), do_print_results=self.do_print_results_option)

        # Prepare modelparts
        self.CreateModelParts()

        if self.do_print_results_option:
            self.SetGraphicalOutput()
        self.report = DEM_procedures.Report()
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.translational_scheme = self.SetTranslationalScheme()
        self.rotational_scheme = self.SetRotationalScheme()

        # Define control variables
        self.p_frequency = 100   # activate every 100 steps
        self.step_count = 0
        self.p_count = self.p_frequency

        #self._solver = self._GetSolver()
        self.SetFinalTime()
        self.AddVariables()

        super().__init__(model, self.DEM_parameters)

    def CreateModelParts(self):
        self.spheres_model_part = self.model.CreateModelPart("SpheresPart")
        self.rigid_face_model_part = self.model.CreateModelPart("RigidFacePart")
        self.cluster_model_part = self.model.CreateModelPart("ClusterPart")
        self.dem_inlet_model_part = self.model.CreateModelPart("DEMInletPart")
        self.mapping_model_part = self.model.CreateModelPart("MappingPart")
        self.contact_model_part = self.model.CreateModelPart("ContactPart")

        mp_list = []
        mp_list.append(self.spheres_model_part)
        mp_list.append(self.rigid_face_model_part)
        mp_list.append(self.cluster_model_part)
        mp_list.append(self.dem_inlet_model_part)
        mp_list.append(self.mapping_model_part)
        mp_list.append(self.contact_model_part)

        self.all_model_parts = DEM_procedures.SetOfModelParts(mp_list)

    def IsCountStep(self):
        self.step_count += 1
        if self.step_count == self.p_count:
            self.p_count += self.p_frequency
            return True

        return False

    def SetAnalyticWatchers(self):
        self.SurfacesAnalyzerClass = analytic_data_procedures.SurfacesAnalyzerClass(self.rigid_face_model_part.SubModelParts, self.main_path)

        if self.post_normal_impact_velocity_option:
            self.ParticlesAnalyzerClass = analytic_data_procedures.ParticlesAnalyzerClass(self.analytic_model_part)

    def MakeAnalyticsMeasurements(self):
        self.SurfacesAnalyzerClass.MakeAnalyticsMeasurements()

        if self.post_normal_impact_velocity_option:
            self.ParticlesAnalyzerClass.MakeAnalyticsMeasurements()

    def SetFinalTime(self):
        self.end_time = self.DEM_parameters["FinalTime"].GetDouble()

    def SetProcedures(self):
        return DEM_procedures.Procedures(self.DEM_parameters)

    def SetDemFemSearch(self):
        return DEM_FEM_Search()

    def GetParticleHistoryWatcher(self):
        return None

    def SetParticleCreatorDestructor(self, creator_destructor_settings):

        self.watcher = self.GetParticleHistoryWatcher()

        if self.watcher is None:
            return ParticleCreatorDestructor(creator_destructor_settings)
        return ParticleCreatorDestructor(self.watcher, creator_destructor_settings)

    def SelectTranslationalScheme(self):
        if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Forward_Euler':
            return ForwardEulerScheme()
        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Symplectic_Euler':
            return SymplecticEulerScheme()
        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Taylor_Scheme':
            return TaylorScheme()
        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet':
            return VelocityVerletScheme()

        return None

    def SelectRotationalScheme(self):
        if self.DEM_parameters["RotationalIntegrationScheme"].GetString() == 'Direct_Integration':
            return self.SelectTranslationalScheme()
        elif self.DEM_parameters["RotationalIntegrationScheme"].GetString() == 'Runge_Kutta':
            return RungeKuttaScheme()
        elif self.DEM_parameters["RotationalIntegrationScheme"].GetString() == 'Quaternion_Integration':
            return QuaternionIntegrationScheme()

        return None

    def SetTranslationalScheme(self):
        translational_scheme = self.SelectTranslationalScheme()

        if translational_scheme is None:
            raise Exception('Error: selected translational integration scheme not defined. Please select a different scheme')

        return translational_scheme

    def SetRotationalScheme(self):
        rotational_scheme = self.SelectRotationalScheme()

        if rotational_scheme is None:
            raise Exception('Error: selected rotational integration scheme not defined. Please select a different scheme')

        return rotational_scheme

    def SetSolver(self):        # TODO why is this still here. -> main_script calls retrocompatibility
        return self._CreateSolver()

    def _CreateSolver(self):
        def SetSolverStrategy():
            strategy_file_name = self.DEM_parameters["solver_settings"]["strategy"].GetString()
            imported_module = import_module("KratosMultiphysics.DEMApplication" + "." + strategy_file_name)
            return imported_module

        return SetSolverStrategy().ExplicitStrategy(self.all_model_parts,
                                                    self.creator_destructor,
                                                    self.dem_fem_search,
                                                    self.DEM_parameters,
                                                    self.procedures)

    def AddVariables(self):
        self.procedures.AddAllVariablesInAllModelParts(self._GetSolver(), self.translational_scheme, self.rotational_scheme, self.all_model_parts, self.DEM_parameters)

    def FillAnalyticSubModelParts(self):
        if not self.spheres_model_part.HasSubModelPart("AnalyticParticlesPart"):
            self.spheres_model_part.CreateSubModelPart('AnalyticParticlesPart')
        self.analytic_model_part = self.spheres_model_part.GetSubModelPart('AnalyticParticlesPart')
        analytic_particle_ids = [elem.Id for elem in self.spheres_model_part.Elements]
        self.analytic_model_part.AddElements(analytic_particle_ids)
        analytic_node_ids = [node.Id for node in self.spheres_model_part.Nodes]
        self.analytic_model_part.AddNodes(analytic_node_ids)

    def FillAnalyticSubModelPartsWithNewParticles(self):
        self.analytic_model_part = self.spheres_model_part.GetSubModelPart('AnalyticParticlesPart')
        self.PreUtilities.FillAnalyticSubModelPartUtility(self.spheres_model_part, self.analytic_model_part)
        #analytic_particle_ids = [elem.Id for elem in self.spheres_model_part.Elements]
        #self.analytic_model_part.AddElements(analytic_particle_ids)

    def Initialize(self):
        self.time = 0.0
        self.time_old_print = 0.0
        self.time_old_update_contact_element = 0.0

        self.ReadModelParts()

        self.SetMaterials()

        self.post_normal_impact_velocity_option = False
        if "PostNormalImpactVelocity" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostNormalImpactVelocity"].GetBool():
                self.post_normal_impact_velocity_option = True
                self.FillAnalyticSubModelParts()

        self.SetAnalyticWatchers()

        # Setting up the buffer size
        self.procedures.SetUpBufferSizeInAllModelParts(self.spheres_model_part, 1, self.cluster_model_part, 1, self.dem_inlet_model_part, 1, self.rigid_face_model_part, 1)

        self.KratosPrintInfo("Initializing Problem...")

        self.GraphicalOutputInitialize()

        # Perform a partition to balance the problem
        self.SetSearchStrategy()

        self.SolverBeforeInitialize()

        self.parallelutils.Repart(self.spheres_model_part)

        #Setting up the BoundingBox
        self.bounding_box_time_limits = self.procedures.SetBoundingBoxLimits(self.all_model_parts, self.creator_destructor)

        self.creator_destructor.SetMaxNodeId(self.all_model_parts.MaxNodeId)

        self.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(self.DEM_parameters, self.graphs_path, self.spheres_model_part, self.rigid_face_model_part)

        self.DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(self.DEM_parameters, self.spheres_model_part, self.cluster_model_part, self.graphs_path, "EnergyPlot.grf")

        self.DEMPropertiesMeasureUtility = DEMPropertiesMeasureUtility(self.DEM_parameters, self.spheres_model_part, self.contact_model_part, self.graphs_path)

        self.KratosPrintInfo("Initialization Complete")

        self.report.Prepare(timer, self.DEM_parameters["ControlTime"].GetDouble())

        self.post_utils = DEM_procedures.PostUtils(self.DEM_parameters, self.spheres_model_part)
        self.report.total_steps_expected = int(self.end_time / self._GetSolver().dt)

        super().Initialize()

        self.seed = self.DEM_parameters["seed"].GetInt()
        #Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation
        #Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
        self.SetInlet()

        self.SetInitialNodalValues()

        self.KratosPrintInfo(self.report.BeginReport(timer))

        if self.DEM_parameters["output_configuration"]["print_number_of_neighbours_histogram"].GetBool():
            self.PreUtilities.PrintNumberOfNeighboursHistogram(self.spheres_model_part, os.path.join(self.graphs_path, "number_of_neighbours_histogram.txt"))

        self.BoundingBoxMinX_update = self.DEM_parameters["BoundingBoxMinX"].GetDouble()
        self.BoundingBoxMinY_update = self.DEM_parameters["BoundingBoxMinY"].GetDouble()
        self.BoundingBoxMinZ_update = self.DEM_parameters["BoundingBoxMinZ"].GetDouble()
        self.BoundingBoxMaxX_update = self.DEM_parameters["BoundingBoxMaxX"].GetDouble()
        self.BoundingBoxMaxY_update = self.DEM_parameters["BoundingBoxMaxY"].GetDouble()
        self.BoundingBoxMaxZ_update = self.DEM_parameters["BoundingBoxMaxZ"].GetDouble()

        self.bounding_box_servo_loading_option = False
        if "BoundingBoxServoLoadingOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["BoundingBoxServoLoadingOption"].GetBool():
                self.bounding_box_servo_loading_option = True
        
        self.bounding_box_move_velocity = [0.0, 0.0, 0.0]

    def SetMaterials(self):

        self.ReadMaterialsFile()

        model_part_import_settings = self.DEM_parameters["solver_settings"]["model_import_settings"]
        input_type = model_part_import_settings["input_type"].GetString()
        if input_type == "rest":
            return

        materials_setter = MaterialsAssignationUtility(self.model, self.spheres_model_part, self.DEM_material_parameters)
        materials_setter.AssignMaterialParametersToProperties()
        materials_setter.AssignPropertiesToEntities()

    def ReadMaterialsFile(self):
        adapted_to_current_os_relative_path = pathlib.Path(self.DEM_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString())
        materials_file_abs_path = os.path.join(self.main_path, str(adapted_to_current_os_relative_path))
        with open(materials_file_abs_path, 'r') as materials_file:
            self.DEM_material_parameters = Parameters(materials_file.read())

    def SetSearchStrategy(self):
        self._GetSolver().search_strategy = self.parallelutils.GetSearchStrategy(self._GetSolver(), self.spheres_model_part)

    def SolverBeforeInitialize(self):
        self._GetSolver().BeforeInitialize()

    def SolverInitialize(self):
        self._GetSolver().Initialize() # Possible modifications of number of elements and number of nodes

    def GetProblemNameWithPath(self):
        return os.path.join(self.mdpas_folder_path, self.DEM_parameters["problem_name"].GetString())

    def GetDiscreteElementsInputFileTag(self):
        return 'DEM'

    def GetDEMInletInputFileTag(self):
        return 'DEM_Inlet'

    def GetDEMWallsInputFileTag(self):
        return 'DEM_FEM_boundary'

    def GetDEMClustersInputFileTag(self):
        return 'DEM_Clusters'

    def GetDiscreteElementsInputFilePath(self):
        return self.GetInputFilePath(self.GetDiscreteElementsInputFileTag())

    def GetDEMInletInputFilePath(self):
        return self.GetInputFilePath(self.GetDEMInletInputFileTag())

    def GetDEMWallsInputFilePath(self):
        return self.GetInputFilePath(self.GetDEMWallsInputFileTag())

    def GetDEMClustersInputFilePath(self):
        return self.GetInputFilePath(self.GetDEMClustersInputFileTag())

    def GetMpFilePath(self):
        return GetInputFilePath('DEM')

    def GetInletFilePath(self):
        return GetInputFilePath('DEM_Inlet')

    def GetFemFilePath(self):
        return GetInputFilePath('DEM_FEM_boundary')

    def GetClusterFilePath(self):
        return GetInputFilePath('DEM_Clusters')

    def GetInputFilePath(self, file_tag=''):
        return self.GetProblemNameWithPath() + file_tag

    def GetProblemTypeFileName(self):
        return self.DEM_parameters["problem_name"].GetString()

    def ReadModelPartsFromRestartFile(self, model_part_import_settings):
        Logger.PrintInfo('DEM', 'Loading model parts from restart file...')
        DEMRestartUtility(self.model, self._GetSolver()._GetRestartSettings(model_part_import_settings)).LoadRestart()
        Logger.PrintInfo('DEM', 'Finished loading model parts from restart file.')

    def ReadModelPartsFromMdpaFile(self, max_node_id, max_elem_id, max_cond_id):
        def UpdateMaxIds(max_node_id, max_elem_id, max_cond_id, model_part):
            max_node_id = max(max_node_id, self.creator_destructor.FindMaxNodeIdInModelPart(model_part))
            max_elem_id = max(max_elem_id, self.creator_destructor.FindMaxElementIdInModelPart(model_part))
            max_cond_id = max(max_cond_id, self.creator_destructor.FindMaxConditionIdInModelPart(model_part))
            return max_node_id, max_elem_id, max_cond_id

        def ReadModelPart(model_part, file_tag, max_node_id, max_elem_id, max_cond_id):
            file_path = self.GetInputFilePath(file_tag)

            if not os.path.isfile(file_path + '.mdpa'):
                self.KratosPrintInfo('Input file ' + file_tag + '.mdpa' + ' not found. Continuing.')
                return

            if model_part.Name == 'SpheresPart':
                model_part_io = self.model_part_reader(file_path, max_node_id, max_elem_id, max_cond_id)

                do_perform_initial_partition = True
                if self.DEM_parameters.Has("do_not_perform_initial_partition"):
                    if self.DEM_parameters["do_not_perform_initial_partition"].GetBool():
                        do_perform_initial_partition = False

                if do_perform_initial_partition:
                    self.parallelutils.PerformInitialPartition(model_part_io)

                model_part_io, model_part = self.parallelutils.SetCommunicator(
                                                model_part,
                                                model_part_io,
                                                file_path)
            else:
                model_part_io = self.model_part_reader(file_path, max_node_id + 1, max_elem_id + 1, max_cond_id + 1)

            model_part_io.ReadModelPart(model_part)

        ReadModelPart(self.spheres_model_part, self.GetDiscreteElementsInputFileTag(), max_node_id, max_elem_id, max_cond_id)
        max_node_id, max_elem_id, max_cond_id = UpdateMaxIds(max_node_id, max_elem_id, max_cond_id, self.spheres_model_part)
        old_max_elem_id_spheres = max_elem_id

        ReadModelPart(self.rigid_face_model_part, self.GetDEMWallsInputFileTag(), max_node_id, max_elem_id, max_cond_id)

        max_node_id, max_elem_id, max_cond_id = UpdateMaxIds(max_node_id, max_elem_id, max_cond_id, self.rigid_face_model_part)

        ReadModelPart(self.cluster_model_part, self.GetDEMClustersInputFileTag(), max_node_id, max_elem_id, max_cond_id)

        max_elem_id = self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)

        # Clusters generate extra spheres, so the following step is necessary
        if max_elem_id != old_max_elem_id_spheres:
            self.creator_destructor.RenumberElementIdsFromGivenValue(self.cluster_model_part, max_elem_id)

        max_node_id, max_elem_id, max_cond_id = UpdateMaxIds(max_node_id, max_elem_id, max_cond_id, self.cluster_model_part)

        ReadModelPart(self.dem_inlet_model_part, self.GetDEMInletInputFileTag(), max_node_id, max_elem_id, max_cond_id)

    def ReadModelParts(self, max_node_id=0, max_elem_id=0, max_cond_id=0):

        model_part_import_settings = self.DEM_parameters["solver_settings"]["model_import_settings"]
        input_type = model_part_import_settings["input_type"].GetString()

        if input_type == "rest":
            self.ReadModelPartsFromRestartFile(model_part_import_settings)
        elif input_type == "mdpa":
            self.ReadModelPartsFromMdpaFile(max_node_id, max_elem_id, max_cond_id)
        else:
            raise Exception('DEM', 'Model part input option \'' + input_type + '\' is not yet implemented.')

        self.model_parts_have_been_read = True
        self.all_model_parts.ComputeMaxIds()
        self.ConvertClusterFileNamesFromRelativePathToAbsolutePath()
        self.CheckConsistencyOfElementsAndNodesInEverySubModelPart()

    def CheckConsistencyOfElementsAndNodesInEverySubModelPart(self):
        def ErrorMessage(name):
            raise Exception(" ModelPart (or SubModelPart) "+ name + " has a different number of nodes and elements (particles)! \n")

        if self.spheres_model_part.NumberOfNodes(0) != self.spheres_model_part.NumberOfElements(0):
            ErrorMessage(self.spheres_model_part.Name)
        if self.cluster_model_part.NumberOfNodes(0) != self.cluster_model_part.NumberOfElements(0):
            ErrorMessage(self.cluster_model_part.Name)

        for submp in self.spheres_model_part.SubModelParts:
            if submp.NumberOfNodes(0) != submp.NumberOfElements(0):
                ErrorMessage(submp.Name)

        for submp in self.cluster_model_part.SubModelParts:
            if submp.NumberOfNodes(0) != submp.NumberOfElements(0):
                ErrorMessage(submp.Name)

    def ConvertClusterFileNamesFromRelativePathToAbsolutePath(self):
        for properties in self.cluster_model_part.Properties:
            if properties.Has(CLUSTER_FILE_NAME):
                cluster_file_name = properties[CLUSTER_FILE_NAME]
                properties[CLUSTER_FILE_NAME] = os.path.join(self.main_path, cluster_file_name)

        for submp in self.dem_inlet_model_part.SubModelParts:
            if submp.Has(CLUSTER_FILE_NAME):
                cluster_file_name = submp[CLUSTER_FILE_NAME]
                submp[CLUSTER_FILE_NAME] = os.path.join(self.main_path, cluster_file_name)

    def RunAnalytics(self, time):
        self.MakeAnalyticsMeasurements()
        if self.IsTimeToPrintPostProcess():
            self.SurfacesAnalyzerClass.MakeAnalyticsPipeLine(time)

        if self.post_normal_impact_velocity_option and self.IsTimeToPrintPostProcess():
            self.ParticlesAnalyzerClass.SetNodalMaxImpactVelocities()
            self.ParticlesAnalyzerClass.SetNodalMaxFaceImpactVelocities()

    def IsTimeToPrintPostProcess(self):
        return self.do_print_results_option and self.DEM_parameters["OutputTimeStep"].GetDouble() - (self.time - self.time_old_print) < 1e-2 * self._GetSolver().dt

    def IsTimeToUpdateContactElementForServo(self):
        return self.bounding_box_servo_loading_option and self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingFrequency"].GetInt() * \
                self._GetSolver().dt - (self.time - self.time_old_update_contact_element) < 1e-2 * self._GetSolver().dt

    def PrintResults(self):
        #### GiD IO ##########################################
        if self.IsTimeToPrintPostProcess():
            self.PrintResultsForGid(self.time)
            self.time_old_print = self.time

    def SolverSolve(self):
        self._GetSolver().SolveSolutionStep()

    def SetInlet(self):
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            #Constructing the inlet and initializing it (must be done AFTER the self.spheres_model_part Initialize)
            self.DEM_inlet = DEM_Inlet(self.dem_inlet_model_part, self.DEM_parameters["dem_inlets_settings"], self.seed)
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor, self._GetSolver().continuum_type)

    def SetInitialNodalValues(self):
        self.procedures.SetInitialNodalValues(self.spheres_model_part, self.cluster_model_part, self.dem_inlet_model_part, self.rigid_face_model_part)

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        if self.post_normal_impact_velocity_option:
            if self.IsCountStep():
                self.FillAnalyticSubModelPartsWithNewParticles()
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self.UpdateIsTimeToPrintInModelParts(self.IsTimeToPrintPostProcess())
        if self.bounding_box_servo_loading_option:
            self.UpdateIsTimeToUpdateContactElementForServo(self.IsTimeToUpdateContactElementForServo())

        if self.DEM_parameters["Dimension"].GetInt() == 2:
            self.spheres_model_part.ProcessInfo[IMPOSED_Z_STRAIN_OPTION] = self.DEM_parameters["ImposeZStrainIn2DOption"].GetBool()
            if not self.DEM_parameters["ImposeZStrainIn2DWithControlModule"].GetBool():
                if self.spheres_model_part.ProcessInfo[IMPOSED_Z_STRAIN_OPTION]:
                    self.spheres_model_part.ProcessInfo.SetValue(IMPOSED_Z_STRAIN_VALUE, eval(self.DEM_parameters["ZStrainValue"].GetString()))

        bounding_box_servo_loading_option = False
        if "BoundingBoxServoLoadingOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["BoundingBoxServoLoadingOption"].GetBool():
                bounding_box_servo_loading_option = True

        if "BoundingBoxMoveOption" in self.DEM_parameters.keys():
            if self.DEM_parameters["BoundingBoxMoveOption"].GetBool():
                time_step = self.spheres_model_part.ProcessInfo[TIME_STEPS]
                if bounding_box_servo_loading_option:
                    NStepSearch = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingFrequency"].GetInt()
                    if (time_step + 1) % NStepSearch == 0 and (time_step > 0):
                        measured_global_stress = self.MeasureSphereForGettingGlobalStressTensor()
                        self.CalculateBoundingBoxMoveVelocity(measured_global_stress)
                        self.UpdateSearchStartegyAndCPlusPlusStrategy(self.bounding_box_move_velocity)
                        self.procedures.UpdateBoundingBox(self.spheres_model_part, self.creator_destructor, self.bounding_box_move_velocity)
                else:
                    NStepSearch = self.DEM_parameters["NeighbourSearchFrequency"].GetInt()
                    if (time_step + 1) % NStepSearch == 0 and (time_step > 0):
                        bounding_box_move_velocity = self.DEM_parameters["BoundingBoxMoveVelocity"].GetVector()
                        self.UpdateSearchStartegyAndCPlusPlusStrategy(bounding_box_move_velocity)
                        self.procedures.UpdateBoundingBox(self.spheres_model_part, self.creator_destructor, bounding_box_move_velocity)

    def CalculateBoundingBoxMoveVelocity(self, measured_global_stress):

        # note: servo_loading_type = "isotropic" or "anisotropic"
        servo_loading_type = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingType"].GetString()

        if servo_loading_type == "isotropic":
            self.CalculateBoundingBoxMoveVelocityIsotropic(measured_global_stress)
        elif servo_loading_type == "anisotropic":
            self.CalculateBoundingBoxMoveVelocityAnisotropic(measured_global_stress)
        elif servo_loading_type == "triaxail_loading":
            self.CalculateBoundingBoxMoveVelocityTriaxailLoading(measured_global_stress)
    
    def CalculateBoundingBoxMoveVelocityIsotropic(self, measured_global_stress):

        servo_loading_stress = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingStress"].GetVector()
        servo_loading_factor = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingFactor"].GetDouble()
        d50 = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["MeanParticleDiameterD50"].GetDouble()
        Youngs_modulus = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["ParticleYoungsModulus"].GetDouble()

        servo_loading_stress_mean = (servo_loading_stress[0] + servo_loading_stress[1] + servo_loading_stress[2]) / 3.0
        measured_global_stress_mean = (measured_global_stress[0][0] + measured_global_stress[1][1] + measured_global_stress[2][2]) / 3.0

        delta_time = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        servo_loading_coefficient = servo_loading_factor * d50 / (delta_time * Youngs_modulus)
        self.bounding_box_move_velocity[0] = (servo_loading_stress_mean - measured_global_stress_mean) * servo_loading_coefficient

        bonuding_box_serva_loading_velocity_max = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingVelocityMax"].GetDouble()
        if abs(self.bounding_box_move_velocity[0]) > bonuding_box_serva_loading_velocity_max:
            self.bounding_box_move_velocity[0] = bonuding_box_serva_loading_velocity_max * np.sign(self.bounding_box_move_velocity[0])
        
        self.bounding_box_move_velocity[1] = self.bounding_box_move_velocity[0]
        self.bounding_box_move_velocity[2] = self.bounding_box_move_velocity[0]

    def CalculateBoundingBoxMoveVelocityAnisotropic(self, measured_global_stress):

        servo_loading_stress = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingStress"].GetVector()
        servo_loading_factor = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingFactor"].GetDouble()
        d50 = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["MeanParticleDiameterD50"].GetDouble()
        Youngs_modulus = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["ParticleYoungsModulus"].GetDouble()
        
        delta_time = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        servo_loading_coefficient = servo_loading_factor * d50 / (delta_time * Youngs_modulus)
        self.bounding_box_move_velocity[0] = (servo_loading_stress[0] - measured_global_stress[0][0]) * servo_loading_coefficient
        self.bounding_box_move_velocity[1] = (servo_loading_stress[1] - measured_global_stress[1][1]) * servo_loading_coefficient
        self.bounding_box_move_velocity[2] = (servo_loading_stress[2] - measured_global_stress[2][2]) * servo_loading_coefficient

        bonuding_box_move_velocity_max = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingVelocityMax"].GetDouble()
        if abs(self.bounding_box_move_velocity[0]) > bonuding_box_move_velocity_max:
            self.bounding_box_move_velocity[0] = bonuding_box_move_velocity_max * np.sign(self.bounding_box_move_velocity[0])
        if abs(self.bounding_box_move_velocity[1]) > bonuding_box_move_velocity_max:
            self.bounding_box_move_velocity[1] = bonuding_box_move_velocity_max * np.sign(self.bounding_box_move_velocity[1])
        if abs(self.bounding_box_move_velocity[2]) > bonuding_box_move_velocity_max:
            self.bounding_box_move_velocity[2] = bonuding_box_move_velocity_max * np.sign(self.bounding_box_move_velocity[2])

    def CalculateBoundingBoxMoveVelocityTriaxailLoading(self, measured_global_stress):

        servo_loading_stress = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingStress"].GetVector()
        servo_loading_factor = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingFactor"].GetDouble()
        d50 = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["MeanParticleDiameterD50"].GetDouble()
        Youngs_modulus = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["ParticleYoungsModulus"].GetDouble()
        
        delta_time = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)
        servo_loading_coefficient = servo_loading_factor * d50 / (delta_time * Youngs_modulus)
        stress_difference = [servo_loading_stress[0] - measured_global_stress[0][0], servo_loading_stress[1] - measured_global_stress[1][1], servo_loading_stress[2] - measured_global_stress[2][2]]
        self.bounding_box_move_velocity[0] = stress_difference[0] * servo_loading_coefficient
        self.bounding_box_move_velocity[2] = stress_difference[2] * servo_loading_coefficient

        # the loading velocity is defined by the user
        triaxial_loading_velocity = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["TriaxialLoadingVelocity"].GetDouble()
        self.bounding_box_move_velocity[1] = triaxial_loading_velocity

        bonuding_box_move_velocity_max = self.DEM_parameters["BoundingBoxServoLoadingSettings"]["BoundingBoxServoLoadingVelocityMax"].GetDouble()
        if abs(self.bounding_box_move_velocity[0]) > bonuding_box_move_velocity_max:
            self.bounding_box_move_velocity[0] = bonuding_box_move_velocity_max * np.sign(self.bounding_box_move_velocity[0])
        if abs(self.bounding_box_move_velocity[2]) > bonuding_box_move_velocity_max:
            self.bounding_box_move_velocity[2] = bonuding_box_move_velocity_max * np.sign(self.bounding_box_move_velocity[2])

    def UpdateSearchStartegyAndCPlusPlusStrategy(self, move_velocity):

        delta_time = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)

        # note: control_bool_vector = [MoveMinX, MoveMinY, MoveMinZ, MoveMaxX, MoveMaxY, MoveMaxZ]
        control_bool_vector = self.DEM_parameters["BoundingBoxMoveOptionDetail"].GetVector()
        if control_bool_vector[0]:
            self.BoundingBoxMinX_update += delta_time * move_velocity[0]
        if control_bool_vector[1]:
            self.BoundingBoxMinY_update += delta_time * move_velocity[1]
        if control_bool_vector[2]:
            self.BoundingBoxMinZ_update += delta_time * move_velocity[2]
        if control_bool_vector[3]:
            self.BoundingBoxMaxX_update -= delta_time * move_velocity[0]
        if control_bool_vector[4]:
            self.BoundingBoxMaxY_update -= delta_time * move_velocity[1]
        if control_bool_vector[5]:
            self.BoundingBoxMaxZ_update -= delta_time * move_velocity[2]

        self._GetSolver().search_strategy = OMP_DEMSearch(self.BoundingBoxMinX_update,
                                                        self.BoundingBoxMinY_update,
                                                        self.BoundingBoxMinZ_update,
                                                        self.BoundingBoxMaxX_update,
                                                        self.BoundingBoxMaxY_update,
                                                        self.BoundingBoxMaxZ_update)
        self._GetSolver().UpdateCPlusPlusStrategy()

    def UpdateIsTimeToPrintInModelParts(self, is_time_to_print):
        self.UpdateIsTimeToPrintInOneModelPart(self.spheres_model_part, is_time_to_print)
        self.UpdateIsTimeToPrintInOneModelPart(self.cluster_model_part, is_time_to_print)
        self.UpdateIsTimeToPrintInOneModelPart(self.dem_inlet_model_part, is_time_to_print)
        self.UpdateIsTimeToPrintInOneModelPart(self.rigid_face_model_part, is_time_to_print)

    def UpdateIsTimeToUpdateContactElementForServo(self, is_time_to_update_contact_element):
        self.spheres_model_part.ProcessInfo[IS_TIME_TO_UPDATE_CONTACT_ELEMENT] = is_time_to_update_contact_element
        self.cluster_model_part.ProcessInfo[IS_TIME_TO_UPDATE_CONTACT_ELEMENT] = is_time_to_update_contact_element
        self.dem_inlet_model_part.ProcessInfo[IS_TIME_TO_UPDATE_CONTACT_ELEMENT] = is_time_to_update_contact_element
        self.rigid_face_model_part.ProcessInfo[IS_TIME_TO_UPDATE_CONTACT_ELEMENT] = is_time_to_update_contact_element
        self.time_old_update_contact_element = self.time

    def UpdateIsTimeToPrintInOneModelPart(self, model_part, is_time_to_print):
        model_part.ProcessInfo[IS_TIME_TO_PRINT] = is_time_to_print

    def BeforePrintingOperations(self, time):
        pass

    def PrintAnalysisStageProgressInformation(self):
        step = self.spheres_model_part.ProcessInfo[TIME_STEPS]
        stepinfo = self.report.StepiReport(timer, self.time, step)
        if stepinfo:
            self.KratosPrintInfo(stepinfo)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        #Phantom Walls
        self.RunAnalytics(self.time)

        ##### adding DEM elements by the inlet ######
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            self.DEM_inlet.CreateElementsFromInletMesh(self.spheres_model_part, self.cluster_model_part, self.creator_destructor)  # After solving, to make sure that neighbours are already set.

    def OutputSolutionStep(self):
        #### PRINTING GRAPHS ####
        self.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.time, self.graphs_path)
        self.DEMFEMProcedures.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintBallsGraph(self.time)
        self.DEMFEMProcedures.PrintAdditionalGraphs(self.time, self._GetSolver())
        self.DEMEnergyCalculator.CalculateEnergyAndPlot(self.time)
        self.BeforePrintingOperations(self.time)
        self.PrintResults()

        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                output_process.PrintOutput()

    def BreakSolutionStepsLoop(self):
        return False

    def TheSimulationMustGoOn(self):
        it_must_or_not = self.time < self.end_time
        it_must_or_not = it_must_or_not and not self.BreakSolutionStepsLoop()
        return it_must_or_not

    def __SafeDeleteModelParts(self):
        self.model.DeleteModelPart(self.cluster_model_part.Name)
        self.model.DeleteModelPart(self.rigid_face_model_part.Name)
        self.model.DeleteModelPart(self.dem_inlet_model_part.Name)
        self.model.DeleteModelPart(self.mapping_model_part.Name)
        self.model.DeleteModelPart(self.spheres_model_part.Name)

    def Finalize(self):
        self.KratosPrintInfo("Finalizing execution...")
        super().Finalize()
        if self.do_print_results_option:
            self.GraphicalOutputFinalize()
        self.DEMFEMProcedures.FinalizeGraphs(self.rigid_face_model_part)
        self.DEMFEMProcedures.FinalizeBallsGraphs(self.spheres_model_part)
        self.DEMEnergyCalculator.FinalizeEnergyPlot()

        self.CleanUpOperations()

    def __SafeDeleteModelParts(self):
        self.model.DeleteModelPart(self.cluster_model_part.Name)
        self.model.DeleteModelPart(self.rigid_face_model_part.Name)
        self.model.DeleteModelPart(self.dem_inlet_model_part.Name)
        self.model.DeleteModelPart(self.mapping_model_part.Name)
        self.model.DeleteModelPart(self.spheres_model_part.Name)

    def CleanUpOperations(self):

        self.procedures.DeleteFiles()

        self.KratosPrintInfo(self.report.FinalReport(timer))

        if self.post_normal_impact_velocity_option:
            del self.analytic_model_part

        del self.KratosPrintInfo
        del self.all_model_parts
        if self.do_print_results_option:
            del self.demio
        del self.procedures
        del self.creator_destructor
        del self.dem_fem_search
        #del self._solver
        del self.DEMFEMProcedures
        del self.post_utils
        self.__SafeDeleteModelParts()
        del self.cluster_model_part
        del self.rigid_face_model_part
        del self.spheres_model_part
        del self.dem_inlet_model_part
        del self.mapping_model_part
        del self.contact_model_part

        if self.DEM_parameters["dem_inlet_option"].GetBool():
            del self.DEM_inlet

    def SetGraphicalOutput(self):
        self.demio = DEM_procedures.DEMIo(self.model, self.DEM_parameters, self.post_path, self.all_model_parts)
        if self.DEM_parameters["post_vtk_option"].GetBool():
            import KratosMultiphysics.DEMApplication.dem_vtk_output as dem_vtk_output
            self.vtk_output = dem_vtk_output.VtkOutput(self.main_path, self.problem_name, self.spheres_model_part, self.contact_model_part, self.rigid_face_model_part, self.DEM_parameters)

    def GraphicalOutputInitialize(self):
        if self.do_print_results_option:
            self.demio.Initialize(self.DEM_parameters)
            self.demio.InitializeMesh(self.all_model_parts)

    def PrintResultsForGid(self, time):
        if self._GetSolver().poisson_ratio_option:
            self.DEMFEMProcedures.PrintPoisson(self.spheres_model_part, self.DEM_parameters, "Poisson_ratio.txt", time)

        if self.DEM_parameters["PostEulerAngles"].GetBool():
            self.post_utils.PrintEulerAngles(self.spheres_model_part, self.cluster_model_part)

        self.demio.PrintMultifileLists(time, self.post_path)
        self._GetSolver().PrepareElementsForPrinting()
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self._GetSolver().PrepareContactElementsForPrinting()

        if "post_gid_option" in self.DEM_parameters.keys():
            if self.DEM_parameters["post_gid_option"].GetBool() != False:
                self.demio.ShowPrintingResultsOnScreen(self.all_model_parts, 'GID')
                self.demio.PrintResults(self.all_model_parts, self.creator_destructor, self.dem_fem_search, time, self.bounding_box_time_limits)

        if "post_vtk_option" in self.DEM_parameters.keys():
            if self.DEM_parameters["post_vtk_option"].GetBool():
                self.demio.ShowPrintingResultsOnScreen(self.all_model_parts, 'VTK')
                self.vtk_output.WriteResults(self.time)

        self.file_msh = self.demio.GetMultiFileListName(self.problem_name + "_" + "%.12g"%time + ".post.msh")
        self.file_res = self.demio.GetMultiFileListName(self.problem_name + "_" + "%.12g"%time + ".post.res")

    def GraphicalOutputFinalize(self):
        self.demio.FinalizeMesh()
        self.demio.CloseMultifiles()

    #these functions are needed for coupling, so that single time loops can be done

    def InitializeTime(self):
        self.time = 0.0
        self.time_old_print = 0.0

    def FinalizeSingleTimeStep(self):
        message = 'Warning!'
        message += '\nFunction \'FinalizeSingleTimeStep\' is deprecated. Use FinalizeSolutionStep instead.'
        message += '\nIt will be removed after 10/31/2019.\n'
        Logger.PrintWarning("DEM_analysis_stage.py", message)
        ##### adding DEM elements by the inlet ######
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            self.DEM_inlet.CreateElementsFromInletMesh(self.spheres_model_part, self.cluster_model_part, self.creator_destructor)  # After solving, to make sure that neighbours are already set.
        stepinfo = self.report.StepiReport(timer, self.time, self.step)
        if stepinfo:
            self.KratosPrintInfo(stepinfo)

    def OutputSingleTimeLoop(self):
        #### PRINTING GRAPHS ####
        os.chdir(self.graphs_path)
        self.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.time, self.graphs_path)
        self.DEMFEMProcedures.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintBallsGraph(self.time)
        self.DEMEnergyCalculator.CalculateEnergyAndPlot(self.time)
        self.BeforePrintingOperations(self.time)
        #### GiD IO ##########################################
        time_to_print = self.time - self.time_old_print
        if self.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self._GetSolver().dt:
            self.PrintResultsForGid(self.time)
            self.time_old_print = self.time

    def MeasureSphereForGettingPackingProperties(self, radius, center_x, center_y, center_z, type):        
        domain_size = [0.0, 0.0, 0.0]
        domain_size[0] = self.BoundingBoxMaxX_update - self.BoundingBoxMinX_update
        domain_size[1] = self.BoundingBoxMaxY_update - self.BoundingBoxMinY_update
        domain_size[2] = self.BoundingBoxMaxZ_update - self.BoundingBoxMinZ_update
        return self.DEMPropertiesMeasureUtility.MeasureSphereForGettingPackingProperties(radius, center_x, center_y, center_z, type, domain_size)
    
    def MeasureSphereForGettingGlobalStressTensor(self):
        Lx = self.BoundingBoxMaxX_update - self.BoundingBoxMinX_update
        Ly = self.BoundingBoxMaxY_update - self.BoundingBoxMinY_update
        Lz = self.BoundingBoxMaxZ_update - self.BoundingBoxMinZ_update
        return self.DEMPropertiesMeasureUtility.MeasureSphereForGettingGlobalStressTensor(Lx, Ly, Lz)
    
    def MeasureSphereForGettingRadialDistributionFunction(self, radius, center_x, center_y, center_z, delta_r, d_mean):
        self.DEMPropertiesMeasureUtility.MeasureSphereForGettingRadialDistributionFunction(radius, center_x, center_y, center_z, delta_r, d_mean)

    def MeasureCubicForGettingPackingProperties(self, side_length, center_x, center_y, center_z, type):
        domain_size = [0.0, 0.0, 0.0]
        domain_size[0] = self.BoundingBoxMaxX_update - self.BoundingBoxMinX_update
        domain_size[1] = self.BoundingBoxMaxY_update - self.BoundingBoxMinY_update
        domain_size[2] = self.BoundingBoxMaxZ_update - self.BoundingBoxMinZ_update
        return self.DEMPropertiesMeasureUtility.MeasureCubicForGettingPackingProperties(side_length, center_x, center_y, center_z, type, domain_size)
   
        '''
        if type == "stress_tensor":
            
            if self.DEM_parameters["PostStressStrainOption"].GetBool() and self.DEM_parameters["ContactMeshOption"].GetBool():
                
                measure_cubic_volume = side_length ** 3
                total_tensor        = np.empty((3, 3))
                total_tensor[:]     = 0.0

                Lx = self.BoundingBoxMaxX_update - self.BoundingBoxMinX_update
                Ly = self.BoundingBoxMaxY_update - self.BoundingBoxMinY_update
                Lz = self.BoundingBoxMaxZ_update - self.BoundingBoxMinZ_update
                
                for element in self.contact_model_part.Elements:
            
                    x_0 = element.GetNode(0).X
                    x_1 = element.GetNode(1).X
                    y_0 = element.GetNode(0).Y
                    y_1 = element.GetNode(1).Y
                    z_0 = element.GetNode(0).Z
                    z_1 = element.GetNode(1).Z
                    
                    x_mid = 0.5 * (x_0 + x_1)
                    y_mid = 0.5 * (y_0 + y_1)
                    z_mid = 0.5 * (z_0 + z_1)

                    sphere_1_is_inside = False
                    sphere_2_is_inside = False
                    contact_is_inside = False

                    if (center_x - 0.5 * side_length) < x_0 < (center_x + 0.5 * side_length):
                        if (center_y - 0.5 * side_length) < y_0 < (center_y + 0.5 * side_length):
                            if (center_z - 0.5 * side_length) < z_0 < (center_z + 0.5 * side_length):
                                sphere_1_is_inside = True

                    if (center_x - 0.5 * side_length) < x_1 < (center_x + 0.5 * side_length):
                        if (center_y - 0.5 * side_length) < y_1 < (center_y + 0.5 * side_length):
                            if (center_z - 0.5 * side_length) < z_1 < (center_z + 0.5 * side_length):
                                sphere_2_is_inside = True
                    
                    if (center_x - 0.5 * side_length) < x_mid < (center_x + 0.5 * side_length):
                        if (center_y - 0.5 * side_length) < y_mid < (center_y + 0.5 * side_length):
                            if (center_z - 0.5 * side_length) < z_mid < (center_z + 0.5 * side_length):
                                contact_is_inside = True

                    if (sphere_1_is_inside or sphere_2_is_inside) and contact_is_inside:
                
                        dx = x_0 - x_1
                        if dx > 0.5 * Lx:
                            dx -= Lx
                        elif dx < -0.5 * Lx:
                            dx += Lx
                        
                        dy = y_0 - y_1
                        if dy > 0.5 * Ly:
                            dy -= Ly
                        elif dy < -0.5 * Ly:
                            dy += Ly

                        dz = z_0 - z_1
                        if dz > 0.5 * Lz:
                            dz -= Lz
                        elif dz < -0.5 * Lz:
                            dz += Lz
                        
                        local_contact_force_X = element.GetValue(GLOBAL_CONTACT_FORCE)[0]
                        local_contact_force_Y = element.GetValue(GLOBAL_CONTACT_FORCE)[1]
                        local_contact_force_Z = element.GetValue(GLOBAL_CONTACT_FORCE)[2]
                        contact_force_vector = np.array([local_contact_force_X , local_contact_force_Y, local_contact_force_Z])
                        vector_l = np.array([dx , dy, dz])
                        tensor = np.outer(contact_force_vector, vector_l)
                        total_tensor += tensor

                total_tensor = total_tensor / measure_cubic_volume

                return total_tensor
            
            else:
                
                raise Exception('The \"PostStressStrainOption\" and \"ContactMeshOption\" in the [ProjectParametersDEM.json] should be [True].')
            '''
if __name__ == "__main__":
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    DEMAnalysisStage(model, project_parameters).Run()
