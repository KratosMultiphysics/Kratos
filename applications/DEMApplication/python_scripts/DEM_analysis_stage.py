from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.DEMApplication.DEM_restart_utility import DEMRestartUtility

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

    @classmethod
    def GetDefaultInputParameters(self):
        import KratosMultiphysics.DEMApplication.dem_default_input_parameters
        return KratosMultiphysics.DEMApplication.dem_default_input_parameters.GetDefaultInputParameters()

    @classmethod
    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid, IO.SKIP_TIMER)

    @classmethod
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
        self.creator_destructor = self.SetParticleCreatorDestructor()
        self.dem_fem_search = self.SetDemFemSearch()
        self.procedures = self.SetProcedures()
        self.SetAnalyticParticleWatcher()
        self.PreUtilities = PreUtilities()
        self.aux = AuxiliaryUtilities()

        # Set the print function TO_DO: do this better...
        self.KratosPrintInfo = self.procedures.KratosPrintInfo

        # Creating necessary directories:
        self.problem_name = self.GetProblemTypeFileName()

        [self.post_path,
        self.data_and_results,
        self.graphs_path] = self.procedures.CreateDirectories(str(self.main_path),
                                                              str(self.problem_name),
                                                              do_print_results=self.do_print_results_option)[:-1]

        # Prepare modelparts
        self.CreateModelParts()

        if self.do_print_results_option:
            self.SetGraphicalOutput()
        self.report = DEM_procedures.Report()
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest = DEM_procedures.MaterialTest()
        self.translational_scheme = self.SetTranslationalScheme()
        self.rotational_scheme = self.SetRotationalScheme()

        # Define control variables
        self.p_frequency = 100   # activate every 100 steps
        self.step_count = 0
        self.p_count = self.p_frequency

        #self._solver = self._GetSolver()
        self.SetDt()
        self.SetFinalTime()
        self.AddVariables()
        super(DEMAnalysisStage, self).__init__(model, self.DEM_parameters)

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

    def SetAnalyticParticleWatcher(self):
        from KratosMultiphysics.DEMApplication.analytic_tools import analytic_data_procedures
        self.particle_watcher = AnalyticParticleWatcher()

        # is this being used? TODO
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(analytic_particle_watcher=self.particle_watcher, path=self.main_path)


    def SetAnalyticFaceWatcher(self):
        from KratosMultiphysics.DEMApplication.analytic_tools import analytic_data_procedures
        self.FaceAnalyzerClass = analytic_data_procedures.FaceWatcherAnalyzer
        self.face_watcher_dict = dict()
        self.face_watcher_analysers = dict()
        for sub_part in self.rigid_face_model_part.SubModelParts:
            if sub_part[IS_GHOST] == True:
                name = sub_part.Name
                self.face_watcher_dict[sub_part.Name] = AnalyticFaceWatcher(sub_part)
                self.face_watcher_analysers[sub_part.Name] = analytic_data_procedures.FaceWatcherAnalyzer(name=name, analytic_face_watcher=self.face_watcher_dict[sub_part.Name], path=self.main_path)

    def MakeAnalyticsMeasurements(self):
        for face_watcher in self.face_watcher_dict.values():
            face_watcher.MakeMeasurements()

    def SetFinalTime(self):
        self.end_time = self.DEM_parameters["FinalTime"].GetDouble()

    def SetDt(self):
        self._GetSolver().dt = self.DEM_parameters["MaxTimeStep"].GetDouble()

    def SetProcedures(self):
        return DEM_procedures.Procedures(self.DEM_parameters)

    @classmethod
    def SetDemFemSearch(self):
        return DEM_FEM_Search()

    @classmethod
    def GetParticleHistoryWatcher(self):
        return None

    def SetParticleCreatorDestructor(self):

        self.watcher = self.GetParticleHistoryWatcher()

        if self.watcher is None:
            return ParticleCreatorDestructor()
        return ParticleCreatorDestructor(self.watcher)

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
            self.KratosPrintWarning('Error: selected translational integration scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return translational_scheme

    def SetRotationalScheme(self):
        rotational_scheme = self.SelectRotationalScheme()

        if rotational_scheme is None:
            self.KratosPrintWarning('Error: selected rotational integration scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
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

    def FillAnalyticSubModelPartsWithNewParticles(self):
        self.analytic_model_part = self.spheres_model_part.GetSubModelPart('AnalyticParticlesPart')
        self.PreUtilities.FillAnalyticSubModelPartUtility(self.spheres_model_part, self.analytic_model_part)
        #analytic_particle_ids = [elem.Id for elem in self.spheres_model_part.Elements]
        #self.analytic_model_part.AddElements(analytic_particle_ids)

    def Initialize(self):
        self.time = 0.0
        self.time_old_print = 0.0

        self.ReadModelParts()

        self.SetAnalyticFaceWatcher()

        self.post_normal_impact_velocity_option = False
        if "PostNormalImpactVelocity" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostNormalImpactVelocity"].GetBool():
                self.post_normal_impact_velocity_option = True
                self.FillAnalyticSubModelParts()

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

        #Finding the max id of the nodes... (it is necessary for anything that will add spheres to the self.spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.z
        #max_Id = self.procedures.FindMaxNodeIdAccrossModelParts(self.creator_destructor, self.all_model_parts)   # TODO this seems not be longer required
        #self.creator_destructor.SetMaxNodeId(max_Id)
        self.creator_destructor.SetMaxNodeId(self.all_model_parts.MaxNodeId)  #TODO check functionalities

        self.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(self.DEM_parameters, self.graphs_path, self.spheres_model_part, self.rigid_face_model_part)

        self.DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(self.DEM_parameters, self.spheres_model_part, self.cluster_model_part, self.graphs_path, "EnergyPlot.grf")

        self.materialTest.Initialize(self.DEM_parameters, self.procedures, self._GetSolver(), self.graphs_path, self.post_path, self.spheres_model_part, self.rigid_face_model_part)

        self.KratosPrintInfo("Initialization Complete")

        self.report.Prepare(timer, self.DEM_parameters["ControlTime"].GetDouble())

        self.materialTest.PrintChart()
        self.materialTest.PrepareDataForGraph()

        self.post_utils = DEM_procedures.PostUtils(self.DEM_parameters, self.spheres_model_part)
        self.report.total_steps_expected = int(self.end_time / self._GetSolver().dt)

        super(DEMAnalysisStage, self).Initialize()

        #Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation
        #Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
        self.SetInlet()

        self.SetInitialNodalValues()

        self.KratosPrintInfo(self.report.BeginReport(timer))

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

    def RunAnalytics(self, time, is_time_to_print=True):
        for sp in (sp for sp in self.rigid_face_model_part.SubModelParts if sp[IS_GHOST]):
            self.MakeAnalyticsMeasurements()
            if is_time_to_print:
                self.FaceAnalyzerClass.CreateNewFile()
                for sp in (sp for sp in self.rigid_face_model_part.SubModelParts if sp[IS_GHOST]):
                    self.face_watcher_analysers[sp.Name].UpdateDataFiles(time)
                self.FaceAnalyzerClass.RemoveOldFile()

    def IsTimeToPrintPostProcess(self):
        return self.do_print_results_option and self.DEM_parameters["OutputTimeStep"].GetDouble() - (self.time - self.time_old_print) < 1e-2 * self._GetSolver().dt

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
            self.DEM_inlet = DEM_Inlet(self.dem_inlet_model_part)
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor, self._GetSolver().continuum_type)

    def SetInitialNodalValues(self):
        self.procedures.SetInitialNodalValues(self.spheres_model_part, self.cluster_model_part, self.dem_inlet_model_part, self.rigid_face_model_part)

    def InitializeSolutionStep(self):
        super(DEMAnalysisStage, self).InitializeSolutionStep()
        if self.post_normal_impact_velocity_option:
            if self.IsCountStep():
                self.FillAnalyticSubModelPartsWithNewParticles()
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self.UpdateIsTimeToPrintInModelParts(self.IsTimeToPrintPostProcess())

        if self.DEM_parameters["Dimension"].GetInt() == 2:
            self.spheres_model_part.ProcessInfo[IMPOSED_Z_STRAIN_OPTION] = self.DEM_parameters["ImposeZStrainIn2DOption"].GetBool()
            if not self.DEM_parameters["ImposeZStrainIn2DWithControlModule"].GetBool():
                if self.spheres_model_part.ProcessInfo[IMPOSED_Z_STRAIN_OPTION]:
                    t = self.time
                    self.spheres_model_part.ProcessInfo.SetValue(IMPOSED_Z_STRAIN_VALUE, eval(self.DEM_parameters["ZStrainValue"].GetString()))


    def UpdateIsTimeToPrintInModelParts(self, is_time_to_print):
        self.UpdateIsTimeToPrintInOneModelPart(self.spheres_model_part, is_time_to_print)
        self.UpdateIsTimeToPrintInOneModelPart(self.cluster_model_part, is_time_to_print)
        self.UpdateIsTimeToPrintInOneModelPart(self.dem_inlet_model_part, is_time_to_print)
        self.UpdateIsTimeToPrintInOneModelPart(self.rigid_face_model_part, is_time_to_print)

    @classmethod
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
        super(DEMAnalysisStage, self).FinalizeSolutionStep()
        if self.post_normal_impact_velocity_option:
            self.particle_watcher.MakeMeasurements(self.analytic_model_part)
            if self.IsTimeToPrintPostProcess():
                self.particle_watcher.SetNodalMaxImpactVelocities(self.analytic_model_part)
                self.particle_watcher.SetNodalMaxFaceImpactVelocities(self.analytic_model_part)

        #Phantom Walls
        self.RunAnalytics(self.time, self.IsTimeToPrintPostProcess())

        ##### adding DEM elements by the inlet ######
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            self.DEM_inlet.CreateElementsFromInletMesh(self.spheres_model_part, self.cluster_model_part, self.creator_destructor)  # After solving, to make sure that neighbours are already set.

    def OutputSolutionStep(self):
        #### PRINTING GRAPHS ####
        self.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.time, self.graphs_path)
        self.materialTest.MeasureForcesAndPressure()
        self.materialTest.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintBallsGraph(self.time)
        self.DEMEnergyCalculator.CalculateEnergyAndPlot(self.time)
        self.BeforePrintingOperations(self.time)
        self.PrintResults()

        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                output_process.PrintOutput()

        self.FinalizeTimeStep(self.time)

    def AfterSolveOperations(self):
        message = 'Warning!'
        message += '\nFunction \'AfterSolveOperations\' is deprecated.'
        message += '\nIt will be removed after 10/31/2019.\n'
        Logger.PrintWarning("DEM_analysis_stage.py", message)
        if self.post_normal_impact_velocity_option:
            self.particle_watcher.MakeMeasurements(self.analytic_model_part)
            if self.IsTimeToPrintPostProcess():
                self.particle_watcher.SetNodalMaxImpactVelocities(self.analytic_model_part)
                self.particle_watcher.SetNodalMaxFaceImpactVelocities(self.analytic_model_part)

        #Phantom Walls
        self.RunAnalytics(self.time, self.IsTimeToPrintPostProcess())

    def FinalizeTimeStep(self, time):
        pass

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
        if self.do_print_results_option:
            self.GraphicalOutputFinalize()
        self.materialTest.FinalizeGraphs()
        self.DEMFEMProcedures.FinalizeGraphs(self.rigid_face_model_part)
        self.DEMFEMProcedures.FinalizeBallsGraphs(self.spheres_model_part)
        self.DEMEnergyCalculator.FinalizeEnergyPlot()

        self.AdditionalFinalizeOperations()

        self.CleanUpOperations()

    def AdditionalFinalizeOperations(self):
        pass

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

        if self.DEM_parameters["dem_inlet_option"].GetBool():
            del self.DEM_inlet

    def SetGraphicalOutput(self):
        self.demio = DEM_procedures.DEMIo(self.model, self.DEM_parameters, self.post_path, self.all_model_parts)
        if self.DEM_parameters["post_vtk_option"].GetBool():
            import KratosMultiphysics.DEMApplication.dem_vtk_output as dem_vtk_output
            self.vtk_output = dem_vtk_output.VtkOutput(self.main_path, self.problem_name, self.spheres_model_part, self.rigid_face_model_part)

    def GraphicalOutputInitialize(self):
        if self.do_print_results_option:
            self.demio.Initialize(self.DEM_parameters)
            self.demio.InitializeMesh(self.all_model_parts)

    def PrintResultsForGid(self, time):
        if self._GetSolver().poisson_ratio_option:
            self.DEMFEMProcedures.PrintPoisson(self.spheres_model_part, self.DEM_parameters, "Poisson_ratio.txt", time)

        if self.DEM_parameters["PostEulerAngles"].GetBool():
            self.post_utils.PrintEulerAngles(self.spheres_model_part, self.cluster_model_part)

        self.demio.ShowPrintingResultsOnScreen(self.all_model_parts)

        self.demio.PrintMultifileLists(time, self.post_path)
        self._GetSolver().PrepareElementsForPrinting()
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self._GetSolver().PrepareContactElementsForPrinting()

        self.demio.PrintResults(self.all_model_parts, self.creator_destructor, self.dem_fem_search, time, self.bounding_box_time_limits)
        if "post_vtk_option" in self.DEM_parameters.keys():
            if self.DEM_parameters["post_vtk_option"].GetBool():
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
        self.materialTest.MeasureForcesAndPressure()
        self.materialTest.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintBallsGraph(self.time)
        self.DEMEnergyCalculator.CalculateEnergyAndPlot(self.time)
        self.BeforePrintingOperations(self.time)
        #### GiD IO ##########################################
        time_to_print = self.time - self.time_old_print
        if self.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self._GetSolver().dt:
            self.PrintResultsForGid(self.time)
            self.time_old_print = self.time
        self.FinalizeTimeStep(self.time)


if __name__ == "__main__":
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    DEMAnalysisStage(model, project_parameters).Run()
