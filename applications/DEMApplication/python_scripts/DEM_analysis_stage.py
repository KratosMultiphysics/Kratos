from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
sys.path.insert(0, '')
from analysis_stage import AnalysisStage

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ or "I_MPI_INFO_NUMA_NODE_NUM" in os.environ:
    if "DO_NOT_PARTITION_DOMAIN" in os.environ:
        Logger.PrintInfo("DEM", "Running under MPI........")
        from KratosMultiphysics.mpi import *
        import DEM_procedures_mpi_no_partitions as DEM_procedures
        import DEM_material_test_script
    else:
        Logger.PrintInfo("DEM", "Running under OpenMP........")
        from KratosMultiphysics.MetisApplication import *
        from KratosMultiphysics.MPISearchApplication import *
        from KratosMultiphysics.mpi import *
        import DEM_procedures_mpi as DEM_procedures
        import DEM_material_test_script_mpi as DEM_material_test_script
else:
    Logger.PrintInfo("DEM", "Running under OpenMP........")
    import DEM_procedures
    import DEM_material_test_script

class DEMAnalysisStage(AnalysisStage):

    def GetParametersFileName(self):
        return "ProjectParametersDEM.json"

    def GetInputParameters(self):
        parameters_file_name = self.GetParametersFileName()
        parameters_file = open(parameters_file_name, 'r')
        return Parameters(parameters_file.read())

    def LoadParametersFile(self):
        self.DEM_parameters = self.GetInputParameters()
        self.project_parameters = self.GetInputParameters()
        default_input_parameters = self.GetDefaultInputParameters()
        self.DEM_parameters.ValidateAndAssignDefaults(default_input_parameters)

    @classmethod
    def GetDefaultInputParameters(self):
        import dem_default_input_parameters
        return dem_default_input_parameters.GetDefaultInputParameters()

    @classmethod
    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid, IO.SKIP_TIMER)

    @classmethod
    def GetMainPath(self):
        return os.getcwd()

    def __init__(self, model, parameters):
        self.model = model
        self.main_path = self.GetMainPath()
        self.LoadParametersFile()
        self.solver_strategy = self.SetSolverStrategy()
        self.creator_destructor = self.SetParticleCreatorDestructor()
        self.dem_fem_search = self.SetDemFemSearch()
        self.procedures = self.SetProcedures()
        self.SetAnalyticParticleWatcher()
        self.PreUtilities = PreUtilities()
        self.aux = AuxiliaryUtilities()

        # Set the print function TO_DO: do this better...
        self.KRATOSprint = self.procedures.KRATOSprint

        # Creating necessary directories:
        self.problem_name = self.GetProblemTypeFilename()

        [self.post_path,
        self.data_and_results,
        self.graphs_path,
        MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(self.problem_name))

        # Prepare modelparts
        self.CreateModelParts()

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

        self.solver = self.SetSolver()
        self.SetDt()
        self.SetFinalTime()
        super(DEMAnalysisStage, self).__init__(model, self.DEM_parameters)

    def CreateModelParts(self):
        self.spheres_model_part = self.model.CreateModelPart("SpheresPart")
        self.rigid_face_model_part = self.model.CreateModelPart("RigidFacePart")
        self.cluster_model_part = self.model.CreateModelPart("ClusterPart")
        self.DEM_inlet_model_part = self.model.CreateModelPart("DEMInletPart")
        self.mapping_model_part = self.model.CreateModelPart("MappingPart")
        self.contact_model_part = self.model.CreateModelPart("ContactPart")

        mp_list = []
        mp_list.append(self.spheres_model_part)
        mp_list.append(self.rigid_face_model_part)
        mp_list.append(self.cluster_model_part)
        mp_list.append(self.DEM_inlet_model_part)
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
        from analytic_tools import analytic_data_procedures
        self.particle_watcher = AnalyticParticleWatcher()

        # is this being used? TODO
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(analytic_particle_watcher=self.particle_watcher, path=self.main_path)


    def SetAnalyticFaceWatcher(self):
        from analytic_tools import analytic_data_procedures
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
        self.solver.dt = self.DEM_parameters["MaxTimeStep"].GetDouble()

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
            self.KRATOSprint('Error: selected translational integration scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return translational_scheme

    def SetRotationalScheme(self):
        rotational_scheme = self.SelectRotationalScheme()

        if rotational_scheme is None:
            self.KRATOSprint('Error: selected rotational integration scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return rotational_scheme

    def SetSolverStrategy(self):

        # TODO: Ugly fix. Change it. I don't like this to be in the main...
        # Strategy object
        if self.DEM_parameters["ElementType"].GetString() == "SphericPartDEMElement3D" or self.DEM_parameters["ElementType"].GetString() == "CylinderPartDEMElement2D":
            import sphere_strategy as SolverStrategy
        elif self.DEM_parameters["ElementType"].GetString() == "SphericContPartDEMElement3D" or self.DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D":
            import continuum_sphere_strategy as SolverStrategy
        elif self.DEM_parameters["ElementType"].GetString() == "ThermalSphericContPartDEMElement3D":
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif self.DEM_parameters["ElementType"].GetString() == "ThermalSphericPartDEMElement3D":
            import thermal_sphere_strategy as SolverStrategy
        elif self.DEM_parameters["ElementType"].GetString() == "SinteringSphericConPartDEMElement3D":
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif self.DEM_parameters["ElementType"].GetString() == "IceContPartDEMElement3D":
            import ice_continuum_sphere_strategy as SolverStrategy
        else:
            self.KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')

        return SolverStrategy

    def SetSolver(self):
        return self.solver_strategy.ExplicitStrategy(self.all_model_parts,
                                                     self.creator_destructor,
                                                     self.dem_fem_search,
                                                     self.DEM_parameters,
                                                     self.procedures)


    def AddVariables(self):
        self.procedures.AddAllVariablesInAllModelParts(self.solver, self.translational_scheme, self.rotational_scheme, self.all_model_parts, self.DEM_parameters)

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
        self.step = 0
        self.time = 0.0

        self.AddVariables()

        self.ReadModelParts()

        self.SetAnalyticFaceWatcher()  # TODO check order

        self.post_normal_impact_velocity_option = False
        if "PostNormalImpactVelocity" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostNormalImpactVelocity"].GetBool():
                self.post_normal_impact_velocity_option = True
                self.FillAnalyticSubModelParts()

        # Setting up the buffer size
        self.procedures.SetUpBufferSizeInAllModelParts(self.spheres_model_part, 1, self.cluster_model_part, 1, self.DEM_inlet_model_part, 1, self.rigid_face_model_part, 1)
        # Adding dofs
        self.AddAllDofs()

        #-----------os.chdir(self.main_path)
        self.KRATOSprint("Initializing Problem...")

        self.GraphicalOutputInitialize()

        # Perform a partition to balance the problem
        self.SetSearchStrategy()

        self.SolverBeforeInitialize()

        self.parallelutils.Repart(self.spheres_model_part)

        #Setting up the BoundingBox
        self.bounding_box_time_limits = self.procedures.SetBoundingBoxLimits(self.all_model_parts, self.creator_destructor)

        #Finding the max id of the nodes... (it is necessary for anything that will add spheres to the self.spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.z
        max_Id = self.procedures.FindMaxNodeIdAccrossModelParts(self.creator_destructor, self.all_model_parts)
        #self.creator_destructor.SetMaxNodeId(max_Id)
        self.creator_destructor.SetMaxNodeId(self.all_model_parts.MaxNodeId)  #TODO check functionalities

        #Strategy Initialization
        #-------------os.chdir(self.main_path)

        self.SolverInitialize()

        #Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation
        #Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
        self.SetInlet()

        self.SetInitialNodalValues()

        self.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(self.DEM_parameters, self.graphs_path, self.spheres_model_part, self.rigid_face_model_part)

        #------------os.chdir(self.graphs_path)
        self.DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(self.DEM_parameters, self.spheres_model_part, self.cluster_model_part, self.graphs_path, "EnergyPlot.grf")

        self.materialTest.Initialize(self.DEM_parameters, self.procedures, self.solver, self.graphs_path, self.post_path, self.spheres_model_part, self.rigid_face_model_part)

        self.KRATOSprint("Initialization Complete")

        self.report.Prepare(timer, self.DEM_parameters["ControlTime"].GetDouble())

        #self.procedures.ModelData(self.spheres_model_part, self.solver) #check link with ModelDataInfo = "OFF"

        self.materialTest.PrintChart()
        self.materialTest.PrepareDataForGraph()

        self.post_utils = DEM_procedures.PostUtils(self.DEM_parameters, self.spheres_model_part)
        self.report.total_steps_expected = int(self.end_time / self.solver.dt)
        self.KRATOSprint(self.report.BeginReport(timer))
        #-----os.chdir(self.main_path)

    def AddAllDofs(self):
        self.solver.AddDofs(self.spheres_model_part)
        self.solver.AddDofs(self.cluster_model_part)
        self.solver.AddDofs(self.DEM_inlet_model_part)

    def SetSearchStrategy(self):
        self.solver.search_strategy = self.parallelutils.GetSearchStrategy(self.solver, self.spheres_model_part)

    def SolverBeforeInitialize(self):
        self.solver.BeforeInitialize()

    def SolverInitialize(self):
        self.solver.Initialize() # Possible modifications of number of elements and number of nodes

    def GetProblemNameWithPath(self):
        return self.DEM_parameters["problem_name"].GetString()

    def GetMpFilename(self):
        return self.GetProblemNameWithPath() + "DEM"

    def GetInletFilename(self):
        return self.GetProblemNameWithPath() + "DEM_Inlet"

    def GetFemFilename(self):
        return self.GetProblemNameWithPath() + "DEM_FEM_boundary"

    def GetClusterFilename(self):
        return self.GetProblemNameWithPath() + "DEM_Clusters"

    def GetProblemTypeFilename(self):
        return self.DEM_parameters["problem_name"].GetString()

    def ReadModelParts(self, max_node_Id=0, max_elem_Id=0, max_cond_Id=0):
        #-----os.chdir(self.main_path)

        # Reading the model_part
        spheres_mp_filename = self.GetMpFilename()
        model_part_io_spheres = self.model_part_reader(spheres_mp_filename, max_node_Id, max_elem_Id, max_cond_Id)

        if "do_not_perform_initial_partition" in self.DEM_parameters.keys() and self.DEM_parameters["do_not_perform_initial_partition"].GetBool():
            pass
        else:
            self.parallelutils.PerformInitialPartition(model_part_io_spheres)

        #-----os.chdir(self.main_path)
        [model_part_io_spheres, self.spheres_model_part, MPICommSetup] = self.parallelutils.SetCommunicator(self.spheres_model_part, model_part_io_spheres, spheres_mp_filename)
        model_part_io_spheres.ReadModelPart(self.spheres_model_part)

        max_node_Id = max(max_node_Id, self.creator_destructor.FindMaxNodeIdInModelPart(self.spheres_model_part))
        max_elem_Id = max(max_elem_Id, self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part))
        old_max_elem_Id_spheres = max_elem_Id
        max_cond_Id = max(max_cond_Id, self.creator_destructor.FindMaxConditionIdInModelPart(self.spheres_model_part))
        rigidFace_mp_filename = self.GetFemFilename()
        model_part_io_fem = self.model_part_reader(rigidFace_mp_filename, max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)
        model_part_io_fem.ReadModelPart(self.rigid_face_model_part)

        max_node_Id = max(max_node_Id, self.creator_destructor.FindMaxNodeIdInModelPart(self.rigid_face_model_part))
        max_elem_Id = max(max_elem_Id, self.creator_destructor.FindMaxElementIdInModelPart(self.rigid_face_model_part))
        max_cond_Id = max(max_cond_Id, self.creator_destructor.FindMaxConditionIdInModelPart(self.rigid_face_model_part))
        clusters_mp_filename = self.GetClusterFilename()
        model_part_io_clusters = self.model_part_reader(clusters_mp_filename, max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)
        model_part_io_clusters.ReadModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        if max_elem_Id != old_max_elem_Id_spheres:
            self.creator_destructor.RenumberElementIdsFromGivenValue(self.cluster_model_part, max_elem_Id)

        max_node_Id = max(max_node_Id, self.creator_destructor.FindMaxNodeIdInModelPart(self.cluster_model_part))
        max_elem_Id = max(max_elem_Id, self.creator_destructor.FindMaxElementIdInModelPart(self.cluster_model_part))
        max_cond_Id = max(max_cond_Id, self.creator_destructor.FindMaxConditionIdInModelPart(self.cluster_model_part))
        DEM_Inlet_filename = self.GetInletFilename()
        model_part_io_demInlet = self.model_part_reader(DEM_Inlet_filename, max_node_Id + 1, max_elem_Id + 1, max_cond_Id + 1)
        model_part_io_demInlet.ReadModelPart(self.DEM_inlet_model_part)

        self.model_parts_have_been_read = True
        self.all_model_parts.ComputeMaxIds()

    def RunMainTemporalLoop(self): # deprecated
        self.RunSolutionLoop()

    def RunSolutionLoop(self):

        self.step = 0
        self.time = 0.0
        self.time_old_print = 0.0

        while self.time < self.end_time:
            self.step, self.time = self._GetSolver().AdvanceInTime(self.step, self.time)

            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            self._GetSolver().SolveSolutionStep()
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
            if self.BreakSolutionStepsLoop():
                break

    def RunAnalytics(self, time, is_time_to_print=True):
        for sp in (sp for sp in self.rigid_face_model_part.SubModelParts if sp[IS_GHOST]):
            self.MakeAnalyticsMeasurements()
            if is_time_to_print:
                self.FaceAnalyzerClass.CreateNewFile()
                for sp in (sp for sp in self.rigid_face_model_part.SubModelParts if sp[IS_GHOST]):
                    self.face_watcher_analysers[sp.Name].UpdateDataFiles(time)
                self.FaceAnalyzerClass.RemoveOldFile()

    def IsTimeToPrintPostProcess(self, time):
        return self.DEM_parameters["OutputTimeStep"].GetDouble() - (time - self.time_old_print) < 1e-2 * self.solver.dt

    def PrintResults(self):
        #### GiD IO ##########################################
        if self.IsTimeToPrintPostProcess(self.time):
            self.PrintResultsForGid(self.time)
            self.time_old_print = self.time


    def UpdateTimeInModelParts(self):
        self.DEMFEMProcedures.UpdateTimeInModelParts(self.all_model_parts, self.time, self.solver.dt, self.step, self.IsTimeToPrintPostProcess(self.time))

    def UpdateTimeInOneModelPart(self):
        pass

    def SolverSolve(self):
        self.solver.SolveSolutionStep()

    def _GetSolver(self):
        return self.solver

    def SetInlet(self):
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            #Constructing the inlet and initializing it (must be done AFTER the self.spheres_model_part Initialize)
            self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor, self.solver.continuum_type)

    def SetInitialNodalValues(self):
        self.procedures.SetInitialNodalValues(self.spheres_model_part, self.cluster_model_part, self.DEM_inlet_model_part, self.rigid_face_model_part)

    def InitializeTimeStep(self): # deprecated
        self.InitializeSolutionStep()

    def InitializeSolutionStep(self):

        self.UpdateTimeInModelParts()

        self.BeforeSolveOperations(self.time)

    def BeforeSolveOperations(self, time):
        if self.post_normal_impact_velocity_option:
            if self.IsCountStep():
                self.FillAnalyticSubModelPartsWithNewParticles()

    def BeforePrintingOperations(self, time):
        pass

    def FinalizeSolutionStep(self):
        super(DEMAnalysisStage, self).FinalizeSolutionStep()
        self.AfterSolveOperations()

        self.DEMFEMProcedures.MoveAllMeshes(self.all_model_parts, self.time, self.solver.dt)

        ##### adding DEM elements by the inlet ######
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            self.DEM_inlet.CreateElementsFromInletMesh(self.spheres_model_part, self.cluster_model_part, self.creator_destructor)  # After solving, to make sure that neighbours are already set.

        stepinfo = self.report.StepiReport(timer, self.time, self.step)
        if stepinfo:
            self.KRATOSprint(stepinfo)

    def OutputSolutionStep(self):
        #### PRINTING GRAPHS ####
        #-------os.chdir(self.graphs_path)
        self.post_utils.ComputeMeanVelocitiesInTrap("Average_Velocity.txt", self.time, self.graphs_path)
        self.materialTest.MeasureForcesAndPressure()
        self.materialTest.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintGraph(self.time)
        self.DEMFEMProcedures.PrintBallsGraph(self.time)
        self.DEMEnergyCalculator.CalculateEnergyAndPlot(self.time)
        self.BeforePrintingOperations(self.time)
        self.PrintResults()
        self.FinalizeTimeStep(self.time)

    def AfterSolveOperations(self):
        if self.post_normal_impact_velocity_option:
            self.particle_watcher.MakeMeasurements(self.analytic_model_part)
            if self.IsTimeToPrintPostProcess(self.time):
                self.particle_watcher.SetNodalMaxImpactVelocities(self.analytic_model_part)
                self.particle_watcher.SetNodalMaxFaceImpactVelocities(self.analytic_model_part)

        #Phantom Walls
        self.RunAnalytics(self.time, self.IsTimeToPrintPostProcess(self.time))

    def FinalizeTimeStep(self, time):
        pass

    def BreakSolutionStepsLoop(self):
        return False

    def Finalize(self):

        self.KRATOSprint("Finalizing execution...")
        self.GraphicalOutputFinalize()
        self.materialTest.FinalizeGraphs()
        self.DEMFEMProcedures.FinalizeGraphs(self.rigid_face_model_part)
        self.DEMFEMProcedures.FinalizeBallsGraphs(self.spheres_model_part)
        self.DEMEnergyCalculator.FinalizeEnergyPlot()
        self.CleanUpOperations()

        #------os.chdir(self.main_path)

    def CleanUpOperations(self):

        self.procedures.DeleteFiles()

        self.KRATOSprint(self.report.FinalReport(timer))

        if self.post_normal_impact_velocity_option:
            del self.analytic_model_part

        del self.KRATOSprint
        del self.all_model_parts
        del self.demio
        del self.procedures
        del self.creator_destructor
        del self.dem_fem_search
        del self.solver
        del self.DEMFEMProcedures
        del self.post_utils
        del self.cluster_model_part
        del self.rigid_face_model_part
        del self.spheres_model_part
        del self.DEM_inlet_model_part
        del self.mapping_model_part

        if self.DEM_parameters["dem_inlet_option"].GetBool():
            del self.DEM_inlet

    def SetGraphicalOutput(self):
        self.demio = DEM_procedures.DEMIo(self.model, self.DEM_parameters, self.post_path, self.all_model_parts)
        if self.DEM_parameters["post_vtk_option"].GetBool():
            import dem_vtk_output
            self.vtk_output = dem_vtk_output.VtkOutput(self.main_path, self.problem_name, self.spheres_model_part, self.rigid_face_model_part)

    def GraphicalOutputInitialize(self):
        self.demio.Initialize(self.DEM_parameters)

        #-------------os.chdir(self.post_path)
        self.demio.InitializeMesh(self.all_model_parts)

    def PrintResultsForGid(self, time):
        if self.solver.poisson_ratio_option:
            self.DEMFEMProcedures.PrintPoisson(self.spheres_model_part, self.DEM_parameters, "Poisson_ratio.txt", time)

        if self.DEM_parameters["PostEulerAngles"].GetBool():
            self.post_utils.PrintEulerAngles(self.spheres_model_part, self.cluster_model_part)

        self.demio.ShowPrintingResultsOnScreen(self.all_model_parts)

        #------os.chdir(self.data_and_results)
        self.demio.PrintMultifileLists(time, self.post_path)
        self.solver.PrepareElementsForPrinting()
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self.solver.PrepareContactElementsForPrinting()

        #os.chdir(self.post_path)

        self.demio.PrintResults(self.all_model_parts, self.creator_destructor, self.dem_fem_search, time, self.bounding_box_time_limits)
        #------os.chdir(self.main_path)
        if "post_vtk_option" in self.DEM_parameters.keys():
            if self.DEM_parameters["post_vtk_option"].GetBool():
                self.vtk_output.WriteResults(self.time)

    def GraphicalOutputFinalize(self):
        self.demio.FinalizeMesh()
        self.demio.CloseMultifiles()


    #these functions are needed for coupling, so that single time loops can be done

    def InitializeTime(self):
        self.step = 0
        self.time = 0.0
        self.time_old_print = 0.0

    def UpdateTimeParameters(self):
        self.InitializeSolutionStep()
        self.step, self.time = self._GetSolver().AdvanceInTime(self.step, self.time)
        self.DEMFEMProcedures.UpdateTimeInModelParts(self.all_model_parts, self.time, self.solver.dt, self.step)

    def FinalizeSingleTimeStep(self):
        self.DEMFEMProcedures.MoveAllMeshes(self.all_model_parts, self.time, self.solver.dt)
        #DEMFEMProcedures.MoveAllMeshesUsingATable(rigid_face_model_part, time, dt)
        ##### adding DEM elements by the inlet ######
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            self.DEM_inlet.CreateElementsFromInletMesh(self.spheres_model_part, self.cluster_model_part, self.creator_destructor)  # After solving, to make sure that neighbours are already set.
        print(self.time,self.step)
        stepinfo = self.report.StepiReport(timer, self.time, self.step)
        if stepinfo:
            self.KRATOSprint(stepinfo)

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
        if self.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.solver.dt:
            self.PrintResultsForGid(self.time)
            self.time_old_print = self.time
        self.FinalizeTimeStep(self.time)


if __name__ == "__main__":
    model = Model()
    Solution(model).Run()
