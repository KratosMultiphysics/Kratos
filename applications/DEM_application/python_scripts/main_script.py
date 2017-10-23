from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
import shutil
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
sys.path.insert(0,'')

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ or "I_MPI_INFO_NUMA_NODE_NUM" in os.environ:
    print("Running under MPI...........")
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
else:
    print("Running under OpenMP........")
    import DEM_procedures
    import DEM_material_test_script

class Solution(object):

    def LoadParametersFile(self):
        parameters_file = open("ProjectParametersDEM.json",'r')
        self.DEM_parameters = Parameters(parameters_file.read())

    def model_part_reader(self, modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)

    def __init__(self):

        self.LoadParametersFile()
        self.solver_strategy = self.SetSolverStrategy()
        self.creator_destructor = self.SetParticleCreatorDestructor()
        self.dem_fem_search = self.SetDemFemSearch()
        self.procedures = self.SetProcedures()
        self.SetAnalyticParticleWatcher()
        self.PreUtilities = PreUtilities()


        # Creating necessary directories:
        self.main_path = os.getcwd()
        problem_name = self.GetProblemTypeFilename()
        [self.post_path, self.data_and_results, self.graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(problem_name))

        self.SetGraphicalOutput()
        self.report        = DEM_procedures.Report()
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest  = DEM_procedures.MaterialTest()
        self.scheme = self.SetScheme()

        # Define control variables
        self.p_frequency = 100   # activate every 100 steps
        self.step_count = 0
        self.p_count = self.p_frequency

        # Set the print function TO_DO: do this better...
        self.KRATOSprint   = self.procedures.KRATOSprint

        # Prepare modelparts
        self.spheres_model_part    = ModelPart("SpheresPart")
        self.rigid_face_model_part = ModelPart("RigidFacePart")
        self.cluster_model_part    = ModelPart("ClusterPart")
        self.DEM_inlet_model_part  = ModelPart("DEMInletPart")
        self.mapping_model_part    = ModelPart("MappingPart")
        self.contact_model_part    = ModelPart("ContactPart")

        mp_list = []
        mp_list.append(self.spheres_model_part)
        mp_list.append(self.rigid_face_model_part)
        mp_list.append(self.cluster_model_part)
        mp_list.append(self.DEM_inlet_model_part)
        mp_list.append(self.mapping_model_part)
        mp_list.append(self.contact_model_part)

        self.all_model_parts = DEM_procedures.SetOfModelParts(mp_list)
        self.solver = self.SetSolver()
        #self.final_time = DEM_parameters.FinalTime
        #self.dt = DEM_parameters.MaxTimeStep
        self.Setdt()
        self.SetFinalTime()

    def IsCountStep(self):
        self.step_count += 1
        if self.step_count == self.p_count:
           self.p_count += self.p_frequency
           return True
        else:
            return False

    def SetAnalyticParticleWatcher(self):
        self.main_path = os.getcwd()  #revisar
        from analytic_tools import analytic_data_procedures
        self.particle_watcher = AnalyticParticleWatcher()
        self.particle_watcher_analyser = analytic_data_procedures.ParticleWatcherAnalyzer(analytic_particle_watcher = self.particle_watcher, path = self.main_path)

    def SetFinalTime(self):
        self.final_time = self.DEM_parameters["FinalTime"].GetDouble()

    def Setdt(self):
        self.dt = self.DEM_parameters["MaxTimeStep"].GetDouble()

    def SetProcedures(self):
        return DEM_procedures.Procedures(self.DEM_parameters)

    def SetDemFemSearch(self):
        return DEM_FEM_Search()

    def GetParticleHistoryWatcher(self):
        return None

    def SetParticleCreatorDestructor(self):

        self.watcher = self.GetParticleHistoryWatcher()

        if self.watcher == None:
            return ParticleCreatorDestructor()
        else:
            return ParticleCreatorDestructor(self.watcher)

    def SelectScheme(self):
        if (self.DEM_parameters["IntegrationScheme"].GetString() == 'Forward_Euler'):
            return ForwardEulerScheme()
        elif (self.DEM_parameters["IntegrationScheme"].GetString() == 'Symplectic_Euler'):
            return SymplecticEulerScheme()
        elif (self.DEM_parameters["IntegrationScheme"].GetString() == 'Taylor_Scheme'):
            return TaylorScheme()
        elif (self.DEM_parameters["IntegrationScheme"].GetString() == 'Newmark_Beta_Method'):
            return NewmarkBetaScheme(0.5, 0.25)
        elif (self.DEM_parameters["IntegrationScheme"].GetString() == 'Verlet_Velocity'):
            return VerletVelocityScheme()
        else:
            return None

    def SetScheme(self):
        scheme = self.SelectScheme()

        if scheme == None:
            self.KRATOSprint('Error: selected scheme not defined. Please select a different scheme')
            sys.exit("\nExecution was aborted.\n")
        return scheme

    def SetSolverStrategy(self):

        # TODO: Ugly fix. Change it. I don't like this to be in the main...
        # Strategy object
        if (self.DEM_parameters["ElementType"].GetString() == "SphericPartDEMElement3D" or self.DEM_parameters["ElementType"].GetString() == "CylinderPartDEMElement2D"):
            import sphere_strategy as SolverStrategy
        elif (self.DEM_parameters["ElementType"].GetString() == "SphericContPartDEMElement3D" or self.DEM_parameters["ElementType"].GetString() == "CylinderContPartDEMElement2D"):
            import continuum_sphere_strategy as SolverStrategy
        elif (self.DEM_parameters["ElementType"].GetString() == "ThermalSphericContPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (self.DEM_parameters["ElementType"].GetString() == "ThermalSphericPartDEMElement3D"):
            import thermal_sphere_strategy as SolverStrategy
        elif (self.DEM_parameters["ElementType"].GetString() == "SinteringSphericConPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (self.DEM_parameters["ElementType"].GetString() == "IceContPartDEMElement3D"):
            import ice_continuum_sphere_strategy as SolverStrategy
        else:
            self.KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')

        return SolverStrategy


    def SetSolver(self):
        return self.solver_strategy.ExplicitStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, self.DEM_parameters, self.procedures)

    def Run(self):
        self.Initialize()

        self.RunMainTemporalLoop()

        self.Finalize()

        self.CleanUpOperations()

    def AddVariables(self):
        self.procedures.AddAllVariablesInAllModelParts(self.solver, self.scheme, self.all_model_parts, self.DEM_parameters)

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
        self.AddVariables()

        self.ReadModelParts()

        self.FillAnalyticSubModelParts()

        # Setting up the buffer size
        self.procedures.SetUpBufferSizeInAllModelParts(self.spheres_model_part, 1, self.cluster_model_part, 1, self.DEM_inlet_model_part, 1, self.rigid_face_model_part, 1)

        # Adding dofs
        self.solver.AddDofs(self.spheres_model_part)
        self.solver.AddDofs(self.cluster_model_part)
        self.solver.AddDofs(self.DEM_inlet_model_part)

        os.chdir(self.main_path)
        self.KRATOSprint("\nInitializing Problem...")

        self.GraphicalOutputInitialize()

        # Perform a partition to balance the problem
        self.solver.search_strategy = self.parallelutils.GetSearchStrategy(self.solver, self.spheres_model_part)
        self.solver.BeforeInitialize()
        self.parallelutils.Repart(self.spheres_model_part)

        #Setting up the BoundingBox
        self.bounding_box_time_limits = self.procedures.SetBoundingBoxLimits(self.all_model_parts, self.creator_destructor)

        #Finding the max id of the nodes... (it is necessary for anything that will add spheres to the self.spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.z
        max_Id = self.procedures.FindMaxNodeIdAccrossModelParts(self.creator_destructor, self.all_model_parts)

        self.creator_destructor.SetMaxNodeId(self.all_model_parts.MaxNodeId)

        #Strategy Initialization
        os.chdir(self.main_path)
        self.solver.Initialize() # Possible modifications of number of elements and number of nodes

        #Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation
        #Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
        self.SetInlet()

        self.SetInitialNodalValues()

        self.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(self.DEM_parameters, self.graphs_path, self.spheres_model_part, self.rigid_face_model_part)

        os.chdir(self.graphs_path)
        self.DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(self.DEM_parameters, self.spheres_model_part, self.cluster_model_part, "EnergyPlot.grf")

        self.materialTest.Initialize(self.DEM_parameters, self.procedures, self.solver, self.graphs_path, self.post_path, self.spheres_model_part, self.rigid_face_model_part)

        self.KRATOSprint("Initialization Complete" + "\n")

        self.report.Prepare(timer, self.DEM_parameters["ControlTime"].GetDouble())

        #self.procedures.ModelData(self.spheres_model_part, self.solver) #check link with ModelDataInfo = "OFF"

        self.materialTest.PrintChart()
        self.materialTest.PrepareDataForGraph()

        self.post_utils = DEM_procedures.PostUtils(self.DEM_parameters, self.spheres_model_part)
        self.report.total_steps_expected = int(self.final_time / self.dt)
        self.KRATOSprint(self.report.BeginReport(timer))

    def GetMpFilename(self):
        return self.DEM_parameters["problem_name"].GetString() + "DEM"

    def GetInletFilename(self):
        return self.DEM_parameters["problem_name"].GetString() + "DEM_Inlet"

    def GetFemFilename(self):
        return self.DEM_parameters["problem_name"].GetString() + "DEM_FEM_boundary"

    def GetClusterFilename(self):
        return self.DEM_parameters["problem_name"].GetString() + "DEM_Clusters"

    def GetProblemTypeFilename(self):
        return self.DEM_parameters["problem_name"].GetString()

    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):
        os.chdir(self.main_path)

        # Reading the model_part
        spheres_mp_filename   = self.GetMpFilename()
        model_part_io_spheres = self.model_part_reader(spheres_mp_filename, max_node_Id, max_elem_Id, max_cond_Id)

        if "do_not_perform_initial_partition" in self.DEM_parameters.keys() and self.DEM_parameters["do_not_perform_initial_partition"].GetBool():
            pass
        else:
            self.parallelutils.PerformInitialPartition(model_part_io_spheres)

        os.chdir(self.main_path)
        [model_part_io_spheres, self.spheres_model_part, MPICommSetup] = self.parallelutils.SetCommunicator(self.spheres_model_part, model_part_io_spheres, spheres_mp_filename)
        model_part_io_spheres.ReadModelPart(self.spheres_model_part)

        max_node_Id += self.creator_destructor.FindMaxNodeIdInModelPart(self.spheres_model_part)
        max_elem_Id += self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        old_max_elem_Id_spheres = max_elem_Id
        max_cond_Id += self.creator_destructor.FindMaxConditionIdInModelPart(self.spheres_model_part)
        rigidFace_mp_filename   = self.GetFemFilename()
        model_part_io_fem = self.model_part_reader(rigidFace_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_fem.ReadModelPart(self.rigid_face_model_part)

        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.rigid_face_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.rigid_face_model_part)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.rigid_face_model_part)
        clusters_mp_filename   = self.GetClusterFilename()
        model_part_io_clusters = self.model_part_reader(clusters_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_clusters.ReadModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        if (max_elem_Id != old_max_elem_Id_spheres):
            self.creator_destructor.RenumberElementIdsFromGivenValue(self.cluster_model_part, max_elem_Id)

        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.cluster_model_part)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.cluster_model_part)
        DEM_Inlet_filename = self.GetInletFilename()
        model_part_io_demInlet = self.model_part_reader(DEM_Inlet_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_demInlet.ReadModelPart(self.DEM_inlet_model_part)

        self.model_parts_have_been_read = True

        self.all_model_parts.ComputeMaxIds()


    def RunMainTemporalLoop(self):

        self.step           = 0
        self.time           = 0.0
        self.time_old_print = 0.0
        while (self.time < self.final_time):

            self.InitializeTimeStep()
            self.time  = self.time + self.dt
            self.step += 1

            self.DEMFEMProcedures.UpdateTimeInModelParts(self.all_model_parts, self.time,self.dt,self.step)

            self.BeforeSolveOperations(self.time)

            #### SOLVE #########################################
            self.solver.Solve()
            ####################################################

            self.AfterSolveOperations()

            self.DEMFEMProcedures.MoveAllMeshes(self.all_model_parts, self.time, self.dt)
            #DEMFEMProcedures.MoveAllMeshesUsingATable(rigid_face_model_part, time, dt)

            ##### adding DEM elements by the inlet ######
            if self.DEM_parameters["dem_inlet_option"].GetBool():
                self.DEM_inlet.CreateElementsFromInletMesh(self.spheres_model_part, self.cluster_model_part, self.creator_destructor)  # After solving, to make sure that neighbours are already set.

            stepinfo = self.report.StepiReport(timer,self.time,self.step)
            if stepinfo:
                self.KRATOSprint(stepinfo)

            #### PRINTING GRAPHS ####
            os.chdir(self.graphs_path)
            self.post_utils.ComputeMeanVelocitiesinTrap("Average_Velocity.txt", self.time)

            self.materialTest.MeasureForcesAndPressure()
            self.materialTest.PrintGraph(self.time)

            self.DEMFEMProcedures.PrintGraph(self.time)
            self.DEMFEMProcedures.PrintBallsGraph(self.time)

            self.DEMEnergyCalculator.CalculateEnergyAndPlot(self.time)

            self.BeforePrintingOperations(self.time)

            #### GiD IO ##########################################
            time_to_print = self.time - self.time_old_print

            if (self.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.dt):

                self.PrintResultsForGid(self.time)
                self.time_old_print = self.time

            self.FinalizeTimeStep(self.time)


    def SetInlet(self):
        if self.DEM_parameters["dem_inlet_option"].GetBool():
            #Constructing the inlet and initializing it (must be done AFTER the self.spheres_model_part Initialize)
            self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor, self.solver.continuum_type)

    def SetInitialNodalValues(self):
        self.procedures.SetInitialNodalValues(self.spheres_model_part, self.cluster_model_part, self.DEM_inlet_model_part, self.rigid_face_model_part)

    def InitializeTimeStep(self):
        pass

    def BeforeSolveOperations(self, time):
       if "PostNormalImpactVelocity" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostNormalImpactVelocity"].GetBool():
                if self.IsCountStep():
                    #time_to_print = self.time - self.time_old_print    # add new particles to analytic mp each time an output is generated
                    #if (self.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.dt):
                        self.FillAnalyticSubModelPartsWithNewParticles()

    def BeforePrintingOperations(self, time):
        pass

    def AfterSolveOperations(self):
        if "PostNormalImpactVelocity" in self.DEM_parameters.keys():
            if self.DEM_parameters["PostNormalImpactVelocity"].GetBool():
                self.particle_watcher.MakeMeasurements(self.analytic_model_part)
                time_to_print = self.time - self.time_old_print
                if (self.DEM_parameters["OutputTimeStep"].GetDouble() - time_to_print < 1e-2 * self.dt):
                    self.particle_watcher.SetNodalMaxImpactVelocities(self.analytic_model_part)
                    self.particle_watcher.SetNodalMaxFaceImpactVelocities(self.analytic_model_part)

    def FinalizeTimeStep(self, time):
        pass

    def Finalize(self):

        self.KRATOSprint("Finalizing execution...")

        self.GraphicalOutputFinalize()

        self.materialTest.FinalizeGraphs()
        self.DEMFEMProcedures.FinalizeGraphs(self.rigid_face_model_part)
        self.DEMFEMProcedures.FinalizeBallsGraphs(self.spheres_model_part)
        self.DEMEnergyCalculator.FinalizeEnergyPlot()


        os.chdir(self.main_path)

    def CleanUpOperations(self):

        objects_to_destroy = [self.demio, self.procedures, self.creator_destructor, self.dem_fem_search, self.solver, self.DEMFEMProcedures, self.post_utils,
                              self.cluster_model_part, self.rigid_face_model_part, self.spheres_model_part, self.DEM_inlet_model_part, self.mapping_model_part]

        if self.DEM_parameters["dem_inlet_option"].GetBool():
            objects_to_destroy.append(self.DEM_inlet)

        for obj in objects_to_destroy:
            del obj

        self.procedures.DeleteFiles()

        self.KRATOSprint(self.report.FinalReport(timer))

    def SetGraphicalOutput(self):
        self.demio         = DEM_procedures.DEMIo(self.DEM_parameters, self.post_path)

    def GraphicalOutputInitialize(self):
        self.demio.Initialize(self.DEM_parameters)

        os.chdir(self.post_path)
        self.demio.InitializeMesh(self.all_model_parts)

    def PrintResultsForGid(self, time):
        if self.solver.poisson_ratio_option:
            self.DEMFEMProcedures.PrintPoisson(self.spheres_model_part, self.DEM_parameters, "Poisson_ratio.txt", time)

        if self.DEM_parameters["PostEulerAngles"].GetBool():
            self.post_utils.PrintEulerAngles(self.spheres_model_part, self.cluster_model_part)

        self.demio.ShowPrintingResultsOnScreen(self.all_model_parts)

        os.chdir(self.data_and_results)
        self.demio.PrintMultifileLists(time, self.post_path)
        os.chdir(self.post_path)

        self.solver.PrepareElementsForPrinting()
        if self.DEM_parameters["ContactMeshOption"].GetBool():
            self.solver.PrepareContactElementsForPrinting()

        self.demio.PrintResults(self.all_model_parts, self.creator_destructor, self.dem_fem_search, time, self.bounding_box_time_limits)
        os.chdir(self.main_path)

    def GraphicalOutputFinalize(self):
        self.demio.FinalizeMesh()
        self.demio.CloseMultifiles()

if __name__ == "__main__":
    Solution().Run()
