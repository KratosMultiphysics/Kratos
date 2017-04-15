from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import time as timer
import os
import sys
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

sys.path.insert(0,'')
import DEM_explicit_solver_var as DEM_parameters

# Import MPI modules if needed. This way to do this is only valid when using OpenMPI. For other implementations of MPI it will not work.
if "OMPI_COMM_WORLD_SIZE" in os.environ or "I_MPI_INFO_NUMA_NODE_NUM" in os.environ:
    print("Running under MPI...........")
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.MPISearchApplication import *
    from KratosMultiphysics.mpi import *
    import DEM_procedures_mpi as DEM_procedures
    import DEM_material_test_script_mpi as DEM_material_test_script
    def model_part_reader(modelpart, nodeid=0, elemid=0, condid=0):
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)

else:
    print("Running under OpenMP........")
    import DEM_procedures
    import DEM_material_test_script
    #def model_part_reader(modelpart, a=0, b=0, c=0):
    #    return ModelPartIO(modelpart)
    def model_part_reader(modelpart, nodeid=0, elemid=0, condid=0):
        #return ModelPartIO(modelpart)
        return ReorderConsecutiveFromGivenIdsModelPartIO(modelpart, nodeid, elemid, condid)


class Solution(object):

    def __init__(self):
        self.solver_strategy = self.SetSolverStrategy()
        self.creator_destructor = self.SetParticleCreatorDestructor()
        self.dem_fem_search = self.SetDemFemSearch()
        self.procedures = self.SetProcedures()

        self.procedures.CheckInputParameters(DEM_parameters)

        # Creating necessary directories:
        self.main_path = os.getcwd()
        [self.post_path, self.data_and_results, self.graphs_path, MPI_results] = self.procedures.CreateDirectories(str(self.main_path), str(DEM_parameters.problem_name))

        self.demio         = DEM_procedures.DEMIo(DEM_parameters, self.post_path)
        self.report        = DEM_procedures.Report()
        self.parallelutils = DEM_procedures.ParallelUtils()
        self.materialTest  = DEM_procedures.MaterialTest()
        self.scheme = self.SetScheme()

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


    def SetProcedures(self):
        return DEM_procedures.Procedures(DEM_parameters)

    def SetDemFemSearch(self):
        return DEM_FEM_Search()

    def SetParticleCreatorDestructor(self):
        return ParticleCreatorDestructor()

    def SelectScheme(self):
        if (DEM_parameters.IntegrationScheme == 'Forward_Euler'):
            return ForwardEulerScheme()
        elif (DEM_parameters.IntegrationScheme == 'Symplectic_Euler'):
            return SymplecticEulerScheme()
        elif (DEM_parameters.IntegrationScheme == 'Taylor_Scheme'):
            return TaylorScheme()
        elif (DEM_parameters.IntegrationScheme == 'Newmark_Beta_Method'):
            return NewmarkBetaScheme(0.5, 0.25)
        elif (DEM_parameters.IntegrationScheme == 'Verlet_Velocity'):
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
        if (DEM_parameters.ElementType == "SphericPartDEMElement3D" or DEM_parameters.ElementType == "CylinderPartDEMElement2D"):
            import sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "SphericContPartDEMElement3D" or DEM_parameters.ElementType == "CylinderContPartDEMElement2D"):
            import continuum_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "ThermalSphericContPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "ThermalSphericPartDEMElement3D"):
            import thermal_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "SinteringSphericConPartDEMElement3D"):
            import thermal_continuum_sphere_strategy as SolverStrategy
        elif (DEM_parameters.ElementType == "IceContPartDEMElement3D"):
            import ice_continuum_sphere_strategy as SolverStrategy
        else:
            self.KRATOSprint('Error: Strategy unavailable. Select a different scheme-element')

        return SolverStrategy


    def SetSolver(self):
        return self.solver_strategy.ExplicitStrategy(self.all_model_parts, self.creator_destructor, self.dem_fem_search, self.scheme, DEM_parameters, self.procedures)


    def Run(self):

        self.Initialize()

        self.RunMainTemporalLoop()

        self.Finalize()

    def AddVariables(self):
        self.procedures.AddAllVariablesInAllModelParts(self.solver, self.scheme, self.all_model_parts, DEM_parameters)

    def Initialize(self):
        self.AddVariables()

        self.ReadModelParts()

        # Setting up the buffer size
        self.procedures.SetUpBufferSizeInAllModelParts(self.spheres_model_part, 1, self.cluster_model_part, 1, self.DEM_inlet_model_part, 1, self.rigid_face_model_part, 1)
        # Adding dofs

        self.solver.AddDofs(self.spheres_model_part)
        self.solver.AddDofs(self.cluster_model_part)

        self.solver.AddDofs(self.DEM_inlet_model_part)

        os.chdir(self.main_path)

        self.KRATOSprint("\nInitializing Problem...")

        self.demio.Initialize(DEM_parameters)

        os.chdir(self.post_path)
        self.demio.InitializeMesh(self.all_model_parts)

        # Perform a partition to balance the problem
        self.solver.search_strategy = self.parallelutils.GetSearchStrategy(self.solver, self.spheres_model_part)
        self.solver.BeforeInitialize()
        self.parallelutils.Repart(self.spheres_model_part)

        #Setting up the BoundingBox
        self.bounding_box_time_limits = self.procedures.SetBoundingBoxLimits(self.all_model_parts, self.creator_destructor)

        self.dt = DEM_parameters.MaxTimeStep

        #Finding the max id of the nodes... (it is necessary for anything that will add spheres to the self.spheres_model_part, for instance, the INLETS and the CLUSTERS read from mdpa file.z
        max_Id = self.procedures.FindMaxNodeIdAccrossModelParts(self.creator_destructor, self.all_model_parts)
        
        self.creator_destructor.SetMaxNodeId(self.all_model_parts.MaxNodeId)

        #Strategy Initialization
        os.chdir(self.main_path)
        self.solver.Initialize() # Possible modifications of number of elements and number of nodes
        self.dt = min(DEM_parameters.MaxTimeStep, self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME)) # under revision. linked to automatic timestep? Possible modifications of DELTA_TIME
        #Constructing a model part for the DEM inlet. It contains the DEM elements to be released during the simulation
        #Initializing the DEM solver must be done before creating the DEM Inlet, because the Inlet configures itself according to some options of the DEM model part
        if DEM_parameters.dem_inlet_option:
            #Constructing the inlet and initializing it (must be done AFTER the self.spheres_model_part Initialize)
            self.DEM_inlet = DEM_Inlet(self.DEM_inlet_model_part)
            self.DEM_inlet.InitializeDEM_Inlet(self.spheres_model_part, self.creator_destructor, self.solver.continuum_type)


        self.DEMFEMProcedures = DEM_procedures.DEMFEMProcedures(DEM_parameters, self.graphs_path, self.spheres_model_part, self.rigid_face_model_part)

        os.chdir(self.graphs_path)
        self.DEMEnergyCalculator = DEM_procedures.DEMEnergyCalculator(DEM_parameters, self.spheres_model_part, self.cluster_model_part, "EnergyPlot.grf")

        self.materialTest.Initialize(DEM_parameters, self.procedures, self.solver, self.graphs_path, self.post_path, self.spheres_model_part, self.rigid_face_model_part)

        self.KRATOSprint("Initialization Complete" + "\n")

        self.report.Prepare(timer, DEM_parameters.ControlTime)

        self.procedures.ModelData(self.spheres_model_part, self.solver)

        self.materialTest.PrintChart()
        self.materialTest.PrepareDataForGraph()

        self.post_utils = DEM_procedures.PostUtils(DEM_parameters, self.spheres_model_part)

        self.report.total_steps_expected = int(DEM_parameters.FinalTime / self.dt)
        self.KRATOSprint(self.report.BeginReport(timer))


    def ReadModelParts(self, max_node_Id = 0, max_elem_Id = 0, max_cond_Id = 0):

        os.chdir(self.main_path)
        # Reading the model_part
        spheres_mp_filename   = DEM_parameters.problem_name + "DEM"
        model_part_io_spheres = model_part_reader(spheres_mp_filename, max_node_Id, max_elem_Id, max_cond_Id)

        if (hasattr(DEM_parameters, "do_not_perform_initial_partition") and DEM_parameters.do_not_perform_initial_partition == 1):
            pass
        else:
            self.parallelutils.PerformInitialPartition(model_part_io_spheres)

        os.chdir(self.main_path)
        [model_part_io_spheres, self.spheres_model_part, MPICommSetup] = self.parallelutils.SetCommunicator(self.spheres_model_part, model_part_io_spheres, spheres_mp_filename)

        model_part_io_spheres.ReadModelPart(self.spheres_model_part)

        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.spheres_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        old_max_elem_Id_spheres = max_elem_Id
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.spheres_model_part)
        rigidFace_mp_filename = DEM_parameters.problem_name + "DEM_FEM_boundary"
        model_part_io_fem = model_part_reader(rigidFace_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_fem.ReadModelPart(self.rigid_face_model_part)

        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.rigid_face_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.rigid_face_model_part)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.rigid_face_model_part)
        clusters_mp_filename = DEM_parameters.problem_name + "DEM_Clusters"
        model_part_io_clusters = model_part_reader(clusters_mp_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_clusters.ReadModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.spheres_model_part)
        if (max_elem_Id != old_max_elem_Id_spheres):
            self.creator_destructor.RenumberElementIdsFromGivenValue(self.cluster_model_part, max_elem_Id)

        max_node_Id = self.creator_destructor.FindMaxNodeIdInModelPart(self.cluster_model_part)
        max_elem_Id = self.creator_destructor.FindMaxElementIdInModelPart(self.cluster_model_part)
        max_cond_Id = self.creator_destructor.FindMaxConditionIdInModelPart(self.cluster_model_part)
        DEM_Inlet_filename = DEM_parameters.problem_name + "DEM_Inlet"
        
        
        model_part_io_demInlet = model_part_reader(DEM_Inlet_filename,max_node_Id+1, max_elem_Id+1, max_cond_Id+1)
        model_part_io_demInlet.ReadModelPart(self.DEM_inlet_model_part)

        self.model_parts_have_been_read = True
        
        self.all_model_parts.ComputeMaxIds()


    def RunMainTemporalLoop(self):

        self.step           = 0
        self.time           = 0.0
        time_old_print = 0.0

        while (self.time < DEM_parameters.FinalTime):

            self.InitializeTimeStep()

            self.dt    = self.spheres_model_part.ProcessInfo.GetValue(DELTA_TIME) # Possible modifications of DELTA_TIME
            self.time  = self.time + self.dt
            self.step += 1

            self.DEMFEMProcedures.UpdateTimeInModelParts(self.all_model_parts, self.time,self.dt,self.step)

            self.BeforeSolveOperations()

            #### SOLVE #########################################
            self.solver.Solve()
            ####################################################

            self.AfterSolveOperations()

            self.DEMFEMProcedures.MoveAllMeshes(self.all_model_parts, self.time, self.dt)
            #DEMFEMProcedures.MoveAllMeshesUsingATable(rigid_face_model_part, time, dt)

            ##### adding DEM elements by the inlet ######
            if (DEM_parameters.dem_inlet_option):
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

            #### GiD IO ##########################################
            time_to_print = self.time - time_old_print

            if (DEM_parameters.OutputTimeStep - time_to_print < 1e-2 * self.dt):

                self.PrintResultsForGid(self.time)
                time_old_print = self.time

            self.FinalizeTimeStep()


    def PrintResultsForGid(self, time):

        if self.solver.poisson_ratio_option:
            self.DEMFEMProcedures.PrintPoisson(self.spheres_model_part, DEM_parameters, "Poisson_ratio.txt", time)

        if DEM_parameters.PostEulerAngles:
            self.post_utils.PrintEulerAngles(self.spheres_model_part, self.cluster_model_part)

        self.demio.ShowPrintingResultsOnScreen(self.all_model_parts)

        os.chdir(self.data_and_results)
        self.demio.PrintMultifileLists(time, self.post_path)
        os.chdir(self.post_path)

        self.solver.PrepareElementsForPrinting()
        if (DEM_parameters.ContactMeshOption == "ON"):
            self.solver.PrepareContactElementsForPrinting()

        self.demio.PrintResults(self.all_model_parts, self.creator_destructor, self.dem_fem_search, time, self.bounding_box_time_limits)
        os.chdir(self.main_path)


    def InitializeTimeStep(self):
        pass


    def BeforeSolveOperations(self):
        pass


    def AfterSolveOperations(self):
        pass


    def FinalizeTimeStep(self):
        pass


    def Finalize(self):

        self.KRATOSprint("Finalizing execution...")

        self.demio.FinalizeMesh()
        self.materialTest.FinalizeGraphs()
        self.DEMFEMProcedures.FinalizeGraphs(self.rigid_face_model_part)
        self.DEMFEMProcedures.FinalizeBallsGraphs(self.spheres_model_part)
        self.DEMEnergyCalculator.FinalizeEnergyPlot()
        self.demio.CloseMultifiles()

        os.chdir(self.main_path)

        objects_to_destroy = [self.demio, self.procedures, self.creator_destructor, self.dem_fem_search, self.solver, self.DEMFEMProcedures, self.post_utils,
                              self.cluster_model_part, self.rigid_face_model_part, self.spheres_model_part, self.DEM_inlet_model_part, self.mapping_model_part]

        if (DEM_parameters.dem_inlet_option):
            objects_to_destroy.append(self.DEM_inlet)

        for obj in objects_to_destroy:
            del obj

        self.procedures.DeleteFiles()

        self.KRATOSprint(self.report.FinalReport(timer))
