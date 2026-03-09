

import KratosMultiphysics
import KratosMultiphysics.FemToDemApplication as KratosFemDem
import KratosMultiphysics.StructuralMechanicsApplication as KratosSMA

import KratosMultiphysics.FemToDemApplication.MainDEM_for_coupling as DEM
import KratosMultiphysics.FemToDemApplication.MainFEM_for_coupling as FEM
import KratosMultiphysics.FemToDemApplication.FEMDEMParticleCreatorDestructor as PCD
import math
import os
import KratosMultiphysics.DEMApplication as KratosDEM
import KratosMultiphysics.DemStructuresCouplingApplication as DemFem
import KratosMultiphysics.FemToDemApplication.fem_dem_coupled_gid_output as gid_output

def Wait():
    input("Press Something")

#============================================================================================================================
class MainCoupledFemDem_Solution:
#============================================================================================================================
    def __init__(self, Model, path = ""):
        self.model = Model
        # Initialize solutions

        if path == "":
            DEMProjectParametersFile = open("ProjectParametersDEM.json", 'r')
        else:
            DEMProjectParametersFile = open(os.path.join(path, "ProjectParametersDEM.json"), 'r')
        DEM_project_parameters = KratosMultiphysics.Parameters(DEMProjectParametersFile.read())

        self.FEM_Solution = FEM.FEM_for_coupling_Solution(Model, path)
        self.DEM_Solution = DEM.DEM_for_coupling_Solution(Model, DEM_project_parameters)

        self.domain_size = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # self.InitializePlotsFiles()
        self.echo_level = 0
        self.is_slave = False

#============================================================================================================================
    def Run(self):
        self.Initialize()
        self.RunMainTemporalLoop()
        self.Finalize()

#============================================================================================================================
    def SetNonHistoricalVariables(self):
        nodes = self.FEM_Solution.main_model_part.Nodes
        utils = KratosMultiphysics.VariableUtils()
        # Initialize the "flag" IS_DEM in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.IS_DEM, False, nodes)
        # Initialize the "flag" NODAL_FORCE_APPLIED in all the nodes
        utils.SetNonHistoricalVariable(KratosFemDem.NODAL_FORCE_APPLIED, False, nodes)
        # Initialize the "flag" RADIUS in all the nodes
        utils.SetNonHistoricalVariable(KratosMultiphysics.RADIUS, 0.0, nodes)

        # Initialize the var to track volume erased for each pressure
        utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_VOLUME, 0.0, nodes)
        utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_INITIAL_VOLUME, 0.0, nodes)

#============================================================================================================================

    def InitializeProcessesAndVariables(self):
        self.SetNonHistoricalVariables()

        self.SpheresModelPart = self.DEM_Solution.spheres_model_part
        self.DEMParameters = self.DEM_Solution.DEM_parameters
        self.DEMProperties = self.SpheresModelPart.GetProperties()[1]
        self.ParticleCreatorDestructor = PCD.FemDemParticleCreatorDestructor(self.SpheresModelPart,
                                                                           self.DEMProperties,
                                                                           self.DEMParameters)
        if self.domain_size == 3:
            KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.FEM_Solution.main_model_part).Execute()
            KratosMultiphysics.GenericFindElementalNeighboursProcess(self.FEM_Solution.main_model_part).Execute()

        femdem_custom_settings = self.FEM_Solution.ProjectParameters["fem_dem_settings"]
        femdem_default_settings = KratosMultiphysics.Parameters("""
            {
                "transfer_dem_contact_forces" : true,
                "pressure_load_extrapolation" : true,
                "DEM_FEM_contact"             : true,
                "tangent_operator"            : 1,
                "create_initial_skin"         : false,
                "do_stabilization_solve"      : false,
                "smoothing_of_stresses"       : true,
                "maximum_damage_erase"        : 0.98
            }""")
        femdem_custom_settings.ValidateAndAssignDefaults(femdem_default_settings)
        process_info = self.FEM_Solution.main_model_part.ProcessInfo
        self.TransferDEMContactForcesToFEM = femdem_custom_settings["transfer_dem_contact_forces"].GetBool()
        self.PressureLoad = femdem_custom_settings["pressure_load_extrapolation"].GetBool()
        self.DEMFEM_contact = femdem_custom_settings["DEM_FEM_contact"].GetBool()
        process_info[KratosFemDem.TANGENT_CONSTITUTIVE_TENSOR] = femdem_custom_settings["tangent_operator"].GetInt()
        self.CreateInitialSkin = femdem_custom_settings["create_initial_skin"].GetBool()
        self.do_stabilization_solve = femdem_custom_settings["do_stabilization_solve"].GetBool()
        process_info[KratosFemDem.SMOOTHING_OF_STRESSES] = femdem_custom_settings["smoothing_of_stresses"].GetBool()
        process_info[KratosFemDem.MAX_DAMAGE_ERASE] = femdem_custom_settings["maximum_damage_erase"].GetDouble()

        # Initialize IP variables to zero
        self.InitializeIntegrationPointsVariables()

        if self.PressureLoad:
            KratosFemDem.AssignPressureIdProcess(self.FEM_Solution.main_model_part).Execute()
            KratosFemDem.ComputeInitialVolumeProcess(self.FEM_Solution.main_model_part).Execute()

        self.SkinDetectionProcessParameters = KratosMultiphysics.Parameters("""
        {
            "name_auxiliar_model_part" : "SkinDEMModelPart",
            "name_auxiliar_condition"  : "Condition",
            "echo_level"               : 0
        }""")


        # for the dem contact forces coupling
        self.InitializeDummyNodalForces()

        # Just to find neighbours the 1st time
        self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.GENERATE_DEM] = True
        if self.domain_size == 3:
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECOMPUTE_NEIGHBOURS] = True

        if self.domain_size == 3: # only in 3D
            # We assign the flag to recompute neighbours inside the 3D elements the 1st time
            utils = KratosMultiphysics.VariableUtils()
            utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

        if self.CreateInitialSkin:
            self.ComputeSkinSubModelPart()
            if self.DEMFEM_contact:
                self.TransferFEMSkinToDEM()
            KratosFemDem.GenerateInitialSkinDEMProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart).Execute()

        # Initialize the coupled post process
        if not self.is_slave:
            self.InitializePostProcess()

        self.FindNeighboursIfNecessary()
#============================================================================================================================
    def Initialize(self):

        self.FEM_Solution.Initialize()
        self.DEM_Solution.Initialize()
        self.InitializeProcessesAndVariables()

        self.FEM_Solution.KratosPrintInfo("")
        self.FEM_Solution.KratosPrintInfo("    ______                 ___    ____                 ")
        self.FEM_Solution.KratosPrintInfo("   / ____/___   ____ ___  |__ \  / __ \ ___   ____ ___ ")
        self.FEM_Solution.KratosPrintInfo("  / /_   / _ \ / __ `__ \ __/ / / / / // _ \ / __ `__ \ ")
        self.FEM_Solution.KratosPrintInfo(" / __/  /  __// / / / / // __/ / /_/ //  __// / / / / /")
        self.FEM_Solution.KratosPrintInfo("/_/     \___//_/ /_/ /_//____//_____/ \___//_/ /_/ /_/ Application")
        self.FEM_Solution.KratosPrintInfo("                           Developed by Alejandro Cornejo")
        self.FEM_Solution.KratosPrintInfo("")


#============================================================================================================================
    def RunMainTemporalLoop(self):
        # Solving the problem (time integration)
        self.DEM_Solution.step           = 0
        self.DEM_Solution.time           = 0.0
        self.DEM_Solution.time_old_print = 0.0

        # Temporal loop
        while self.FEM_Solution.time <= self.FEM_Solution.end_time:
            self.InitializeSolutionStep()
            self.SolveSolutionStep()
            self.FinalizeSolutionStep()

#============================================================================================================================
    def InitializeSolutionStep(self):
        # Modified for the remeshing
        self.FEM_Solution.delta_time = self.ComputeDeltaTime()
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME] = self.FEM_Solution.delta_time
        self.FEM_Solution.time = self.FEM_Solution.time + self.FEM_Solution.delta_time
        self.FEM_Solution.main_model_part.CloneTimeStep(self.FEM_Solution.time)
        self.FEM_Solution.step = self.FEM_Solution.step + 1
        self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = self.FEM_Solution.step

        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeSolutionStep of the FEM part")

        self.FEM_Solution.InitializeSolutionStep()
        self.DEM_Solution._GetSolver().AdvanceInTime(self.FEM_Solution.time)
        self.DEM_Solution._GetSolver().Predict()

#============================================================================================================================
    def SolveSolutionStep(self):  # Method to perform the coupling FEM <-> DEM

        self.FEM_Solution.clock_time = self.FEM_Solution.StartTimeMeasuring()

        #### SOLVE FEM #########################################
        self.FEM_Solution.solver.Solve()
        ########################################################

        self.ExecuteBeforeGeneratingDEM()
        self.GenerateDEM() # we create the new DEM of this time step
        self.ExecuteAfterGeneratingDEM()
        self.BeforeSolveDEMOperations()

        #### SOLVE DEM #########################################
        self.DEM_Solution.SolverSolve()
        ########################################################


#============================================================================================================================
    def FinalizeSolutionStep(self):

        self.DEM_Solution.FinalizeSolutionStep()

        # to print DEM with the FEM coordinates
        self.UpdateDEMVariables()

        # Transfer the contact forces of the DEM to the FEM nodes
        if self.TransferDEMContactForcesToFEM:
            self.TransferNodalForcesToFEM()

        self.FEM_Solution.StopTimeMeasuring(self.FEM_Solution.clock_time,"Solving", False)

        # Print required info
        # self.PrintPlotsFiles()

        # MODIFIED FOR THE REMESHING
        self.FEM_Solution.GraphicalOutputExecuteFinalizeSolutionStep()

        # processes to be executed at the end of the solution step
        self.FEM_Solution.model_processes.ExecuteFinalizeSolutionStep()

        # processes to be executed before witting the output
        self.FEM_Solution.model_processes.ExecuteBeforeOutputStep()

        # processes to be executed after writing the output
        self.FEM_Solution.model_processes.ExecuteAfterOutputStep()

        if not self.is_slave:
            self.PrintResults()

#============================================================================================================================
    def Finalize(self):
        self.FEM_Solution.Finalize()
        self.DEM_Solution.Finalize()

#InitializeIntegrationPointsVariables============================================================================================================================
    def InitializeIntegrationPointsVariables(self):
        utils = KratosMultiphysics.VariableUtils()
        elements = self.FEM_Solution.main_model_part.Elements
        nodes = self.FEM_Solution.main_model_part.Nodes

        utils.SetNonHistoricalVariable(KratosFemDem.GENERATE_DEM, False, elements)
        utils.SetNonHistoricalVariable(KratosFemDem.STRESS_THRESHOLD, 0.0, elements)
        utils.SetNonHistoricalVariable(KratosFemDem.DAMAGE_ELEMENT, 0.0, elements)
        utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_EXPANDED, 0, elements)
        utils.SetNonHistoricalVariable(KratosFemDem.IS_SKIN, 0, elements)
        utils.SetNonHistoricalVariable(KratosFemDem.SMOOTHING, 0, elements)
        utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, elements)

        if self.domain_size == 3:
            utils.SetNonHistoricalVariable(KratosFemDem.VOLUME_COUNTED, False, elements)
            utils.SetNonHistoricalVariable(KratosFemDem.FEMDEM_STRESS_VECTOR, [0.0,0.0,0.0,0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.FEMDEM_STRAIN_VECTOR, [0.0,0.0,0.0,0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR_INTEGRATED, [0.0,0.0,0.0,0.0,0.0,0.0], elements)
        else: # 2D
            utils.SetNonHistoricalVariable(KratosFemDem.FEMDEM_STRESS_VECTOR, [0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.FEMDEM_STRAIN_VECTOR, [0.0,0.0,0.0], elements)
            utils.SetNonHistoricalVariable(KratosFemDem.STRESS_VECTOR_INTEGRATED, [0.0, 0.0, 0.0], elements)

        if self.PressureLoad:
            utils.SetNonHistoricalVariable(KratosFemDem.PRESSURE_ID, 0, nodes)

#InitializeDummyNodalForces============================================================================================================================
    def InitializeDummyNodalForces(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: InitializeDummyNodalForces")

        # we fill the submodel part with the nodes and dummy conditions
        max_id = self.GetMaximumConditionId()
        props = self.FEM_Solution.main_model_part.Properties[0]
        self.FEM_Solution.main_model_part.CreateSubModelPart("ContactForcesDEMConditions")
        for node in self.FEM_Solution.main_model_part.Nodes:
            self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").AddNode(node, 0)
            max_id += 1
            cond = self.FEM_Solution.main_model_part.GetSubModelPart("ContactForcesDEMConditions").CreateNewCondition(
                                                                            "PointLoadCondition3D1N",
                                                                            max_id,
                                                                            [node.Id],
                                                                            props)
            self.FEM_Solution.main_model_part.GetSubModelPart("computing_domain").AddCondition(cond)
            self.FEM_Solution.main_model_part.GetCondition(max_id).SetValue(KratosSMA.POINT_LOAD, [0.0,0.0,0.0])

#FindNeighboursIfNecessary===================================================================================================================================
    def FindNeighboursIfNecessary(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ComputeNeighboursIfNecessary")

        if self.domain_size == 3:
            KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.FEM_Solution.main_model_part).Execute()
            KratosMultiphysics.GenericFindElementalNeighboursProcess(self.FEM_Solution.main_model_part).Execute()
        else: # 2D
            KratosMultiphysics.GenericFindElementalNeighboursProcess(self.FEM_Solution.main_model_part).Execute()


#ComputeSkinSubModelPart============================================================================================================================
    def ComputeSkinSubModelPart(self):
        # Search the skin nodes for the remeshing
        if self.domain_size == 2:
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess2D(self.FEM_Solution.main_model_part,
                                                                               self.SkinDetectionProcessParameters)
        else: # 3D
            skin_detection_process = KratosMultiphysics.SkinDetectionProcess3D(self.FEM_Solution.main_model_part,
                                                                               self.SkinDetectionProcessParameters)
        skin_detection_process.Execute()

#ComputeDeltaTime============================================================================================================================
    def ComputeDeltaTime(self):
        if self.FEM_Solution.ProjectParameters["problem_data"].Has("time_step"):
            return self.FEM_Solution.ProjectParameters["problem_data"]["time_step"].GetDouble()

        elif self.FEM_Solution.ProjectParameters["problem_data"].Has("time_step_table"):

            current_time = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.TIME]

            tb = KratosMultiphysics.PiecewiseLinearTable()
            time_step_table = self.FEM_Solution.ProjectParameters["problem_data"]["time_step_table"].GetMatrix()
            for interval in range(time_step_table.Size1()):
                tb.AddRow(time_step_table[interval, 0], time_step_table[interval, 1])
            return tb.GetValue(current_time)
            raise Exception("::[MechanicalSolver]:: Time stepping not well defined!")
        else:
            raise Exception("::[MechanicalSolver]:: Time stepping not defined!")

#ExpandWetNodes============================================================================================================================
    def ExpandWetNodes(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ExpandWetNodes")

        if self.PressureLoad:
            # This must be called before Generating DEM
            self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] = 0 # It is modified inside
            extend_wet_nodes_process = KratosFemDem.ExpandWetNodesProcess(self.FEM_Solution.main_model_part)
            extend_wet_nodes_process.Execute()

#GenerateDEM============================================================================================================================
    def GenerateDEM(self): # This method creates the DEM elements and remove the damaged FEM, Additionally remove the isolated elements
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: GenerateDEM")

        if KratosFemDem.FEMDEMCouplingUtilities().IsGenerateDEMRequired(self.FEM_Solution.main_model_part):

            if self.PressureLoad:
                self.ExpandWetNodes()
                KratosFemDem.UpdatePressureVolumeProcess(self.FEM_Solution.main_model_part).Execute()
                self.ExpandWetNodes()

            self.UpdateDEMVariables()

            dem_generator_process = KratosFemDem.GenerateDemProcess(self.FEM_Solution.main_model_part, self.SpheresModelPart)
            dem_generator_process.Execute()

            # We remove the inactive DEM associated to fem_nodes
            element_eliminator = KratosMultiphysics.AuxiliarModelPartUtilities(self.FEM_Solution.main_model_part)
            element_eliminator.RemoveElementsAndBelongings(KratosMultiphysics.TO_ERASE)

            self.FindNeighboursIfNecessary()

            if self.domain_size == 3:
                # We assign the flag to recompute neighbours inside the 3D elements
                utils = KratosMultiphysics.VariableUtils()
                utils.SetNonHistoricalVariable(KratosFemDem.RECOMPUTE_NEIGHBOURS, True, self.FEM_Solution.main_model_part.Elements)

            # We update the skin for the DE-FE contact
            self.ComputeSkinSubModelPart()
            if self.DEMFEM_contact:
                self.TransferFEMSkinToDEM()

            # We reset the flag
            utils = KratosMultiphysics.VariableUtils()
            elements = self.FEM_Solution.main_model_part.Elements
            utils.SetNonHistoricalVariable(KratosFemDem.GENERATE_DEM, False, elements)

            self.ExtrapolatePressureLoad()

            if self.do_stabilization_solve:
                self.FEM_Solution.KratosPrintInfo("FEM-DEM:: Stabilization Calculation after removing FE...")
                self.FEM_Solution.solver.Solve()
                self.ExecuteAfterGeneratingDEM()


#ExtrapolatePressureLoad============================================================================================================================
    def ExtrapolatePressureLoad(self):
        if self.echo_level > 0:
            self.FEM_Solution.KratosPrintInfo("FEM-DEM:: ExtrapolatePressureLoad")

        if self.PressureLoad:
            # we reconstruct the pressure load if necessary
            if self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.RECONSTRUCT_PRESSURE_LOAD] == 1:
                self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] = 1
                while self.FEM_Solution.main_model_part.ProcessInfo[KratosFemDem.INTERNAL_PRESSURE_ITERATION] > 0:
                    if self.domain_size == 2:
                        KratosFemDem.ExtendPressureConditionProcess2D(self.FEM_Solution.main_model_part).Execute()
                    else:
                        KratosFemDem.ExtendPressureConditionProcess3D(self.FEM_Solution.main_model_part).Execute()

#UpdateDEMVariables============================================================================================================================
    def UpdateDEMVariables(self):
        update_de_kinematics_process = KratosFemDem.UpdateDemKinematicsProcess(self.FEM_Solution.main_model_part)
        update_de_kinematics_process.Execute()

#TransferNodalForcesToFEM============================================================================================================================
    def TransferNodalForcesToFEM(self):
        tranfer_nodal_forces_process = KratosFemDem.TransferNodalForcesToFem(self.FEM_Solution.main_model_part, False)
        tranfer_nodal_forces_process.Execute()

#WritePostListFile============================================================================================================================
    def WritePostListFile(self):
        pass
        # post_file_name = self.FEM_Solution.problem_name + ".post.lst"
        # time_label = round(self.FEM_Solution.step, 0)
        # PostListFile = open(post_file_name, "w")
        # PostListFile.write("Merge\n\n")
        # PostListFile.write(self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.res\n")
        # PostListFile.write(self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.msh\n")
        # PostListFile.write(os.path.join(self.FEM_Solution.problem_name + "_Post_Files", self.FEM_Solution.problem_name + "_" + str(time_label) + ".post.bin"))
        # PostListFile.close()

#InitializePlotsFiles============================================================================================================================
    def InitializePlotsFiles(self):
        # open general Displ/Reaction File
        if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size() != 0:
            self.PlotFile = open("PlotFile.txt","w")
            self.PlotFile.write("This File Plots the SUM of the displacement and reactions of the nodes selected in the lists!\n\n")
            if self.domain_size == 2:
                self.PlotFile.write("       time           displ_x        displ_y      Reaction_x     Reaction_y    \n")
            else:
                self.PlotFile.write("       time           displ_x        displ_y      displ_z        Reaction_x     Reaction_y     Reaction_z    \n")
            self.PlotFile.close()

        # self.PlotFileIter = open("iterations.txt","w")
        # self.PlotFileIter.write("This file prints the number of iterations at each time step\n\n")
        # self.PlotFileIter.write("       time           ITER\n")
        # self.PlotFileIter.close()

        self.TimePreviousPlotting = 0.0
        self.plot_files_nodes_list    = []
        self.plot_files_elements_list = []
        self.plot_files_nodes_id_list    = []
        self.plot_files_elements_id_list = []

        # open plots for nodes selected
        if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
            number_nodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
            for node in range(0, number_nodes):

                Id = self.FEM_Solution.ProjectParameters["watch_nodes_list"][node].GetInt()
                i_plot_file_node = open("PlotNode_" + str(Id) + ".txt","w")
                i_plot_file_node.write("\n")
                if self.domain_size == 2:
                    i_plot_file_node.write("       time          displ_x        displ_y         vel_x           vel_y         acc_x          acc_y        Reaction_x     Reaction_y    \n")
                else:
                    i_plot_file_node.write("       time          displ_x        displ_y        displ_z         vel_x           vel_y         vel_z           acc_x          acc_y          acc_z       Reaction_x     Reaction_y     Reaction_Z    \n")
                i_plot_file_node.close()
                self.plot_files_nodes_list.append(i_plot_file_node)
                self.plot_files_nodes_id_list.append(Id)

        # open plots for elements selected
        if self.FEM_Solution.ProjectParameters["watch_elements_list"].size() != 0:

            number_elems = self.FEM_Solution.ProjectParameters["watch_elements_list"].size()
            for elem in range(0, number_elems):
                Id = self.FEM_Solution.ProjectParameters["watch_elements_list"][elem].GetInt()
                i_plot_file_elem = open("PlotElement_" + str(Id) + ".txt","w")
                i_plot_file_elem.write("\n")
                if self.domain_size == 2:
                    i_plot_file_elem.write("          time                       Sxx                   Syy                      Sxy                    Exx                     Eyy                   Exy                Damage  \n")
                else:
                    i_plot_file_elem.write("       time             Sxx           Syy             Szz           Sxy            Syz            Sxz            Exx            Eyy            Ezz             Exy           Eyz            Exz          Damage  \n")
                i_plot_file_elem.close()
                self.plot_files_elements_list.append(i_plot_file_elem)
                self.plot_files_elements_id_list.append(Id)

#PrintPlotsFiles============================================================================================================================
    def PrintPlotsFiles(self):
        # Print the general file
        time = self.FEM_Solution.time
        total_reaction_x     = 0.0
        total_displacement_x = 0.0
        total_reaction_y     = 0.0
        total_displacement_y = 0.0
        total_displacement_z = 0.0
        total_reaction_z     = 0.0
        interval = self.FEM_Solution.ProjectParameters["interval_of_watching"].GetDouble()

        # self.PlotFileIter = open("iterations.txt", "a")
        # max_iter = self.FEM_Solution.ProjectParameters["solver_settings"]["max_iteration"].GetInt()
        # iterations = self.FEM_Solution.main_model_part.ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER]
        # if iterations < max_iter:
        #     self.PlotFileIter.write("    " + "{0:.4e}".format(time).rjust(11) + "        " + str(iterations) + "\n")
        # else:
        #     self.PlotFileIter.write("    " + "{0:.4e}".format(time).rjust(11) + "        " + str(iterations) + "  MAX iterations reached!" + "\n")
        # self.PlotFileIter.close()

        if self.FEM_Solution.time - self.TimePreviousPlotting >= interval:
            if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size() > 0:
                if self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][0].IsInt():
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
                        id_node = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetInt()
                        node = self.FEM_Solution.main_model_part.GetNode(id_node)
                        total_displacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                        total_displacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                        total_displacement_z += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)
                else:
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"].size()):
                        submodel_name = self.FEM_Solution.ProjectParameters["list_of_nodes_displacement"][index].GetString()
                        for node in self.FEM_Solution.main_model_part.GetSubModelPart(submodel_name).Nodes:
                            total_displacement_x += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
                            total_displacement_y += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y)
                            total_displacement_z += node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z)

                if self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][0].IsInt():
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):
                        id_node = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetInt()
                        node = self.FEM_Solution.main_model_part.GetNode(id_node)
                        total_reaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                        total_reaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                        total_reaction_z += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Z)
                else:
                    for index in range(0, self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"].size()):
                        submodel_name = self.FEM_Solution.ProjectParameters["list_of_nodes_reaction"][index].GetString()
                        for node in self.FEM_Solution.main_model_part.GetSubModelPart(submodel_name).Nodes:
                            total_reaction_x += node.GetSolutionStepValue(KratosMultiphysics.REACTION_X)
                            total_reaction_y += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Y)
                            total_reaction_z += node.GetSolutionStepValue(KratosMultiphysics.REACTION_Z)

                self.PlotFile = open("PlotFile.txt","a")
                if self.domain_size == 2:
                    self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(total_displacement_x).rjust(11) +
                                        "    " + "{0:.4e}".format(total_displacement_y).rjust(11) + "    " + "{0:.4e}".format(total_reaction_x).rjust(11) +
                                        "    " + "{0:.4e}".format(total_reaction_y).rjust(11) + "\n")
                else:
                    self.PlotFile.write("    " + "{0:.4e}".format(time).rjust(11) + "    " + "{0:.4e}".format(total_displacement_x).rjust(11) +
                        "    " + "{0:.4e}".format(total_displacement_y).rjust(11) + "    " + "{0:.4e}".format(total_displacement_z).rjust(11) +
                        "    " + "{0:.4e}".format(total_reaction_x).rjust(11) + "    " + "{0:.4e}".format(total_reaction_y).rjust(11) + "    " +
                        "{0:.4e}".format(total_reaction_z).rjust(11) + "\n")
                self.PlotFile.close()

            # Print the selected nodes files
            if self.FEM_Solution.ProjectParameters["watch_nodes_list"].size() != 0:
                NumNodes = self.FEM_Solution.ProjectParameters["watch_nodes_list"].size()
                for inode in range(0, NumNodes):
                    id_node = self.plot_files_nodes_id_list[inode]
                    node = self.FEM_Solution.main_model_part.GetNode(id_node)
                    self.plot_files_nodes_list[inode] = open("PlotNode_" + str(id_node) + ".txt","a")

                    displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
                    velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                    reaction = node.GetSolutionStepValue(KratosMultiphysics.REACTION)
                    acceleration = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION)

                    dx = displacement[0]
                    dy = displacement[1]
                    dz = displacement[2]
                    Rx = reaction[0]
                    Ry = reaction[1]
                    Rz = reaction[2]
                    vx = velocity[0]
                    vy = velocity[1]
                    vz = velocity[2]
                    ax = acceleration[0]
                    ay = acceleration[1]
                    az = acceleration[2]

                    if self.domain_size == 2:
                        self.plot_files_nodes_list[inode].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                            "{0:.4e}".format(dx).rjust(11) + "    " + "{0:.4e}".format(dy).rjust(11) + "    " +
                            "{0:.4e}".format(vx).rjust(11) + "    " + "{0:.4e}".format(vy).rjust(11) + "    " +
                            "{0:.4e}".format(ax).rjust(11) + "    " + "{0:.4e}".format(ay).rjust(11) + "    " +
                            "{0:.4e}".format(Rx).rjust(11) + "    " + "{0:.4e}".format(Ry).rjust(11) + "\n")
                    else:
                        self.plot_files_nodes_list[inode].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                            "{0:.4e}".format(dx).rjust(11) + "    " + "{0:.4e}".format(dy).rjust(11) + "    " + "{0:.4e}".format(dz).rjust(11) + "    " +
                            "{0:.4e}".format(vx).rjust(11) + "    " + "{0:.4e}".format(vy).rjust(11) + "    " + "{0:.4e}".format(vz).rjust(11) + "    " +
                            "{0:.4e}".format(ax).rjust(11) + "    " + "{0:.4e}".format(ay).rjust(11) + "    " + "{0:.4e}".format(az).rjust(11) + "    " +
                            "{0:.4e}".format(Rx).rjust(11) + "    " + "{0:.4e}".format(Ry).rjust(11) + "    " + "{0:.4e}".format(Rz).rjust(11) + "\n")

                    self.plot_files_nodes_list[inode].close()

            # print the selected element files
            if self.FEM_Solution.ProjectParameters["watch_elements_list"].size() != 0:
                NumElem = self.FEM_Solution.ProjectParameters["watch_elements_list"].size()
                for iElem in range(0, NumElem):
                    Idelem = self.PlotFilesElementsIdList[iElem]
                    Elem = self.FEM_Solution.main_model_part.GetElement(Idelem)
                    self.plot_files_elements_list[iElem] = open("PlotElement_" + str(Idelem) + ".txt","a")

                    stress_tensor = Elem.CalculateOnIntegrationPoints(KratosFemDem.STRESS_VECTOR_INTEGRATED, self.FEM_Solution.main_model_part.ProcessInfo)
                    strain_vector = Elem.GetValue(KratosFemDem.FEMDEM_STRAIN_VECTOR)

                    damage = Elem.GetValue(KratosFemDem.DAMAGE_ELEMENT)

                    if self.domain_size == 2:
                        Sxx = stress_tensor[0][0]
                        Syy = stress_tensor[0][1]
                        Sxy = stress_tensor[0][2]
                        Exx = strain_vector[0]
                        Eyy = strain_vector[1]
                        Exy = strain_vector[2]
                        self.plot_files_elements_list[iElem].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                            "{0:.4e}".format(Sxx).rjust(11) + "    " + "{0:.4e}".format(Syy).rjust(11) + "    " +
                            "{0:.4e}".format(Sxy).rjust(11) + "    " + "{0:.4e}".format(Exx).rjust(11) +
                            "    " + "{0:.4e}".format(Eyy).rjust(11) + "    " + "{0:.4e}".format(Exy).rjust(11) +
                            "   " + "{0:.4e}".format(damage).rjust(11) + "\n")
                    else:
                        Sxx = stress_tensor[0][0]
                        Syy = stress_tensor[0][1]
                        Szz = stress_tensor[0][2]
                        Sxy = stress_tensor[0][3]
                        Syz = stress_tensor[0][4]
                        Sxz = stress_tensor[0][5]
                        Exx = strain_tensor[0]
                        Eyy = strain_tensor[1]
                        Ezz = strain_tensor[2]
                        Exy = strain_tensor[3]
                        Eyz = strain_tensor[4]
                        Exz = strain_tensor[5]
                        self.plot_files_elements_list[iElem].write("    " + "{0:.4e}".format(time).rjust(11) + "    " +
                        "{0:.4e}".format(Sxx).rjust(11) + "    " + "{0:.4e}".format(Syy).rjust(11) + "    " +
                        "{0:.4e}".format(Szz).rjust(11) + "    " + "{0:.4e}".format(Sxy).rjust(11) + "    " +
                        "{0:.4e}".format(Syz).rjust(11) + "    " + "{0:.4e}".format(Sxz).rjust(11) + "    " +
                        "{0:.4e}".format(Exx).rjust(11) +
                        "    " + "{0:.4e}".format(Eyy).rjust(11) + "    " + "{0:.4e}".format(Ezz).rjust(11) +
                        "    " + "{0:.4e}".format(Exy).rjust(11) + "    " + "{0:.4e}".format(Eyz).rjust(11) +
                        "    " + "{0:.4e}".format(Exz).rjust(11) +
                        "   "  + "{0:.4e}".format(damage).rjust(11) + "\n")

                    self.plot_files_elements_list[iElem].close()
            self.TimePreviousPlotting = time

#GetMaximumConditionId============================================================================================================================
    def GetMaximumConditionId(self):
        max_id = 0
        for condition in self.FEM_Solution.main_model_part.Conditions:
            if condition.Id > max_id:
                max_id = condition.Id
        return max_id

#PrintResults============================================================================================================================
    def PrintResults(self):

        print_parameters = self.FEM_Solution.ProjectParameters["output_configuration"]["result_file_configuration"]
        if self.FEM_Solution.step == 1: # always print the 1st step
            self.gid_output.Writeresults(self.FEM_Solution.time)
            self.FEM_Solution.time_old_print = self.FEM_Solution.time
            self.FEM_Solution.step_old_print = self.FEM_Solution.step
        else:
            time_to_print = 0
            if print_parameters["output_control_type"].GetString() == "step":
                time_to_print = self.FEM_Solution.step - self.FEM_Solution.step_old_print
            else:
                time_to_print = self.FEM_Solution.time - self.FEM_Solution.time_old_print

            if print_parameters["output_control_type"].GetString() == "step":
                if print_parameters["output_interval"].GetInt() - time_to_print == 0:
                    self.gid_output.Writeresults(self.FEM_Solution.time)
                    self.FEM_Solution.step_old_print = self.FEM_Solution.step
            else:
                if print_parameters["output_interval"].GetDouble() - time_to_print < 1e-2 * self.FEM_Solution.delta_time:
                    self.gid_output.Writeresults(self.FEM_Solution.time)
                    self.FEM_Solution.time_old_print = self.FEM_Solution.time

#ExecuteBeforeGeneratingDEM============================================================================================================================
    def ExecuteBeforeGeneratingDEM(self):
        """Here the erased are labeled as INACTIVE so you can access to them. After calling
           GenerateDEM they are totally erased """
        pass

#ExecuteAfterGeneratingDEM============================================================================================================================
    def ExecuteAfterGeneratingDEM(self):
        # self.ExtrapolatePressureLoad()
        self.SpheresModelPart = self.ParticleCreatorDestructor.GetSpheresModelPart()
        # We update coordinates, displ and velocities of the DEM according to FEM
        self.UpdateDEMVariables()

#BeforeSolveDEMOperations============================================================================================================================
    def BeforeSolveDEMOperations(self):
        self.DEM_Solution.time = self.FEM_Solution.time
        self.DEM_Solution.step = self.FEM_Solution.step
        # self.DEM_Solution.UpdateTimeInModelParts()

#TransferFEMSkinToDEM============================================================================================================================
    def TransferFEMSkinToDEM(self):
        fem_skin_mp = self.FEM_Solution.main_model_part.GetSubModelPart("SkinDEMModelPart")

        if self.DEM_Solution.rigid_face_model_part.HasSubModelPart("SkinTransferredFromStructure"):
            self.EraseConditionsAndNodesSubModelPart()
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.GetSubModelPart("SkinTransferredFromStructure")
            props = self.DEM_Solution.spheres_model_part.GetProperties()[2]
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)
        else: # have to create it
            # props = self.CreateFEMPropertiesForDEFEContact()
            props = self.DEM_Solution.spheres_model_part.GetProperties()[2]
            dem_walls_mp = self.DEM_Solution.rigid_face_model_part.CreateSubModelPart("SkinTransferredFromStructure")
            # dem_walls_mp.AddProperties(props)
            DemFem.DemStructuresCouplingUtilities().TransferStructuresSkinToDem(fem_skin_mp, dem_walls_mp, props)

    #-----------------------------------
    def EraseConditionsAndNodesSubModelPart(self):
        DEM_sub_model_part = self.DEM_Solution.rigid_face_model_part.GetSubModelPart("SkinTransferredFromStructure")
        self.DEM_Solution.rigid_face_model_part.Conditions.clear()
        self.DEM_Solution.rigid_face_model_part.Nodes.clear()


#============================================================================================================================

    def InitializePostProcess(self):
        mixed_fluid_solid_mp       = self.model.CreateModelPart('mixed_fluid_solid_mp')
        mixed_fluid_solid_balls_mp = self.model.CreateModelPart('mixed_fluid_solid_balls_mp')
        mixed_solid_balls_mp       = self.model.CreateModelPart('mixed_solid_balls_mp')
        dummy_fluid_part           = self.model.CreateModelPart('dummy_fluid_part')

        filename = os.path.join(self.DEM_Solution.post_path, self.DEM_Solution.DEM_parameters["problem_name"].GetString())
        self.gid_output = gid_output.FemDemCoupledGiDOutput(
                            filename,
                            True,
                            "Binary",
                            "Multiples",
                            True,
                            True,
                            self.FEM_Solution.main_model_part,
                            dummy_fluid_part,
                            self.DEM_Solution.spheres_model_part,
                            self.DEM_Solution.cluster_model_part,
                            self.DEM_Solution.rigid_face_model_part,
                            mixed_fluid_solid_mp,
                            mixed_solid_balls_mp,
                            mixed_fluid_solid_balls_mp)

        solid_nodal_results = ["DISPLACEMENT", "ACCELERATION", "VELOCITY"]
        dem_nodal_results = ["TOTAL_FORCES", "RADIUS"]
        fluid_nodal_results = []
        clusters_nodal_results = []
        rigid_faces_nodal_results = []
        mixed_solid_fluid_nodal_results = []
        mixed_solid_balls_nodal_results = []
        mixed_solid_balls_fluid_nodal_results = []

        gp_list = self.FEM_Solution.ProjectParameters["output_configuration"]["result_file_configuration"]["gauss_point_results"]
        gauss_points_results = []
        for i in gp_list.values():
            gauss_points_results.append(i.GetString())

        self.gid_output.initialize_dem_fem_results(solid_nodal_results,
                                                   fluid_nodal_results,
                                                   dem_nodal_results,
                                                   clusters_nodal_results,
                                                   rigid_faces_nodal_results,
                                                   mixed_solid_fluid_nodal_results,
                                                   mixed_solid_balls_nodal_results,
                                                   mixed_solid_balls_fluid_nodal_results,
                                                   gauss_points_results)