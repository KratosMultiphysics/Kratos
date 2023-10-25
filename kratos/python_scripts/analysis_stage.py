# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
import numpy as np

class AnalysisStage(object):
    """The base class for the AnalysisStage-classes in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, project_parameters):
        """The constructor of the AnalysisStage-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedAnalysis, self).__init__(project_parameters)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        project_parameters -- The ProjectParameters used
        """
        if not isinstance(model, KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if not isinstance(project_parameters, KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.project_parameters = project_parameters

        ## Get echo level and parallel type
        self.echo_level = self.project_parameters["problem_data"]["echo_level"].GetInt()
        self.parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()
        is_distributed_run = KratosMultiphysics.IsDistributedRun()

        if self.parallel_type == "OpenMP" and is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"OpenMP" is specified as "parallel_type", but Kratos is running distributed!')
        if self.parallel_type == "MPI" and not is_distributed_run:
            KratosMultiphysics.Logger.PrintWarning("Parallel Type", '"MPI" is specified as "parallel_type", but Kratos is not running distributed!')

        self._GetSolver().AddVariables() # this creates the solver and adds the variables

    def Run(self):
        """This function executes the entire AnalysisStage
        It can be overridden by derived classes
        """
        self.Initialize()
        self.RunSolutionLoop()
        self.Finalize()

    def KeepAdvancingSolutionLoop(self):
        """This function specifies the stopping criteria for breaking the solution loop
        It can be overridden by derived classes
        """
        return self.time < self.end_time

    ######################################################################################################
    def ImposeInitialStateProcess(self):
        imposed_strain = KratosMultiphysics.Vector(3)
        imposed_stress = KratosMultiphysics.Vector(3)
        imposed_def_grad = KratosMultiphysics.Matrix(2,2)
        imposed_strain[0] = 0.0
        imposed_strain[1] = 0.0
        imposed_strain[2] = 0.0
        imposed_stress[0] = 0.0
        imposed_stress[1] = 0.0
        imposed_stress[2] = 0.0
        imposed_def_grad[0,0] = 0.0
        imposed_def_grad[0,1] = 0.0
        imposed_def_grad[1,0] = 0.0
        imposed_def_grad[1,1] = 0.0

        # create process
        KratosMultiphysics.SetInitialStateProcess2D(self._GetSolver().GetComputingModelPart(),
                                                    imposed_strain,
                                                    imposed_stress,
                                                    imposed_def_grad).ExecuteInitializeSolutionStep()
    ######################################################################################################
    
    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        # initialize InitialState
        self.ImposeInitialStateProcess()
        # initialization and process call placed in separate method, see above method

        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()     
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()
        # export DISPLACEMENTS
        self.WriteOutputDisplacements()

    def WriteOutputDisplacements(self):
        displacements = []
        finalconfig_X = []
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            node_displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            displacements.append(node_displacement)
            node_coordinate = [node.X, node.Y, node.Z]
            finalconfig_X.append(node_coordinate)
        np.save("DISPLACEMENTS_U" + str(self.case), displacements)
        np.savetxt("DISPLACEMENTS_U" + str(self.case) + ".txt", displacements)
        np.save("COORDINATES_X" + str(self.case), finalconfig_X)
        np.savetxt("COORDINATES_X" + str(self.case) + ".txt", finalconfig_X)

    def Initialize(self):
        """This function initializes the AnalysisStage
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        This function has to be implemented in deriving classes!
        """
        # Modelers:
        self._CreateModelers()
        self._ModelersSetupGeometryModel()
        self._ModelersPrepareGeometryModel()
        self._ModelersSetupModelPart()

        self._GetSolver().ImportModelPart()
        self._GetSolver().PrepareModelPart()
        self._GetSolver().AddDofs()

        self.ModifyInitialProperties()
        self.ModifyInitialGeometry()

        ##here we initialize user-provided processes
        self.__CreateListOfProcesses() # has to be done after importing and preparing the ModelPart
        for process in self._GetListOfProcesses():
            process.ExecuteInitialize()

        self._GetSolver().Initialize()
        self.Check()

        self.ModifyAfterSolverInitialize()

        for process in self._GetListOfProcesses():
            process.ExecuteBeforeSolutionLoop()

        ## Stepping and time settings
        self.end_time = self.project_parameters["problem_data"]["end_time"].GetDouble()

        if self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self.time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        else:
            self.time = self.project_parameters["problem_data"]["start_time"].GetDouble()
            self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME] = self.time

        ## If the echo level is high enough, print the complete list of settings used to run the simulation
        if self.echo_level > 1:
            with open("ProjectParametersOutput.json", 'w') as parameter_output_file:
                parameter_output_file.write(self.project_parameters.PrettyPrintJsonString())

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -START- ")

    def Finalize(self):
        """This function finalizes the AnalysisStage
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        for process in self._GetListOfProcesses():
            process.ExecuteFinalize()

        self._GetSolver().Finalize()

        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Analysis -END- ")

    def GetFinalData(self):
        """Returns the final data dictionary.

        The main purpose of this function is to retrieve any data (in a key-value format) from outside the stage.
        Note that even though it can be called at any point, it is intended to be called at the end of the stage run.
        """

        return {}

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        self.PrintAnalysisStageProgressInformation()

        self.ApplyBoundaryConditions() #here the processes are called
        self.ChangeMaterialProperties() #this is normally empty
        self._GetSolver().InitializeSolutionStep()


    def PrintAnalysisStageProgressInformation(self):
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "STEP: ", self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP])
        KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "TIME: ", self.time)

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        self._GetSolver().FinalizeSolutionStep()

        for process in self._GetListOfProcesses():
            process.ExecuteFinalizeSolutionStep()

    def OutputSolutionStep(self):
        """This function printed / writes output files after the solution of a step
        """
        execute_was_called = False
        for output_process in self._GetListOfOutputProcesses():
            if output_process.IsOutputStep():
                if not execute_was_called:
                    for process in self._GetListOfProcesses():
                        process.ExecuteBeforeOutputStep()
                    execute_was_called = True

                output_process.PrintOutput()

        if execute_was_called:
            for process in self._GetListOfProcesses():
                process.ExecuteAfterOutputStep()

    def Check(self):
        """This function checks the AnalysisStage
        """
        self._GetSolver().Check()
        for process in self._GetListOfProcesses():
            process.Check()

    def Clear(self):
        """This function clears the AnalysisStage
        """
        self._GetSolver().Clear()

    def ModifyInitialProperties(self):
        """this is the place to eventually modify material properties in the stage """
        pass

    def ModifyInitialGeometry(self):
        """this is the place to eventually modify geometry (for example moving nodes) in the stage """
        ######################################################################################################
        # save deformed configuration nodal positions
        var_utils = KratosMultiphysics.VariableUtils()
        self.deformed_config_coordinates = var_utils.GetCurrentPositionsVector(self._GetSolver().GetComputingModelPart().Nodes, 3)
        self.case = 1
        if self.case == 2:
            self.UpdateMyGeometry("DISPLACEMENTS_U0.txt", True)
            self.UpdateMyGeometry("DISPLACEMENTS_U1.txt", True)
        # print("original coordinates:\n", self.deformed_config_coordinates)
        # for node in self._GetSolver().GetComputingModelPart().Nodes:
        #     print(node.Id, node.X, node.Y, node.Z)
        # new_coordinates = np.load('COORDINATES_X1.npy')
        # new_coord_vector = np.reshape(new_coordinates, -1)
        # print("new coordinates:\n", new_coordinates)
        # print("new coordinates (vector):\n", new_coord_vector)
        # newcoord_KratosVector = KratosMultiphysics.Vector(new_coord_vector)
        # print("kratos vector:\n", newcoord_KratosVector)
        # var_utils.SetCurrentPositionsVector(self._GetSolver().GetComputingModelPart().Nodes, newcoord_KratosVector)
        # self.deformed_config_coordinates = var_utils.GetCurrentPositionsVector(self._GetSolver().GetComputingModelPart().Nodes, 3)
        # print("updated coordinates:\n", self.deformed_config_coordinates)
        # for node in self._GetSolver().GetComputingModelPart().Nodes:
        #     print(node.Id, node.X, node.Y, node.Z)
        # # test GetPositions functions 
        # initial_positions_def = self.initial_unmodified_coordinates
        # current_positions_def = var_utils.GetCurrentPositionsVector(self._GetSolver().GetComputingModelPart().Nodes, 3)
        # print("initial positions (deformed)")
        # print(initial_positions_def)
        # print("current positions (deformed)")
        # print(current_positions_def)
    
        # Calculate Jacobians from deformed configuration
        self.deformed_config_jacobians = []
        for elem in self._GetSolver().GetComputingModelPart().Elements:
            H0 = elem.Properties[KratosMultiphysics.THICKNESS]
            # TODO: See if numpy objects can be replaced all by Kratos-type elements. Check KratosMultiphysics.Matrix available methods
            # idea: instead of list, read length of ModelPart.Elements and create Kratos Vector of this length, add J here
            extracolumn = np.array([[0.0, 0.0, H0/2.0]])
            elm_Jacobian = np.array(elem.GetGeometry().Jacobian(0))
            J_X = np.concatenate((elm_Jacobian, extracolumn.T), axis=1)

            local_axis_1 = elem.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, self._GetSolver().GetComputingModelPart().ProcessInfo)
            local_axis_2 = elem.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_2, self._GetSolver().GetComputingModelPart().ProcessInfo)
            local_axis_3 = elem.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_3, self._GetSolver().GetComputingModelPart().ProcessInfo)
            # print("local axes DEFORMED\n")
            # print("elem", elem.Id,":")
            # print(local_axis_1[0])
            # print(local_axis_2[0])
            # print(local_axis_3[0])
            
            T_elm = KratosMultiphysics.Matrix(3, 3)
            T_elm[0, 0] = local_axis_1[0][0]
            T_elm[0, 1] = local_axis_1[0][1]
            T_elm[0, 2] = local_axis_1[0][2]
            T_elm[1, 0] = local_axis_2[0][0]
            T_elm[1, 1] = local_axis_2[0][1]
            T_elm[1, 2] = local_axis_2[0][2]
            T_elm[2, 0] = local_axis_3[0][0]
            T_elm[2, 1] = local_axis_3[0][1]
            T_elm[2, 2] = local_axis_3[0][2]

            J_X = T_elm.transpose() @ J_X
            self.deformed_config_jacobians.append(J_X)
        # export Jacobians J
        np.save("jacobians_J_" + str(self.case), self.deformed_config_jacobians)
        J0txt = np.reshape(self.deformed_config_jacobians, (-1,3))
        np.savetxt("jacobians_J_" + str(self.case) + ".txt", J0txt)

        # Calculate normals from deformed configuration
        print("\n ::TESTING:: START Calculate normals \n")
        self.deformed_config_normals = []
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()
        normal_calculation_utils.CalculateUnitNormalsNonHistorical(self._GetSolver().GetComputingModelPart(), 0)
    
        print("\n ::TESTING:: Store deformed configuration nodal normals \n")
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            self.deformed_config_normals.append(np.array(node.GetValue(KratosMultiphysics.NORMAL)))
            # print(node.Id, self.deformed_config_normals[node.Id - 1])
        print("\n ::TESTING:: FINISH Calculate normals \n")

        # go back to X0 to calculate J0
        if self.case == 1:
            self.UpdateMyGeometry("DISPLACEMENTS_U0.txt", True)
        if self.case == 2:
            self.UpdateMyGeometry("DISPLACEMENTS_U1.txt", False)
        
        # save flat configuration nodal positions
        self.flattened_coordinates = var_utils.GetInitialPositionsVector(self._GetSolver().GetComputingModelPart().Nodes, 3)
        self.initial_displacements = self.deformed_config_coordinates - self.flattened_coordinates
        # print("initial_displacements:", self.initial_displacements)

        # test GetPositions functions 
        initial_positions_flat = self.flattened_coordinates
        current_positions_flat = var_utils.GetCurrentPositionsVector(self._GetSolver().GetComputingModelPart().Nodes, 3)
        print("initial positions (flat)")
        print(initial_positions_flat)
        print("current positions (flat)")
        print(current_positions_flat)

        # assign initial strains
        self.InitializeMyStrains()
        ######################################################################################################

    def UpdateMyGeometry(self, displacements, boolval):
        # read prescribed displacements
        u_prescribed = np.loadtxt(displacements)
        i_disp = 0
        # stretch flat geometry to induce strains
        if boolval == True:
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                node.X0 += u_prescribed[i_disp,0]
                node.Y0 += u_prescribed[i_disp,1]
                node.Z0 += u_prescribed[i_disp,2]
                node.X += u_prescribed[i_disp,0]
                node.Y += u_prescribed[i_disp,1]
                node.Z += u_prescribed[i_disp,2]
                print(node.Id, node.X0, node.Y0, node.Z0)
                print(node.Id, node.X, node.Y, node.Z)
                i_disp += 1
        if boolval == False:
            for node in self._GetSolver().GetComputingModelPart().Nodes:
                node.X0 -= u_prescribed[i_disp,0]
                node.Y0 -= u_prescribed[i_disp,1]
                node.Z0 -= u_prescribed[i_disp,2]
                node.X -= u_prescribed[i_disp,0]
                node.Y -= u_prescribed[i_disp,1]
                node.Z -= u_prescribed[i_disp,2]
                print(node.Id, node.X0, node.Y0, node.Z0)
                print(node.Id, node.X, node.Y, node.Z)
                i_disp += 1
    
    def InitializeMyStrains(self):
        self.flat_config_jacobians = []
        for element in self._GetSolver().GetComputingModelPart().Elements:
            H0 = element.Properties[KratosMultiphysics.THICKNESS]
            extracolumn = np.array([[0.0, 0.0, H0/2.0]])
            elm_Jacobian = np.array(element.GetGeometry().Jacobian(0))
            J_X0 = np.concatenate((elm_Jacobian, extracolumn.T), axis=1)
            self.flat_config_jacobians.append(J_X0)
        # export Jacobians J0
        np.save("jacobians_J0_" + str(self.case), self.flat_config_jacobians)
        Jtxt = np.reshape(self.flat_config_jacobians, (-1,3))
        np.savetxt("jacobians_J0_" + str(self.case) + ".txt", Jtxt)
            
        F_list = []
        strains_list = []
        # Read pre-calculated strains
        # strains_def_config = np.load('strains_cylinder.npy')
        for J, J0, element in zip(self.deformed_config_jacobians,self.flat_config_jacobians, self._GetSolver().GetComputingModelPart().Elements):

            local_axis_1 = element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_1, self._GetSolver().GetComputingModelPart().ProcessInfo)
            local_axis_2 = element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_2, self._GetSolver().GetComputingModelPart().ProcessInfo)
            local_axis_3 = element.CalculateOnIntegrationPoints(KratosMultiphysics.LOCAL_AXIS_3, self._GetSolver().GetComputingModelPart().ProcessInfo)
            # print("local axes FLAT\n")
            # print("elem", element.Id,":")
            # print(local_axis_1[0])
            # print(local_axis_2[0])
            # print(local_axis_3[0])

            T0_elm = KratosMultiphysics.Matrix(3, 3)
            T0_elm[0, 0] = local_axis_1[0][0]
            T0_elm[0, 1] = local_axis_1[0][1]
            T0_elm[0, 2] = local_axis_1[0][2]
            T0_elm[1, 0] = local_axis_2[0][0]
            T0_elm[1, 1] = local_axis_2[0][1]
            T0_elm[1, 2] = local_axis_2[0][2]
            T0_elm[2, 0] = local_axis_3[0][0]
            T0_elm[2, 1] = local_axis_3[0][1]
            T0_elm[2, 2] = local_axis_3[0][2]

            J0 = T0_elm.transpose() @ J0
            J0_inv = np.linalg.inv(J0)

            # Calculate strain by Jacobians
            JTJ = J.T @ J
            G_hat = J0_inv
            g_hat = JTJ
            C3D = G_hat.T @ (g_hat @ G_hat)
            E3D = 0.5 * (C3D - np.eye(3))

            # Calculate strain by Deformation Gradient F
            defgrad = J @ J0_inv
            C3D_FFT = defgrad.T @ defgrad
            F_list.append(defgrad)
            E3D_FFT = 0.5 * (C3D_FFT - np.eye(3))
            
            # Check matrix operations
            print("F[", element.Id, "]:\n", defgrad)
            # print("C3D:\n", C3D)
            # print("F*FT:\n", C3D_FFT)
            
            # TODO: CHECK STRAINS

            E_voigt = KratosMultiphysics.Vector(3)
            E_voigt[0] = E3D_FFT[0, 0]
            E_voigt[1] = E3D_FFT[1, 1]
            E_voigt[2] = (E3D_FFT[0, 1] + E3D_FFT[1, 0])

            # # selected elements
            # E_voigt = KratosMultiphysics.Vector(3)
            # E_voigt[0] = 0.0
            # E_voigt[1] = 0.0
            # E_voigt[2] = 0.0
            # if element.Id == 1 or element.Id == 2 or element.Id == 21 or element.Id == 22:
            #     E_voigt[0] = 0.01
            #     E_voigt[1] = 0.01
            #     E_voigt[2] = 0.0
            # if element.Id == 15 or element.Id == 19 or element.Id == 39 or element.Id == 40:
            #     E_voigt[0] = 0.01
            #     E_voigt[1] = 0.01
            #     E_voigt[2] = 0.0

            # apply pre-calculated strains from original deformed configuration to flat input
            element.SetValue(KratosMultiphysics.INITIAL_STRAIN_VECTOR, E_voigt)
            strainvalue = np.array(E_voigt)
            print("E[", element.Id, "]:\n", strainvalue)
            strains_list.append(strainvalue)
            # TODO later - Create condition for PointLoad, reference command:
            # self._GetSolver().GetComputingModelPart().CreateCondition()
        # print("strain list:\n", strains_list)
        # export F, E
        np.save("F" + str(self.case) + "_defgrad", F_list)
        Ftxt = np.reshape(F_list, (-1,3))
        np.savetxt("F" + str(self.case) + "_defgrad.txt", Ftxt)
        np.save("E" + str(self.case) + "_strains", strains_list)
        np.savetxt("E" + str(self.case) + "_strains.txt", strains_list)


    def ModifyAfterSolverInitialize(self):
        """this is the place to eventually do any modification that requires the solver to be initialized """
        pass

    def ApplyBoundaryConditions(self):
        """here the boundary conditions is applied, by calling the InitializeSolutionStep function of the processes"""

        for process in self._GetListOfProcesses():
            process.ExecuteInitializeSolutionStep()

        #other operations as needed

    def ChangeMaterialProperties(self):
        """this function is where the user could change material parameters as a part of the solution step """
        pass

    def Save(self, serializer: KratosMultiphysics.StreamSerializer) -> None:
        """Serializes current analysis stage instance

        This method is intended to make the class pure Python (pickable). This means serialize all the Kratos objects,
        that is to say all the objects coming from Pybind, with the provided serializer. After the serialization, it is
        required to assign None value to all the objects in order to make the class pickable.
        """
        pass

    def Load(self, serializer: KratosMultiphysics.StreamSerializer) -> None:
        """Loads current analysis stage instance

        From the given serializer, this method restores current class from a pure Python status (pickable) to the one in the serializer.
        """
        pass

    def _GetSolver(self):
        if not hasattr(self, '_solver'):
            self._solver = self._CreateSolver()
        return self._solver

    def _CreateSolver(self):
        """Create the solver
        """
        raise Exception("Creation of the solver must be implemented in the derived class.")

    def _AdvanceTime(self):
        """ Computes the following time
            The default method simply calls the solver
        """
        return self._GetSolver().AdvanceInTime(self.time)

    ### Modelers
    def _ModelersSetupGeometryModel(self):
        # Import or generate geometry models from external input.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model started.")
            modeler.SetupGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup Geometry Model finished.")

    def _ModelersPrepareGeometryModel(self):
        # Prepare or update the geometry model_part.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model started.")
            modeler.PrepareGeometryModel()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Prepare Geometry Model finished.")

    def _ModelersSetupModelPart(self):
        # Convert the geometry model or import analysis suitable models.
        for modeler in self._GetListOfModelers():
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart started.")
            modeler.SetupModelPart()
            if self.echo_level > 1:
                KratosMultiphysics.Logger.PrintInfo(self._GetSimulationName(), "Modeler: ", str(modeler), " Setup ModelPart finished.")

    ### Modelers
    def _GetListOfModelers(self):
        """ This function returns the list of modelers
        """
        if not hasattr(self, '_list_of_modelers'):
            raise Exception("The list of modelers was not yet created!")
        return self._list_of_modelers

    def _CreateModelers(self):
        """ List of modelers in following format:
        "modelers" : [{
            "modeler_name" : "geometry_import",
            "Parameters" : {
                "echo_level" : 0,
                // settings for this modeler
            }
        },{ ... }]
        """
        self._list_of_modelers = []

        if self.project_parameters.Has("modelers"):
            from KratosMultiphysics.modeler_factory import KratosModelerFactory
            factory = KratosModelerFactory()

            modelers_list = self.project_parameters["modelers"]
            self._list_of_modelers = factory.ConstructListOfModelers(self.model, modelers_list)

    ### Processes
    def _GetListOfProcesses(self):
        """This function returns the list of processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_processes'):
            raise Exception("The list of processes was not yet created!")
        return self._list_of_processes

    def _GetListOfOutputProcesses(self):
        """This function returns the list of output processes involved in this Analysis
        """
        if not hasattr(self, '_list_of_output_processes'):
            raise Exception("The list of output-processes was not yet created!")
        return self._list_of_output_processes

    def _CreateProcesses(self, parameter_name, initialization_order):
        """Create a list of processes
        Format:
        "processes" : {
            initial_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ],
            boundary_processes : [
                { proces_specific_params },
                { proces_specific_params }
            ]
        }
        The order of initialization can be specified by setting it in "initialization_order"
        if e.g. the "boundary_processes" should be constructed before the "initial_processes", then
        initialization_order should be a list containing ["boundary_processes", "initial_processes"]
        see the functions _GetOrderOfProcessesInitialization and _GetOrderOfOutputProcessesInitialization
        """
        list_of_processes = []

        factory = KratosProcessFactory(self.model)

        if self.project_parameters.Has(parameter_name):
            processes_params = self.project_parameters[parameter_name]

            # first initialize the processes that depend on the order
            for processes_names in initialization_order:
                if processes_params.Has(processes_names):
                    list_of_processes += factory.ConstructListOfProcesses(processes_params[processes_names])

            # then initialize the processes that don't depend on the order
            for name, value in processes_params.items():
                if not name in initialization_order:
                    list_of_processes += factory.ConstructListOfProcesses(value) # Does this work? or should it be processes[name]

        return list_of_processes

    def _GetOrderOfProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the processes matters
        """
        return []

    def _GetOrderOfOutputProcessesInitialization(self):
        """This function can be overridden in derived classes if the order of
        initialization for the output-processes matters
        """
        return []

    def _CheckDeprecatedOutputProcesses(self, list_of_processes):
        deprecated_output_processes = []
        for process in list_of_processes:
            if issubclass(type(process), KratosMultiphysics.OutputProcess):
                deprecated_output_processes.append(process)
                msg  = "{} is an OutputProcess. However, it has been constructed as a regular process.\n"
                msg += "Please, define it as an 'output_processes' in the ProjectParameters."
                IssueDeprecationWarning("AnalysisStage", msg.format(process.__class__.__name__))
        return deprecated_output_processes

    def _GetSimulationName(self):
        """Returns the name of the Simulation
        """
        return "Analysis"

    def __CreateListOfProcesses(self):
        """This function creates the processes and the output-processes
        """
        order_processes_initialization = self._GetOrderOfProcessesInitialization()
        self._list_of_processes        = self._CreateProcesses("processes", order_processes_initialization)
        deprecated_output_processes    = self._CheckDeprecatedOutputProcesses(self._list_of_processes)
        order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
        self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
        self._list_of_processes.extend(self._list_of_output_processes) # Adding the output processes to the regular processes
        self._list_of_output_processes.extend(deprecated_output_processes)

    def __CheckIfSolveSolutionStepReturnsAValue(self, is_converged):
        """In case the solver does not return the state of convergence
        (same as the SolvingStrategy does) then issue ONCE a deprecation-warning

        """
        if is_converged is None:
            if not hasattr(self, '_map_ret_val_depr_warnings'):
                self._map_ret_val_depr_warnings = []
            solver_class_name = self._GetSolver().__class__.__name__
            # used to only print the deprecation-warning once
            if not solver_class_name in self._map_ret_val_depr_warnings:
                self._map_ret_val_depr_warnings.append(solver_class_name)
                warn_msg  = 'Solver "{}" does not return '.format(solver_class_name)
                warn_msg += 'the state of convergence from "SolveSolutionStep"'
                IssueDeprecationWarning("AnalysisStage", warn_msg)
