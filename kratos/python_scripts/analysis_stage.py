# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

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
    
    # INITIALIZE FOR INVERSE FORMING. CLEAN UP
    # def Initialize(self):
    #     for elm in self._GetSolver().GetComputingModelPart().Elements:
    #         for node in elm.GetNodes():
    #             node.X0 = node.X
    #             node.Y0 = node.Y
    #             node.Z0 = 0.0
    #             print(node.Id, node.X0, node.Y0, node.Z0)       

    #     # super.().Initialize()

    def InvertMatrix(self,A):
        Ainv = KratosMultiphysics.Matrix(2,2);
        det = A[0,0]*A[1,1]-A[1,0]*A[0,1]

        Ainv[0,0] = A[1,1]/det
        Ainv[0,1] = -A[0,1]/det
        Ainv[1,0] = -A[1,0]/det
        Ainv[1,1] = A[0,0]/det
        return Ainv
    
    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """

        print("/n ::TESTING:: START Calculate normals /n")
        normal_calculation_utils = KratosMultiphysics.NormalCalculationUtils()
        normal_calculation_utils.CalculateUnitNormalsNonHistorical(self._GetSolver().GetComputingModelPart(), 0)
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            normal = node.GetValue(KratosMultiphysics.NORMAL)
            # print(node.Id, normal)
        print("/n ::TESTING:: FINISH Calculate normals /n")

        while self.KeepAdvancingSolutionLoop():
            self.time = self._AdvanceTime()
            self.InitializeSolutionStep()
            self._GetSolver().Predict()
            is_converged = self._GetSolver().SolveSolutionStep()

            ##########################################################

            # print("/n ::TESTING:: START Calculate Jacobian /n")
            # for element in self._GetSolver().GetComputingModelPart().Elements:
            #     J = element.GetGeometry().Jacobian(0)
            #     Jred = KratosMultiphysics.Matrix([[J[0,0], J[0,1]],[J[1,0], J[1,1]]])
            #     print("J: ", element.Id, J)
            #     #print("Jred: ", element.Id, Jred)
            #     # element.GetNodes()[0].X -= 0.3xx
            #     # element.GetNodes()[0].Y -= 0.1
            #     element.GetNodes()[0].Z = 10
            #     element.GetNodes()[1].Z = 10
            #     element.GetNodes()[2].Z = 10
            #     # for node in element.GetNodes():
            #     #     # node.X = node.X0
            #     #     # node.Y = node.Y0
            #     #     node.Z = node.Z0
            #     J0 = element.GetGeometry().Jacobian(0)
            #     J0red = KratosMultiphysics.Matrix([[J0[0,0], J0[0,1]],[J0[1,0], J0[1,1]]])
            #     J0red_inv = self.InvertMatrix(J0red)
            #     print("J0: ", element.Id, J0)
            #     #print("J0red: ", element.Id, J0red)
            #     #print("J0red_inv: ", element.Id, J0red_inv)
            #     F = Jred * J0red_inv
            #     # print("F: ", element.Id, F)
            # print("/n ::TESTING:: FINISH Calculate Jacobian /n")
            # hehehehe

            DN1DXI = -1.0
            DN2DXI = 1.0
            DN3DXI = 0.0
            
            DNXI = KratosMultiphysics.Matrix(3,9)
            # DN1X1
            DNXI[0,0] = DN1DXI
            DNXI[0,1] = 0.0
            DNXI[0,2] = 0.0
            DNXI[1,0] = 0.0
            DNXI[1,1] = DN1DXI
            DNXI[1,2] = 0.0
            DNXI[2,0] = 0.0
            DNXI[2,1] = 0.0
            DNXI[2,2] = DN1DXI
            # DN2XI
            DNXI[0,3] = DN2DXI
            DNXI[0,4] = 0.0
            DNXI[0,5] = 0.0
            DNXI[1,3] = 0.0
            DNXI[1,4] = DN2DXI
            DNXI[1,5] = 0.0
            DNXI[2,3] = 0.0
            DNXI[2,4] = 0.0
            DNXI[2,5] = DN2DXI
            # DN3XI
            DNXI[0,6] = DN3DXI
            DNXI[0,7] = 0.0
            DNXI[0,8] = 0.0
            DNXI[1,6] = 0.0
            DNXI[1,7] = DN3DXI
            DNXI[1,8] = 0.0
            DNXI[2,6] = 0.0
            DNXI[2,7] = 0.0
            DNXI[2,8] = DN3DXI

            # DNXI[0] = [DN1DXI,      0,       0,   DN2DXI,      0,       0,   DN3DXI,   0,     0]
            # DNXI[1] = [     0,  DN1DXI,      0,        0,  DN2DXI,      0,       0,   DN3DXI, 0]
            # DNXI[2] = [     0,       0, DN1DXI,        0,       0, DN2DXI,       0,   0, DN3DXI]

            DN1DETA = -1.0
            DN2DETA = 0.0
            DN3DETA = 1.0
            
            DNETA = KratosMultiphysics.Matrix(3,9)
            # DN1X1
            DNETA[0,0] = DN1DETA
            DNETA[0,1] = 0.0
            DNETA[0,2] = 0.0
            DNETA[1,0] = 0.0
            DNETA[1,1] = DN1DETA
            DNETA[1,2] = 0.0
            DNETA[2,0] = 0.0
            DNETA[2,1] = 0.0
            DNETA[2,2] = DN1DETA
            # DN2XI
            DNETA[0,3] = DN2DETA
            DNETA[0,4] = 0.0
            DNETA[0,5] = 0.0
            DNETA[1,3] = 0.0
            DNETA[1,4] = DN2DETA
            DNETA[1,5] = 0.0
            DNETA[2,3] = 0.0
            DNETA[2,4] = 0.0
            DNETA[2,5] = DN2DETA
            # DN3XI
            DNETA[0,6] = DN3DETA
            DNETA[0,7] = 0.0
            DNETA[0,8] = 0.0
            DNETA[1,6] = 0.0
            DNETA[1,7] = DN3DETA
            DNETA[1,8] = 0.0
            DNETA[2,6] = 0.0
            DNETA[2,7] = 0.0
            DNETA[2,8] = DN3DETA
            
            # DNETA[0] = [DN1DETA,      0,       0,   DN2DETA,      0,       0,   DN3DETA,   0,     0]
            # DNETA[1] = [     0,  DN1DETA,      0,        0,  DN2DETA,      0,       0,   DN3DETA, 0]
            # DNETA[2] = [     0,       0, DN1DETA,        0,       0, DN2DETA,       0,   0, DN3DETA]

            print("DNXI: ", DNXI)
            print("DNETA: ", DNETA)

            for element in self._GetSolver().GetComputingModelPart().Elements:
                # Jacob = KratosMultiphysics.Matrix([DNXI * element.GetNodes().GetCoordinates()], [DNETA * element.GetNodes().GetCoordinates()])
                coordinates = [element.GetNodes()[0].X, element.GetNodes()[0].Y, element.GetNodes()[0].Z,
                               element.GetNodes()[1].X, element.GetNodes()[1].Y, element.GetNodes()[1].Z,
                               element.GetNodes()[2].X, element.GetNodes()[2].Y, element.GetNodes()[2].Z,]
                # print("ELEMENT:", coordinates)

                Jac_xi = DNXI * coordinates
                Jac_eta = DNETA * coordinates
                J_total = KratosMultiphysics.Matrix([Jac_xi,Jac_eta])
                J_total = J_total.transpose()
                # print(element.Id, "J_xi: ", Jac_xi)
                # print(element.Id, "J_eta: ", Jac_eta)
                # print(element.Id, "J_tot: ", J_total)
                J_Kratos = element.GetGeometry().Jacobian(0)
                # print(element.Id, "J_Kratos: ", J_Kratos)
                J_diff = J_total - J_Kratos
                print(element.Id, "J_diff = ", J_diff)
            


            ##########################################################
            
            self.__CheckIfSolveSolutionStepReturnsAValue(is_converged)
            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

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
        for elm in self._GetSolver().GetComputingModelPart().Elements:
            for node in elm.GetNodes():
                node.X0 = node.X
                node.Y0 = node.Y
                node.Z0 = 0.0
                print(node.Id, node.X0, node.Y0, node.Z0) 
        

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
