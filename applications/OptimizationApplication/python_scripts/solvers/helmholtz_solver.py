import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ConvertContainers

class HelmholtzSolver(PythonSolver):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "model_part_name"          : "",
            "reference_model_part_name": "",
            "echo_level"               : 0,
            "linear_solver_settings"   : {}
        }""")

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters):
        super().__init__(model, settings)

        model_part_name = self.settings["model_part_name"].GetString()
        reference_model_part_name = self.settings["reference_model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        self.reference_model_part = self.model[reference_model_part_name]
        self.main_model_part = self.model.CreateModelPart(model_part_name)

        self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = self.reference_model_part[Kratos.DOMAIN_SIZE]

    def AddVariables(self):
        # we add variables to the reference model part
        # then these will be copied to the helmholtz model part
        # with the connectivity preserve modeller
        self.reference_model_part.AddNodalSolutionStepVariable(KratosOA.HELMHOLTZ_SOURCE_DENSITY)
        self.reference_model_part.AddNodalSolutionStepVariable(KratosOA.HELMHOLTZ_VAR_DENSITY)

    def AddDofs(self):
        # we add dofs to the reference model part
        # then these will be copied to the helmholtz model part
        # with the connectivity preserve modeller
        Kratos.VariableUtils().AddDof(KratosOA.HELMHOLTZ_VAR_DENSITY, self.reference_model_part)

    def GetDofVariable(self):
        return KratosOA.HELMHOLTZ_VAR_DENSITY

    def GetElementName(self) -> str:
        return "HelmholtzBulkTopology3D4N"

    def GetConditionName(self) -> str:
        return ""

    def GetDofsList(self):
        """This function creates and returns a list with the DOFs defined in the conditions and elements specifications
        Note that this requires the main_model_part to be already set, that is to say to have already performed the element substitution (see PrepareModelPart).
        """
        return Kratos.SpecificationsUtilities.GetDofsListFromSpecifications(self.GetComputingModelPart())

    def ImportModelPart(self):
        # here we do nothing, since we are not importing any model part
        pass

    def PrepareModelPart(self):
        # generate a new model part with same nodes, but different elements and
        # conditions using same geometries for Helmholtz solving
        connectivity_preserve_modeler = Kratos.ConnectivityPreserveModeler()

        # set the variables list
        KratosOA.OptimizationUtils.CopySolutionStepVariablesList(self.GetComputingModelPart(), self.reference_model_part)

        # generate the duplicate model part
        if self.GetConditionName() != "":
            connectivity_preserve_modeler.GenerateModelPart(
                self.reference_model_part,
                self.GetComputingModelPart(),
                self.GetElementName(),
                self.GetConditionName())
        else:
            connectivity_preserve_modeler.GenerateModelPart(
                self.reference_model_part,
                self.GetComputingModelPart(),
                self.GetElementName())

        self.GetComputingModelPart().ProcessInfo[Kratos.DOMAIN_SIZE] = self.reference_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

        # initialize all element properties with dummy filter radius
        for properties in self.GetComputingModelPart().GetProperties():
            properties[KratosOA.HELMHOLTZ_RADIUS_DENSITY] = 0.1

        Kratos.Logger.PrintInfo(self.__class__.__name__, f"Populated {self.GetComputingModelPart().FullName()} with {self.GetElementName()} and {self.GetConditionName()} according to reference {self.reference_model_part.FullName()}")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        Kratos.ModelPartIO(name_out_file, Kratos.IO.WRITE).WriteModelPart(self.GetComputingModelPart())
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Model export finished.")

    def GetMinimumBufferSize(self):
        return 1

    def Initialize(self):
        # compute the number of neighbour nodes
        dummy_container = KratosOA.ElementContainerVariableDataHolder(self.reference_model_part)
        self.number_of_neighbour_nodes = KratosOA.NodalContainerVariableDataHolder(self.reference_model_part)
        KratosOA.ContainerVariableDataHolderUtils.ComputeNumberOfNeighbourEntities(self.number_of_neighbour_nodes, dummy_container)
        self.number_of_neighbour_nodes.AssignDataToContainerVariable(Kratos.NUMBER_OF_NEIGHBOUR_ELEMENTS)
        self.GetComputingModelPart().ProcessInfo[KratosOA.COMPUTE_CONTROL_DENSITIES] = False

    def AdvanceInTime(self, current_time):
        self.GetComputingModelPart().CloneTimeStep(self.reference_model_part.ProcessInfo[Kratos.TIME] + 1)
        return self.reference_model_part.ProcessInfo[Kratos.TIME]

    def InitializeSolutionStep(self):
        # free the blocked flag
        Kratos.VariableUtils().SetFlag(Kratos.BLOCKED, False, self.GetComputingModelPart().Nodes)

        self._GetSolutionStrategy().InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        # set all the non-fixed dofs to zero
        Kratos.VariableUtils().SetVariable(KratosOA.HELMHOLTZ_VAR_DENSITY, 0.0, self.reference_model_part.Nodes, Kratos.BLOCKED, False)

        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        if not is_converged:
            msg  = "Helmholtz solver did not converge for step " + str(self.reference_model_part.ProcessInfo[Kratos.STEP]) + "\n"
            Kratos.Logger.PrintWarning(self.__class__.__name__, msg)
        return is_converged

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()

    def Check(self):
        self._GetSolutionStrategy().Check()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def GetComputingModelPart(self):
        return self.main_model_part

    def _GetScheme(self):
        if not hasattr(self, '_scheme'):
            self._scheme = self._CreateScheme()
        return self._scheme

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    def _CreateScheme(self):
        return Kratos.ResidualBasedIncrementalUpdateStaticScheme()

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        return linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        return Kratos.ResidualBasedBlockBuilderAndSolver(linear_solver)

    def _CreateSolutionStrategy(self):
        model_part = self.GetComputingModelPart()
        scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        compute_reactions = False
        reform_dofs = False
        calculate_dx_norm = False
        move_mesh = False
        strategy = Kratos.ResidualBasedLinearStrategy(
            model_part,
            scheme,
            builder_and_solver,
            compute_reactions,
            reform_dofs,
            calculate_dx_norm,
            move_mesh)

        strategy.SetEchoLevel(self.echo_level)
        return strategy

    def __CheckInputContainer(self, input_container_variable_data_holder: ContainerVariableDataHolderUnion):
        if self.reference_model_part != input_container_variable_data_holder.GetModelPart().GetRootModelPart():
            raise RuntimeError(f"Root model part mismatch. [ Reference model part = {self.reference_model_part.FullName()}, input model part = {input_container_variable_data_holder.GetModelPart().FullName()}]")

    def SetData(self, input_container_variable_data_holder: ContainerVariableDataHolderUnion):
        self.__CheckInputContainer(input_container_variable_data_holder)

        # free the blocked flag
        input_model_part: Kratos.ModelPart = input_container_variable_data_holder.GetModelPart()
        Kratos.VariableUtils().SetFlag(Kratos.BLOCKED, False, input_model_part.Nodes)

        historical_container_variable_data_holder = KratosOA.HistoricalContainerVariableDataHolder(input_model_part)
        self.ConvertSolverContainer(input_container_variable_data_holder, historical_container_variable_data_holder)

        # apply the variable values
        historical_container_variable_data_holder.AssignDataToContainerVariable(KratosOA.HELMHOLTZ_SOURCE_DENSITY)

    def FixDofs(self, input_container_variable_data_holder: ContainerVariableDataHolderUnion):
        self.__CheckInputContainer(input_container_variable_data_holder)

        # fix the blocked flag
        input_model_part: Kratos.ModelPart = input_container_variable_data_holder.GetModelPart()
        Kratos.VariableUtils().SetFlag(Kratos.BLOCKED, True, input_model_part.Nodes)
        Kratos.VariableUtils().ApplyFixity(KratosOA.HELMHOLTZ_VAR_DENSITY, True, input_model_part.Nodes)

        # get the historical contaier values and set it
        historical_container_variable_data_holder = KratosOA.HistoricalContainerVariableDataHolder(input_model_part)
        self.ConvertSolverContainer(input_container_variable_data_holder, historical_container_variable_data_holder)
        historical_container_variable_data_holder.AssignDataToContainerVariable(KratosOA.HELMHOLTZ_VAR_DENSITY)

    def GetData(self, output_container_variable_data_holder: ContainerVariableDataHolderUnion):
        self.__CheckInputContainer(output_container_variable_data_holder)

        histoical_container_variable_data_holder = KratosOA.HistoricalContainerVariableDataHolder(output_container_variable_data_holder.GetModelPart())
        histoical_container_variable_data_holder.ReadDataFromContainerVariable(KratosOA.HELMHOLTZ_VAR_DENSITY)

        self.ConvertSolverContainer(histoical_container_variable_data_holder, output_container_variable_data_holder)

    def ConvertSolverContainer(self, input_data_container: ContainerVariableDataHolderUnion, output_data_container: ContainerVariableDataHolderUnion):
        ConvertContainers(input_data_container, output_data_container, self.number_of_neighbour_nodes)
