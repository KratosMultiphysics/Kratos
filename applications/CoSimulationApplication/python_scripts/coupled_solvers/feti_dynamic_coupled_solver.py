# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
import KratosMultiphysics.CoSimulationApplication as CoSim

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

# Other imports
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

def Create(settings, models, solver_name):
    return FetiDynamicCoupledSolver(settings, models, solver_name)

class FetiDynamicCoupledSolver(CoSimulationCoupledSolver):
    def __init__(self, settings, models, solver_name):
        super().__init__(settings, models, solver_name)

        if len(self.solver_wrappers.items()) != 2:
            raise Exception("FETI solver only works with two solvers!")

        # Add solution step variables while models are empty
        self.__CheckSolversCompatibility()
        self.__AddNodalSolutionStepVariables()
        self.is_initial_step = True

    def Initialize(self):
        super().Initialize()

    def AdvanceInTime(self, current_time):
        # The CoSimulation runs at the SMALLEST Timestep.
        # The solver with the smallest timestep will dictate the CoSimulation.
        # The solver(s) with the larger timesteps will be called only at times that match their time
        advanced_time = current_time
        self.departing_time = current_time # the time we are departing from - this is used to sync everything

        if self.is_initial_step:
            self._solver_delta_times = {}
            advanced_time = 1E20
            for solver_name, solver in self.solver_wrappers.items():
                self._solver_delta_times[solver_name] = solver.AdvanceInTime(current_time)
                advanced_time = min(advanced_time,self._solver_delta_times[solver_name])

            self.is_initial_step = False

            # Initialize the FETI utilites
            self.__InitializeFetiMethod()
        else:
            for solver_name, solver in self.solver_wrappers.items():
                if self.SolverSolvesAtThisTime(solver_name):
                    solver_time = solver.AdvanceInTime(self.departing_time)
                    if self._solver_origin_dest_dict[solver_name] == CoSim.FetiSolverIndexType.Destination:
                        advanced_time = solver_time #only advance global time finely

        if advanced_time == current_time:
            raise Exception("No solvers advanced any timestep.")
        else:
            return advanced_time

    def SolverSolvesAtThisTime(self, solver_name):
        solver_delta_time = self._solver_delta_times[solver_name]
        # the following only works if timesteps are multiple of each other
        time_error = (self.departing_time % solver_delta_time)
        if time_error < 1E-12 or abs(time_error - solver_delta_time) < 1E-12:
            return True
        else:
            return False

    def InitializeSolutionStep(self):
        for solver_name, solver in self.solver_wrappers.items():
            if self.SolverSolvesAtThisTime(solver_name):
                solver.InitializeSolutionStep()

        for coupling_op_name, coupling_op in self.coupling_operations_dict.items():
            if self.__CouplingOpActsNow(coupling_op_name):
                coupling_op.InitializeSolutionStep()

    def Predict(self):
        for solver_name, solver in self.solver_wrappers.items():
            if self.SolverSolvesAtThisTime(solver_name):
                solver.Predict()

    def FinalizeSolutionStep(self):
        for solver_name, solver in self.solver_wrappers.items():
            if self.SolverSolvesAtThisTime(solver_name):
                solver.FinalizeSolutionStep()

        for coupling_op_name, coupling_op in self.coupling_operations_dict.items():
            if self.__CouplingOpActsNow(coupling_op_name):
                coupling_op.FinalizeSolutionStep()

    def OutputSolutionStep(self):
        for solver_name, solver in self.solver_wrappers.items():
            if self.SolverSolvesAtThisTime(solver_name):
                solver.OutputSolutionStep()

    def SolveSolutionStep(self):
        for coupling_op_name, coupling_op in self.coupling_operations_dict.items():
            if self.__CouplingOpActsNow(coupling_op_name):
                coupling_op.InitializeCouplingIteration()

        for solver_name, solver in self.solver_wrappers.items():
            if self.SolverSolvesAtThisTime(solver_name):
                #self._SynchronizeInputData(solver_name) @phil not needed since corrections are applied within the feti cpp
                solver.SolveSolutionStep()
                #self._SynchronizeOutputData(solver_name) @phil not needed since corrections are applied within the feti cpp
                self.__SendStiffnessMatrixToUtility(solver_name)

        self.feti_coupling.EquilibrateDomains()

        for coupling_op_name, coupling_op in self.coupling_operations_dict.items():
            if self.__CouplingOpActsNow(coupling_op_name):
                coupling_op.FinalizeCouplingIteration()

        return True

    def __AddNodalSolutionStepVariables(self):
        for solver_name, solver in self.solver_wrappers.items():
            structure = solver.model.GetModelPart("Structure")
            structure.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)

    def __InitializeFetiMethod(self):
        # Create vector of solver indices for convenience
        self.__CreateSolverOriginDestDict()

        # Check timestep ratio is valid and add to settings
        timestep_ratio = self._CalculateAndCheckTimestepRatio()
        self.settings.AddInt('timestep_ratio',int(timestep_ratio))

        # get mapper parameters
        self.mapper_parameters = self.data_transfer_operators_dict["mapper"].settings["mapper_settings"]
        mapper_type = self.mapper_parameters["mapper_type"].GetString()

        # get mapper origin and destination modelparts
        origin_modelpart_name = self.mapper_parameters["modeler_parameters"]["origin_interface_sub_model_part_name"].GetString()
        destination_modelpart_name = self.mapper_parameters["modeler_parameters"]["destination_interface_sub_model_part_name"].GetString()
        for solver_name, solver in self.solver_wrappers.items():
            if self._solver_origin_dest_dict[solver_name] == CoSim.FetiSolverIndexType.Origin:
                self.model_part_origin_interface = self.solver_wrappers[solver_name].model.GetModelPart(origin_modelpart_name)
            else:
                self.model_part_destination_interface = self.solver_wrappers[solver_name].model.GetModelPart(destination_modelpart_name)

        # manually create mapper
        mapper_create_fct = KratosMapping.MapperFactory.CreateMapper
        self.mapper = mapper_create_fct(self.model_part_origin_interface, self.model_part_destination_interface, self.mapper_parameters.Clone())

        # get interface modelparts created by the mapper modeler
        self.modelpart_interface_origin_from_mapper = self.mapper.GetInterfaceModelPartOrigin()
        self.modelpart_interface_destination_from_mapper = self.mapper.GetInterfaceModelPartDestination()

        # Create feti class instance
        self.feti_coupling = CoSim.FetiDynamicCouplingUtilities(
            self.modelpart_interface_origin_from_mapper,
            self.modelpart_interface_destination_from_mapper,
            self.settings)

        # set the mapper
        if mapper_type == "coupling_geometry":
            self.feti_coupling.SetMappingMatrix(self.mapper.GetMappingMatrix())
        else:
            raise Exception("Dynamic coupled solver currently only compatible with the coupling_geometry mapper.")

        # The origin and destination interfaces from the mapper submitted above are both
        # stored within the origin modelpart. Now we submit the 'original' origin and destination
        # interface model parts stored on the origin and destination models to get access to the
        # origin and destination models.
        self.feti_coupling.SetOriginAndDestinationDomainsWithInterfaceModelParts(
            self.model_part_origin_interface,
            self.model_part_destination_interface)

        # create the solver
        linear_solver = self._CreateLinearSolver()
        self.feti_coupling.SetLinearSolver(linear_solver)

        # Set origin initial velocities
        self.feti_coupling.SetOriginInitialKinematics()

        # Create output-solver relation dict to ensure mixed timestep ouput is handled properly
        self.__CreateOutputSolverDict()

    def __CreateSolverOriginDestDict(self):
        self._solver_origin_dest_dict = {}
        for solver_index in range(self.settings["coupling_sequence"].size()):
            ordered_solver_name = self.settings["coupling_sequence"][solver_index]["name"].GetString()
            if solver_index == 0:
                self._solver_origin_dest_dict[ordered_solver_name] = CoSim.FetiSolverIndexType.Origin
            else:
                self._solver_origin_dest_dict[ordered_solver_name] = CoSim.FetiSolverIndexType.Destination

    def __SendStiffnessMatrixToUtility(self, solver_name):
        beta = 0.0
        solver_index = self._solver_origin_dest_dict[solver_name]
        if solver_index == CoSim.FetiSolverIndexType.Origin:
            beta = self.settings["origin_newmark_beta"].GetDouble()
        else:
            beta = self.settings["destination_newmark_beta"].GetDouble()

        if abs(beta-0.25) < 1e-9:
                system_matrix = self._GetSolverStrategy(solver_name).GetSystemMatrix()
                self.feti_coupling.SetEffectiveStiffnessMatrixImplicit(system_matrix,solver_index)
        else:
            self.feti_coupling.SetEffectiveStiffnessMatrixExplicit(solver_index)

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
            return linear_solver_factory.ConstructSolver(linear_solver_configuration)
        else:
            KM.Logger.PrintInfo('::[FETISolver]:: No linear solver was specified, using fastest available solver')
            return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    def _CalculateAndCheckTimestepRatio(self):
        # Check timestep ratio is valid
        timesteps = [0.0] * len(self._solver_delta_times)
        for solver_index in range(self.settings["coupling_sequence"].size()):
            ordered_solver_name = self.settings["coupling_sequence"][solver_index]["name"].GetString()
            timesteps[solver_index] = self._solver_delta_times[ordered_solver_name]
        timestep_ratio = timesteps[0] / timesteps[1]
        if timestep_ratio < 0.99 or int(timestep_ratio) % timestep_ratio > 1E-12:
            raise Exception("The timestep ratio between origin and destination domains is invalid. It must be a positive integer greater than 1.")
        return timestep_ratio

    def _GetSolverStrategy(self,solverName):
        # This is a utility method to get access to the solver's strategy, later used to access the system matrix.
        # Provision to expand to other solver wrappers in the future.
        solver_type = str(self.solver_wrappers[solverName]._ClassName())
        if solver_type == "StructuralMechanicsWrapper":
            return self.solver_wrappers[solverName]._analysis_stage._GetSolver().get_mechanical_solution_strategy()
        else:
            raise Exception("_GetSolverStrategy not implemented for solver wrapper = " + solver_type)

    def __CheckSolversCompatibility(self):
        compatible_solver_wrappers = ['StructuralMechanicsWrapper']
        for solver_name, solver in self.solver_wrappers.items():
            solver_type = str(solver._ClassName())
            if solver_type not in compatible_solver_wrappers:
                raise Exception("The coupled solver '" + solver_type + "' is not yet compatible with the FETI coupling")

    def __CreateOutputSolverDict(self):
        # This function links each coupling output entry with a solver (if applicable)
        # This means we can now apply SolverSolvesAtThisTime() to the coupling output
        self.output_solver_dict = {}
        for coupling_op_name in self.coupling_operations_dict.keys():
            coupling_op_type = self.settings["coupling_operations"][coupling_op_name]["type"].GetString()
            if coupling_op_type == "coupling_output":
                coupling_output_solver_name = self.settings["coupling_operations"][coupling_op_name]["solver"].GetString()
                self.output_solver_dict[coupling_op_name] = coupling_output_solver_name

    def __CouplingOpActsNow(self,couplingOpName):
        if couplingOpName in self.output_solver_dict:
            solver_name = self.output_solver_dict[couplingOpName]
            return self.SolverSolvesAtThisTime(solver_name)
        else:
            return True #only restrict coupling operations for outputs

    @classmethod
    def _GetDefaultParameters(cls):
        this_defaults = KM.Parameters("""{
            "origin_newmark_beta" : -1.0,
            "origin_newmark_gamma" : -1.0,
            "destination_newmark_beta" : -1.0,
            "destination_newmark_gamma" : -1.0,
            "equilibrium_variable" : "VELOCITY",
            "is_disable_coupling" : false,
            "is_linear" : false,
            "linear_solver_settings" : {}
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultParameters())

        return this_defaults
