import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.process_factory import KratosProcessFactory
#-------------------------------------------------------------------
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory


class MultiLoadConstraintAnalysis(AnalysisStage):
    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)
        self.model = model
        self.scheme_settings = project_parameters["future_scheme_settings"]
        self.solver_settings = project_parameters["solver_settings"]
        self.combination_list = project_parameters["process_combinations"]["combinations"]
        self.model_part_name = self.solver_settings["model_part_name"].GetString()
        self.combination_by_constraint = self.__GroupCombinationsByConstraints()
        self.combination_solutions = {}

        #TODO: Change the input format so all of these are accessible via their ids?
        if self.project_parameters.Has("process_combinations"):
            if self.project_parameters["process_combinations"].Has("fixity_processes"):
                self.fixity_processes = self.__CreateFixityDictionary()
            if self.project_parameters["process_combinations"].Has("mpc_processes"):
                #TODO: Create a method "self.__CreateMpcDictionary()"
                pass
            if self.project_parameters["process_combinations"].Has("combinations"):
                self.combinations = self.__CreateCombinationDictionary()
                print(self.combinations)

        if self.model_part_name == "":
            raise Exception('Please specify a model_part name!')
        
        #check if model already has model part, otherwise create it and set domain size
        if self.model.HasModelPart(self.model_part_name):
            self.main_model_part = self.model[self.model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(self.model_part_name)
            domain_size = self.solver_settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def Initialize(self):
        super().Initialize()

    def __GroupCombinationsByConstraints(self):
        combinations_by_constraint = {}

        for combination in self.combination_list.values():
            combination_id = combination["combination_id"].GetInt()

            fixity_ids = []
            if combination.Has("fixity_ids"):
                fixity_ids = [item.GetString() for item in combination["fixity_ids"].values()]

            mpc_ids = []
            if combination.Has("mpc_ids"):
                mpc_ids = [item.GetString() for item in combination["mpc_ids"].values()]

            constraint_key = (
                tuple(sorted(fixity_ids)),
                tuple(sorted(mpc_ids))
            )

            #create a new dictionary key with an empty list
            if constraint_key not in combinations_by_constraint:
                combinations_by_constraint[constraint_key] = []
            #store the current combination id to the key
            combinations_by_constraint[constraint_key].append(combination_id)

        return combinations_by_constraint

    def Run(self):
        #self.debug_check()
        self.__StoreReferenceState()
        self.__SolveCombinations()
        #self.__ResetFixities()
        for comb, dx in self.combination_solutions.items():
            print(f"Combination {comb} solution vector dx: {dx}")

    def __SolveCombinations(self):
        for fixity_group, combination_ids in self.combination_by_constraint.items():
            self.__ResetFixities()
            self.__ApplyFixities(fixity_group)
            for combination_id in combination_ids:
                self.__SolveCombination(combination_id)
                self.__OutputSolutionStep(combination_id)
                self.__RestoreReferenceState()
        print(self.combination_solutions)

    def __SolveCombination(self, combination_id):
        combination = self.__GetCombination(combination_id)
        combination_rhs = self.__CreateCombinationRHS(combination)
        lhs = self.__GetRawLHS()
        dx = KratosMultiphysics.SystemVector(combination_rhs.Size())
        linear_system = self.__CreateLinearSystem(lhs, combination_rhs, dx)
        strategy_data = self.__CreateStrategyData()
        scheme = self.__InitializeScheme(strategy_data)
        self.__InitializeStrategyData(strategy_data, linear_system)
        self.__BuildEffectiveSystem(scheme, strategy_data)
        self.__SolveLinearSystem(strategy_data)
        scheme.Update(strategy_data)
        self.__StoreCombinationSolution(combination_id, strategy_data)

    def __OutputSolutionStep(self, combination_id):
        current_combination = self.__GetCombination(combination_id)
        if self.project_parameters.Has("output_processes"):
            if self.project_parameters["output_processes"].Has("vtk_output"):
                if current_combination.Has("combination_label"):
                    label = current_combination["combination_label"].GetString()
                    self.project_parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f"vtk_output/{label}")
                else:    
                    self.project_parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f"vtk_output/combination_{combination_id}")
                order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
                self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
                for process in self._list_of_output_processes:
                    process.ExecuteBeforeOutputStep()
                    process.PrintOutput()

    def GetFinalData(self):
        return self.combination_solutions

    def __StoreCombinationSolution(self, combination_id, strategy_data):
        dx = strategy_data.GetLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx)
        self.combination_solutions[combination_id] = dx.copy()

    def __SolveLinearSystem(self, strategy_data):
        effective_system = strategy_data.GetEffectiveLinearSystem()
        solver = self.__GetLinearSolver()
        solver.Initialize(effective_system)
        solver.InitializeSolutionStep(effective_system)
        solver.PerformSolutionStep(effective_system)
        solver.FinalizeSolutionStep(effective_system)
        solver.Clear()

    def __GetLinearSolver(self):
        return KratosMultiphysics.Future.SkylineLUFactorizationSolver()
        
    def __ApplyFixities(self, fixity_tuple):
        fixities, mpcs = fixity_tuple
        for id in fixities:
            definition = self.__GetFixity(id)
            process = self.__CreateProcess(definition)
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()
            process.ExecuteInitializeSolutionStep()
    
    def __CreateCombinationRHS(self, combination):
        """Build the RHS for one combination by applying the load factors"""
        loads = combination["loads"]
        combination_scale_factor = combination["combination_factor"].GetDouble()
        combination_rhs = None
        for load in loads.values():
            load_id = load["load_process_id"].GetString()
            load_factor = load["load_factor"].GetDouble()
            rhs = self.__GetRHS(load_id)
            if combination_rhs is None:
                combination_rhs = rhs.copy()
                combination_rhs.SetValue(0.0)

            combination_rhs += load_factor * rhs

        combination_rhs *= combination_scale_factor
        combination_id = combination["combination_id"].GetInt()
        print(f"COMBINATION: {combination_id} with rh: {combination_rhs}")
        return combination_rhs
    
    def __InitializeSchemeDofs(self, scheme, dof_set, effective_dof_set):

        scheme.SetUpDofArrays(dof_set, effective_dof_set)
        scheme.SetUpSystemIds(dof_set, effective_dof_set)

    def __BuildEffectiveSystem(self, scheme, strategy_data):

        #scheme.BuildMasterSlaveConstraints(strategy_data)
        scheme.BuildLinearSystemConstraints(strategy_data)
        scheme.ApplyLinearSystemConstraints(strategy_data)

    def __InitializeStrategyData(self, strategy_data, linear_system):

        strategy_data.SetLinearSystem(linear_system)
        #strategy_data.SetDofSet(dof_set)
        #strategy_data.SetEffectiveDofSet(effective_dof_set)

    def __CreateStrategyData(self):

        return KratosMultiphysics.Future.ImplicitStrategyData() #stores all arrays required to setup linear system

    def __InitializeScheme(self, strategy_data):

        scheme = KratosMultiphysics.Future.StaticScheme(self.main_model_part, self.scheme_settings) # scheme settings from json file
        scheme.Initialize(strategy_data) #initialize the scheme (will reset lhs, rhs and dx...)
        return scheme


    def __CreateProcess(self, process_definition):

        factory = KratosProcessFactory(self.model) #instantiate factory
        # the function "ConstructListOfProcesses" expects a Parameters array, so the process_list is stores in one here
        process_list = KratosMultiphysics.Parameters("[]") 
        process_list.Append(process_definition) 
        return factory.ConstructListOfProcesses(process_list)[0] # function returns a list, but we create processes one by one so it only holds one entry
    
    def GetProjectOutputData(self, stage_name, output_data):

        self.lhs = output_data[stage_name]["lhs"]
        self.rhss = output_data[stage_name]["rhss"]

    def debug_check(self):

        print(self.combination_by_constraint)
        for node in self.main_model_part.Nodes:
            print(f"Node {node.Id}")
            print("  UX fixed:", node.GetDof(KratosMultiphysics.DISPLACEMENT_X).IsFixed())
            print("  UY fixed:", node.GetDof(KratosMultiphysics.DISPLACEMENT_Y).IsFixed())
            print("  UZ fixed:", node.GetDof(KratosMultiphysics.DISPLACEMENT_Z).IsFixed())
            print("  RX fixed:", node.GetDof(KratosMultiphysics.ROTATION_X).IsFixed())
            print("  RY fixed:", node.GetDof(KratosMultiphysics.ROTATION_Y).IsFixed())
            print("  RZ fixed:", node.GetDof(KratosMultiphysics.ROTATION_Z).IsFixed())

        for id, rhs in self.rhss.items():
            print(rhs)

    def __CreateFixityDictionary(self):

        fixities = {}
        for fixity_definition in self.project_parameters["process_combinations"]["fixity_processes"].values():
            process_id = fixity_definition["process_id"].GetString()
            fixities[process_id] = fixity_definition

        return fixities
    
    def __CreateMpcDictionary(self):
        #TODO: Implement this
        pass
    
    def __CreateCombinationDictionary(self):

        combinations = {}
        for combination in self.project_parameters["process_combinations"]["combinations"].values():
            combination_id = combination["combination_id"].GetInt()
            combinations[combination_id] = combination

        return combinations
    
    def __ResetFixities(self):

        dof_variables = [
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
            KratosMultiphysics.ROTATION_X,
            KratosMultiphysics.ROTATION_Y,
            KratosMultiphysics.ROTATION_Z
        ]

        for node in self.main_model_part.Nodes:
            for dof_variable in dof_variables:
                if node.HasDofFor(dof_variable):
                    node.Free(dof_variable)
    
    def __GetRawLHS(self):
        return self.lhs.copy()
    
    def __GetRHS(self, id):
        return self.rhss[id].copy()
    
    def __GetFixity(self, id):
        return self.fixity_processes[id]
    
    def __GetCombination(self, id):
        return self.combinations[id]
    
    def __CreateLinearSystem(self, lhs, rhs, dx, name="linear_system"):
        linear_system = KratosMultiphysics.Future.LinearSystem(name)
        linear_system.SetMatrix(lhs, KratosMultiphysics.Future.SparseMatrixTag.LHS)
        linear_system.SetVector(rhs, KratosMultiphysics.Future.DenseVectorTag.RHS)
        linear_system.SetVector(dx, KratosMultiphysics.Future.DenseVectorTag.Dx)

        return linear_system
    
    def __StoreReferenceState(self):
        self.reference_state = {}

        variables = [
            KratosMultiphysics.DISPLACEMENT_X,
            KratosMultiphysics.DISPLACEMENT_Y,
            KratosMultiphysics.DISPLACEMENT_Z,
            KratosMultiphysics.ROTATION_X,
            KratosMultiphysics.ROTATION_Y,
            KratosMultiphysics.ROTATION_Z
        ]

        for node in self.main_model_part.Nodes:
            self.reference_state[node.Id] = {}
            for variable in variables:
                if node.SolutionStepsDataHas(variable):
                    self.reference_state[node.Id][variable.Name()] = node.GetSolutionStepValue(variable)

    def __RestoreReferenceState(self):
        #self.dofset.SetValues(self.reference_state)
        variable_map = {
            "DISPLACEMENT_X": KratosMultiphysics.DISPLACEMENT_X,
            "DISPLACEMENT_Y": KratosMultiphysics.DISPLACEMENT_Y,
            "DISPLACEMENT_Z": KratosMultiphysics.DISPLACEMENT_Z,
            "ROTATION_X": KratosMultiphysics.ROTATION_X,
            "ROTATION_Y": KratosMultiphysics.ROTATION_Y,
            "ROTATION_Z": KratosMultiphysics.ROTATION_Z
        }

        for node in self.main_model_part.Nodes:
            for variable_name, value in self.reference_state[node.Id].items():
                node.SetSolutionStepValue(variable_map[variable_name], value)
