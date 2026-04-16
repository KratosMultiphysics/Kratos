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
        self.fixity_ids = project_parameters["process_combinations"]["fixity_ids"]
        self.model_part_name = self.solver_settings["model_part_name"].GetString()
        #self.combination_by_constraint = self.__GroupCombinationsByConstraints()
        self.combination_solutions = {}
        self.primitive_solutions = {}

        #TODO: Change the input format so all of these are accessible via their ids?
        if self.project_parameters.Has("process_combinations"):
            if self.project_parameters["process_combinations"].Has("fixity_processes"):
                self.fixity_processes = self.__CreateFixityDictionary()
            if self.project_parameters["process_combinations"].Has("mpc_processes"):
                #TODO: Create a method "self.__CreateMpcDictionary()"
                pass
            if self.project_parameters["process_combinations"].Has("combinations"):
                self.combinations = self.__CreateCombinationDictionary()
                #print(self.combinations)

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

    def Run(self):

        self.__StoreReferenceState()
        self.__SolveConstraintState()
        self.__OutputSolutionStep()
        for comb, dx in self.combination_solutions.items():
            print(f"Combination {comb} solution vector dx: {dx}")

    def __SolveConstraintState(self):

        load_ids = self.__GetPrimaryLoads()
        self.__ApplyFixities()
        self.__SolvePrimaryLoads(load_ids)
        self.__SolveCombinations()
        KratosMultiphysics.Logger.PrintInfo("::[MultiLoadConstraintAnalysis]::", "Finished constraint-state solve")

    def __OutputSolutionStep(self):

        for combination_id, solution in self.combination_solutions.items():
            self.dofset.SetValues(solution)
            self.__PrintCombinationOutput(combination_id)
            self.__RestoreReferenceState()

    def __PrintCombinationOutput(self, combination_id):

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


    def __GetPrimaryLoads(self):
        load_ids = set()

        for combination in self.combination_list.values():
            for load in combination["loads"].values():
                load_ids.add(load["load_process_id"].GetString())

        return list(load_ids)
    
    def __SolveCombinations(self):
        KratosMultiphysics.Logger.PrintInfo("::[MultiLoadConstraintAnalysis]::", "Solving combinations")
        for combination in self.combination_list.values():
            combination_factor = combination["combination_factor"].GetDouble()
            combination_id = combination["combination_id"].GetInt()
            combination_Dx = None
            for load in combination["loads"].values():
                load_factor = load["load_factor"].GetDouble()
                load_id = load["load_process_id"].GetString()
                dx = self.primitive_solutions[load_id]
                if combination_Dx == None:
                    combination_Dx = dx.copy()
                    combination_Dx.SetValue(0.0)
                combination_Dx += load_factor * dx
            combination_Dx *= combination_factor
            self.combination_solutions[combination_id] = combination_Dx.copy()

    
    def __SolvePrimaryLoads(self, load_ids):
        KratosMultiphysics.Logger.PrintInfo("::[MultiLoadConstraintAnalysis]::", f"Primitive loads to solve: {load_ids}")
        for id in load_ids:
            self.__RestoreReferenceState()
            rhs = self.__GetRHS(id)
            lhs = self.__GetRawLHS()
            #linear_system = self.init_strategy_data.GetLinearSystem()
            #linear_system.SetVector(rhs, KratosMultiphysics.Future.DenseVectorTag.RHS)
            dx = KratosMultiphysics.SystemVector(rhs.Size())
            linear_system = self.__CreateLinearSystem(lhs, rhs, dx)
            strategy_data = self.__CreateStrategyData()
            scheme = self.__InitializeScheme(strategy_data)
            self.__InitializeStrategyData(strategy_data, linear_system)
            self.__BuildEffectiveSystem(scheme, strategy_data)
            self.__SolveLinearSystem(strategy_data)
            scheme.Update(strategy_data)
            self.__StorePrimarySolution(id, strategy_data)

    def __StorePrimarySolution(self, id, strategy_data):
        dx = strategy_data.GetLinearSystem().GetVector(KratosMultiphysics.Future.DenseVectorTag.Dx)
        self.primitive_solutions[id] = dx.copy()

    #def __OutputSolutionStep(self):
    #    for combination_id, solution in self.combination_solutions.items():
    #        current_combination = self.__GetCombination(combination_id)
    #        if self.project_parameters.Has("output_processes"):
    #            if self.project_parameters["output_processes"].Has("vtk_output"):
    #                if current_combination.Has("combination_label"):
    #                    label = current_combination["combination_label"].GetString()
    #                    self.project_parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f"vtk_output/{label}")
    #                else:    
    #                    self.project_parameters["output_processes"]["vtk_output"][0]["Parameters"]["output_path"].SetString(f"vtk_output/combination_{combination_id}")
    #                order_processes_initialization = self._GetOrderOfOutputProcessesInitialization()
    #                self._list_of_output_processes = self._CreateProcesses("output_processes", order_processes_initialization)
    #                for process in self._list_of_output_processes:
    #                    process.ExecuteBeforeOutputStep()
    #                    process.PrintOutput()

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
        
    def __ApplyFixities(self):
        KratosMultiphysics.Logger.PrintInfo("::[MultiLoadConstraintAnalysis]::", "Applying fixities")
        for id in self.fixity_ids.values():
            id = id.GetString()
            definition = self.__GetFixity(id)
            process = self.__CreateProcess(definition)
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()
            process.ExecuteInitializeSolutionStep()

    def __BuildEffectiveSystem(self, scheme, strategy_data):

        #scheme.BuildMasterSlaveConstraints(strategy_data)
        scheme.BuildLinearSystemConstraints(strategy_data)
        scheme.ApplyLinearSystemConstraints(strategy_data)

    def __InitializeStrategyData(self, strategy_data, linear_system):

        strategy_data.SetLinearSystem(linear_system)
        strategy_data.SetDofSet(self.dofset)

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

        self.init_strategy_data = output_data[stage_name]["strategy_data"]
        self.lhs = self.init_strategy_data.GetLinearSystem().GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS)
        self.rhss = output_data[stage_name]["rhss"]
        self.fixity_processes = output_data[stage_name]["fixities"]

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
        self.dofset = self.init_strategy_data.GetDofSet()
        self.reference_state = self.dofset.GetValues()

    def __RestoreReferenceState(self):
        self.dofset.SetValues(self.reference_state)
