# co simulation imports
from ... base_classes.co_simulation_base_solver import CoSimulationBaseSolver
from ... import co_simulation_tools as cs_tools
data_structure = cs_tools.cs_data_structure
import random

def Create(name, cosim_solver_settings):
    return PipeFlowCoSimulationSolver(name, cosim_solver_settings)


class PipeFlowCoSimulationSolver(CoSimulationBaseSolver):
    pipe_flow_solver_count = 0
    def __init__(self, name, custom_settings):
        super(PipeFlowCoSimulationSolver, self).__init__(name, custom_settings)
        self.name = name
        self.pipe_flow_model_part = self.model.CreateModelPart('pipe_flow_model_part')
        for variable_name in self.data_map.keys():
            if(not data_structure.KratosGlobals.HasVariable(variable_name)):
                self.variable_obj = data_structure.Array1DVariable3(variable_name)
                self.pipe_flow_model_part.AddNodalSolutionStepVariable(self.variable_obj)
            else:
                self.variable_obj = data_structure.KratosGlobals.GetVariable(variable_name)
                self.pipe_flow_model_part.AddNodalSolutionStepVariable(self.variable_obj)

        self._SetUpModelpart()
        PipeFlowCoSimulationSolver.pipe_flow_solver_count = PipeFlowCoSimulationSolver.pipe_flow_solver_count + 1

    def Initialize(self):
        self.InitializeIO()
        ## Setting initial value for the data
        data_value = [PipeFlowCoSimulationSolver.pipe_flow_solver_count, PipeFlowCoSimulationSolver.pipe_flow_solver_count, PipeFlowCoSimulationSolver.pipe_flow_solver_count]
        for node in self.pipe_flow_model_part.Nodes:
            for data_name in self.data_map.keys():
                data_obj = data_structure.KratosGlobals.GetVariable(data_name)
                node.SetSolutionStepValue(data_structure.MESH_DISPLACEMENT, 0,  [data+random.uniform(0,0.99) for data in data_value])

    def PrintInfo(self):
        #cs_tools.PrintInfo( self.data_map.keys() )
        cs_tools.PrintInfo(cs_tools.bcolors.BLUE+"This is an object of PipeFlowCoSimulationSolver with name", self.name + cs_tools.bcolors.ENDC)

    def AdvanceInTime(self, current_time):
        self.PrintInfo()
        return current_time + self.GetDeltaTime()

    def GetDeltaTime(self):
        return 0.2

    def SolveSolutionStep(self):
        ## Setting initial value for the data
        data_value = [PipeFlowCoSimulationSolver.pipe_flow_solver_count, PipeFlowCoSimulationSolver.pipe_flow_solver_count, PipeFlowCoSimulationSolver.pipe_flow_solver_count]
        for node in self.pipe_flow_model_part.Nodes:
            for data_name in self.data_map.keys():
                data_obj = data_structure.KratosGlobals.GetVariable(data_name)
                node.SetSolutionStepValue(data_name, 0, [data+random.uniform(0,0.99) for data in data_value])

        for node in self.pipe_flow_model_part.Nodes:
            for data_name in self.data_map.keys():
                data_obj = data_structure.KratosGlobals.GetVariable(data_name)

    def Check(self):
        cs_tools.PrintInfo(cs_tools.bcolors.GREEN+"Check from pipe_flow co simulation solver", "CHECKED !"+cs_tools.bcolors.ENDC)

    def _GetIOName(self):
        return "pipe_flow"

    def _SetUpModelpart(self):
        #create nodes
        self.pipe_flow_model_part.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
        self.pipe_flow_model_part.CreateNewNode(18, 1.00000, 0.00000, 0.00000)
        #create elements
