# co simulation imports
from base_co_simulation_classes.co_simulation_base_solver import CoSimulationBaseSolver
import co_simulation_tools as tools
data_structure = tools.cs_data_structure
import random


def Create(name, cosim_solver_settings):
    return DummyCoSimulationSolver(name, cosim_solver_settings)


class DummyCoSimulationSolver(CoSimulationBaseSolver):
    dummy_solver_count = 0
    def __init__(self, name, custom_settings):
        super(DummyCoSimulationSolver, self).__init__(name, custom_settings)
        self.name = name
        print("########################### constructor of :: ", self.name)
        self.dummy_model_part = self.model.CreateModelPart('dummy_model_part')
        for variable_name in self.data_list.keys():
            if(not data_structure.KratosGlobals.HasVariable(variable_name)):
                self.variable_obj = data_structure.Array1DVariable3(variable_name)
                self.dummy_model_part.AddNodalSolutionStepVariable(self.variable_obj)
                print("################### :: ", variable_name, " ########## :: ", self.variable_obj)
                print("################### :: ", data_structure.KratosGlobals.HasVariable(variable_name))
            else:
                self.variable_obj = data_structure.KratosGlobals.GetVariable(variable_name)
                self.dummy_model_part.AddNodalSolutionStepVariable(self.variable_obj)

        self._SetUpModelpart()
        DummyCoSimulationSolver.dummy_solver_count = DummyCoSimulationSolver.dummy_solver_count + 1

    def Initialize(self):
        self.InitializeIO()
        ## Setting initial value for the data
        data_value = [DummyCoSimulationSolver.dummy_solver_count, DummyCoSimulationSolver.dummy_solver_count, DummyCoSimulationSolver.dummy_solver_count]
        for node in self.dummy_model_part.Nodes:
            for data_name in self.data_list.keys():
                data_obj = data_structure.KratosGlobals.GetVariable(data_name)
                node.SetSolutionStepValue(data_obj, [data+random.uniform(0,0.99) for data in data_value])

    def PrintInfo(self):
        #print( self.data_list.keys() )
        print(tools.bcolors.BLUE+"This is an object of DummyCoSimulationSolver with name : ", self.name + tools.bcolors.ENDC)

    def AdvanceInTime(self, current_time):
        self.PrintInfo()
        return current_time + self.GetDeltaTime()

    def GetDeltaTime(self):
        return 0.2

    def SolveSolutionStep(self):
        ## Setting initial value for the data
        data_value = [DummyCoSimulationSolver.dummy_solver_count, DummyCoSimulationSolver.dummy_solver_count, DummyCoSimulationSolver.dummy_solver_count]
        for node in self.dummy_model_part.Nodes:
            for data_name in self.data_list.keys():
                data_obj = data_structure.KratosGlobals.GetVariable(data_name)
                node.SetSolutionStepValue(data_obj, [data+random.uniform(0,0.99) for data in data_value])
    def Check(self):
        print(tools.bcolors.GREEN+"Check from dummy co simulation solver : CHECKED !"+tools.bcolors.ENDC)

    def _GetIOName(self):
        return "dummy"

    def _SetUpModelpart(self):
        #create nodes
        self.dummy_model_part.CreateNewNode(1, 0.00000, 1.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(2, 0.00000, 0.50000, 0.00000)
        self.dummy_model_part.CreateNewNode(3, 0.50000, 1.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(4, 0.50000, 0.50000, 0.00000)
        self.dummy_model_part.CreateNewNode(5, 0.00000, 0.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(6, 1.00000, 1.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(7, 1.00000, 0.50000, 0.00000)
        self.dummy_model_part.CreateNewNode(8, 0.50000, 0.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(9, 1.00000, 0.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(10, 1.50000, 1.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(11, 1.50000, 0.50000, 0.00000)
        self.dummy_model_part.CreateNewNode(12, 1.50000, 0.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(13, 2.00000, 1.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(14, 2.00000, 0.50000, 0.00000)
        self.dummy_model_part.CreateNewNode(15, 2.00000, 0.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(16, 1.00000, 1.00000, 0.00000)
        self.dummy_model_part.CreateNewNode(17, 1.00000, 0.50000, 0.00000)
        self.dummy_model_part.CreateNewNode(18, 1.00000, 0.00000, 0.00000)
        #create elements