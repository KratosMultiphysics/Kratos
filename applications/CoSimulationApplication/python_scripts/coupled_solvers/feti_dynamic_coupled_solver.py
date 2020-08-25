# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

def Create(settings, models, solver_name):
    return FetiDynamicCoupledSolver(settings, models, solver_name)

class FetiDynamicCoupledSolver(CoSimulationCoupledSolver):
    def __init__(self, settings, models, solver_name):
        super().__init__(settings, models, solver_name)

        self.mapper = self.__GetDataTransferOperator("coupling_geometry")
        self.modelpart_interface_origin = self.mapper.GetInterfaceModelPart(0)
        self.modelpart_interface_destination = self.mapper.GetInterfaceModelPart(1)

        self.feti_coupling_solver = FetiDynamicCouplingUtilities(
            self.modelpart_interface_origin,
            self.modelpart_interface_destination)

        self.feti_coupling_solver.SetMappingMatrix(self.mapper.pGetMappingMatrix())






    def SolveSolutionStep(self):
        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.InitializeCouplingIteration()

        solver_index = 0
        for solver_name, solver in self.solver_wrappers.items():
            self._SynchronizeInputData(solver_name)
            solver.SolveSolutionStep()

            system_matrix = solver.get_mechanical_solution_strategy().GetSystemMatrix()
            self.feti_coupling_solver.SetEffectiveStiffnessMatrices(system_matrix,solver_index)

            self._SynchronizeOutputData(solver_name)

            solver_index += 1

        self.feti_coupling_solver.EquilibrateDomains()

        for coupling_op in self.coupling_operations_dict.values():
            coupling_op.FinalizeCouplingIteration()

        return True
