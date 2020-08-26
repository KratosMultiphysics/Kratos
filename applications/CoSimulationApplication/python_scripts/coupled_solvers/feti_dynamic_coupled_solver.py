# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver

def Create(settings, models, solver_name):
    return FetiDynamicCoupledSolver(settings, models, solver_name)

class FetiDynamicCoupledSolver(CoSimulationCoupledSolver):
    def __init__(self, settings, models, solver_name):
        super().__init__(settings, models, solver_name)

        #self.mapper = self.__GetDataTransferOperator("coupling_geometry")
        print(self.data_transfer_operators_dict)
        self.mapper = self.data_transfer_operators_dict["mapper"]
        print(self.mapper)
        self.modelpart_interface_origin = self.mapper.GetInterfaceModelPart(0)
        self.modelpart_interface_destination = self.mapper.GetInterfaceModelPart(1)

        # Get time integration parameters
        origin_newmark_beta = self.settings["origin_newmark_beta"].GetDouble()
        origin_newmark_gamma = self.settings["origin_newmark_gamma"].GetDouble()
        destination_newmark_beta = self.settings["destination_newmark_beta"].GetDouble()
        destination_newmark_gamma = self.settings["destination_newmark_gamma"].GetDouble()

        # Create feti class instance
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

    @classmethod
    def _GetDefaultSettings(cls):
        this_defaults = KM.Parameters("""{
            "origin_newmark_beta" : -1.0,
            "origin_newmark_gamma" : -1.0,
            "destination_newmark_beta" : -1.0,
            "destination_newmark_gamma" : -1.0
        }""")
        this_defaults.AddMissingParameters(super()._GetDefaultSettings())

        return this_defaults
