import KratosMultiphysics

from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory

from KratosMultiphysics.FluidDynamicsApplication.distance_reinitialization import DistanceReinitialization 

class TrilinosDistanceReinitialization(DistanceReinitialization):
    def __init__(self,model_part,params, epetra_communicator):
        self._epetra_communicator = epetra_communicator
        super().__init__(model_part, params)

    def _GetLinearSolver(self):
        linear_solver_configuration = self.params["linear_solver_settings"]
        return trilinos_linear_solver_factory.ConstructSolver(linear_solver_configuration)

    def _ConstructVariationalProcess(self, maximum_iterations, linear_solver, process_flag):
        if self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            return KratosMultiphysics.TrilinosVariationalDistanceCalculationProcess2D(
                self._epetra_communicator,
                self.model_part,
                linear_solver,
                maximum_iterations,
                process_flag)
        else:
            return KratosMultiphysics.TrilinosVariationalDistanceCalculationProcess3D(
                self.model_part,
                linear_solver,
                maximum_iterations,
                process_flag)
