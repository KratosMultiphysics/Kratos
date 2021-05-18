import KratosMultiphysics
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

import KratosMultiphysics.FluidDynamicsApplication.apply_two_fluids_inlet_process as apply_two_fluids_inlet_process
import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as trilinos_linear_solver_factory


def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTwoFluidsInletProcessMPI(Model, settings["Parameters"])


class ApplyTwoFluidsInletProcessMPI(apply_two_fluids_inlet_process.ApplyTwoFluidsInletProcess):

    def set_variational_distance_process(self):
        # Construct the variational distance calculation process
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()
        maximum_iterations = 5

        ### for MPI
        trilinos_settings = KratosMultiphysics.Parameters("""
        {
            "linear_solver_settings"   : {
                "solver_type" : "amgcl"
            }
        }
        """)

        self.linear_solver = trilinos_linear_solver_factory.ConstructSolver( trilinos_settings["linear_solver_settings"] )

        if self.complete_model.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosTrilinos.TrilinosVariationalDistanceCalculationProcess2D(
                self.EpetraCommunicator,
                self.complete_model,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE,
                "distance_from_inlet_2D")
        else:
            variational_distance_process = KratosTrilinos.TrilinosVariationalDistanceCalculationProcess3D(
                self.EpetraCommunicator,
                self.complete_model,
                self.linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE,
                "distance_from_inlet_3D")

        return variational_distance_process