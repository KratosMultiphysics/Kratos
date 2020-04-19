import KratosMultiphysics
import KratosMultiphysics.TrilinosApplication as KratosTrilinos

import KratosMultiphysics.FluidDynamicsApplication.apply_mass_conservation_check_process as apply_mass_conservation_check_process
import KratosMultiphysics.TrilinosApplication.trilinos_linear_solver_factory as trilinos_linear_solver_factory


def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMassConservationCheckProcessMPI(Model, settings["Parameters"])


class ApplyMassConservationCheckProcessMPI(apply_mass_conservation_check_process.ApplyMassConservationCheckProcess):

    def _set_levelset_convection_process(self):
        ### for MPI application
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()
        mpi_settings = KratosMultiphysics.Parameters("""{
            "linear_solver_settings"   : {
                "solver_type" : "amgcl"
            }
        }""")

        self.trilinos_linear_solver = trilinos_linear_solver_factory.ConstructSolver(mpi_settings["linear_solver_settings"])

        if self._fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetForwardConvectionProcess2D(
                self.EpetraCommunicator,
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.trilinos_linear_solver)
        else:
            level_set_convection_process = KratosTrilinos.TrilinosLevelSetForwardConvectionProcess3D(
                self.EpetraCommunicator,
                KratosMultiphysics.DISTANCE,
                self._fluid_model_part,
                self.trilinos_linear_solver)

        return level_set_convection_process