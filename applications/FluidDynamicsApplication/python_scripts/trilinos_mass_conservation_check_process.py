from __future__ import absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.mpi as KratosMPI                          # MPI-python interface

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication","MetisApplication","TrilinosApplication")

# Import applications
import KratosMultiphysics.MetisApplication as KratosMetis           # Partitioning
import KratosMultiphysics.TrilinosApplication as KratosTrilinos     # MPI solvers
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid   # Fluid dynamics application

# Import serial monolithic embedded solver
import navier_stokes_two_fluids_solver
import trilinos_import_model_part_utility
import apply_mass_conservation_check_process

def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return TrilinosApplyMassConservationCheckProcess( Model, settings["Parameters"] )

class TrilinosApplyMassConservationCheckProcess(apply_mass_conservation_check_process.ApplyMassConservationCheckProcess):


    def __init__(self, Model, settings):

        apply_mass_conservation_check_process.ApplyMassConservationCheckProcess.__init__(self, Model, settings)
        self._is_printing_rank = (KratosMPI.mpi.rank == 0)

        self._my_log_file = "TrilinosApplyMassConservationCheckProcess.log"

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                        : "default_model_part_name",
            "mass_computation_frequency"             : 20,
            "compare_to_initial_values"              : true,
            "write_to_log_file"                      : true
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self._fluid_model_part = Model[settings["model_part_name"].GetString()]
        self._write_to_log = settings["write_to_log_file"].GetBool()
        self.MassConservationCheckProcess = KratosFluid.MassConservationCheckProcess(self.fluid_model_part, settings)

        # Construct the communicator
        self.EpetraCommunicator = KratosTrilinos.CreateCommunicator()

        # writing first line in file
        if ( self._write_to_log and self._is_printing_rank ):
            with open(self._my_log_file, "w") as logFile:
                logFile.write( "positiveVolume" + "\t" + "negativeVolume" + "\n" )
                logFile.close()

        if self._is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("TrilinosApplyMassConservationCheckProcess","Construction of Trilinos Mass Conservation Check Process finished.")


    def ExecuteFinalizeSolutionStep(self):

        # retrieve information if the values were updated
        updated = int( self.MassConservationCheckProcess.GetUpdateStatus() )
        KratosMPI.mpi.world.barrier()

        allUpdated = self._fluid_model_part.GetCommunicator().SumAll( updated )

        if ( allUpdated == self.EpetraCommunicator.NumProc() ):
            posVol = self.MassConservationCheckProcess.GetPositiveVolume()
            negVol = self.MassConservationCheckProcess.GetNegativeVolume()
            initPosVol = self.MassConservationCheckProcess.GetInitialPositiveVolume()
            initNegVol = self.MassConservationCheckProcess.GetInitialNegativeVolume()

            # barrier to make sure that all processes have finished the computation
            KratosMPI.mpi.world.barrier()

            posVol = self._fluid_model_part.GetCommunicator().SumAll( posVol )
            negVol = self._fluid_model_part.GetCommunicator().SumAll( negVol )
            initPosVol = self._fluid_model_part.GetCommunicator().SumAll( initPosVol )
            initNegVol = self._fluid_model_part.GetCommunicator().SumAll( initNegVol )

            # syntax example: self.outlet_model_part.GetCommunicator().SumAll(outlet_avg_vel_norm)

            if self._is_printing_rank:
                # managing the output to the console
                KratosMultiphysics.Logger.PrintInfo("TrilinosApplyMassConservationCheckProcess", "Positive Volume = " + str(posVol) + "  ( initially " + str(initPosVol) + ")" )
                KratosMultiphysics.Logger.PrintInfo("TrilinosApplyMassConservationCheckProcess", "Negative Volume = " + str(negVol) + "  ( initially " + str(initNegVol) + ")" )
                KratosMultiphysics.Logger.Flush()

                # adds additional lines to the log file
                if ( self.write_to_log ):
                    with open(self._my_log_file, "a+") as logFile:
                        logFile.write( str(posVol) + "\t" + str(negVol) + "\n" )