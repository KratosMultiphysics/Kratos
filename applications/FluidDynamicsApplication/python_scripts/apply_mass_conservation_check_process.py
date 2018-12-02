import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if( not isinstance(settings, KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMassConservationCheckProcess(Model, settings["Parameters"])

class ApplyMassConservationCheckProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                        : "default_model_part_name",
            "mass_computation_frequency"             : 20,
            "compare_to_initial_values"              : true,
            "write_to_log_file"                      : true
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.write_to_log = settings["write_to_log_file"].GetBool()
        self.mass_conservation_check_process = KratosFluid.MassConservationCheckProcess(self.fluid_model_part, settings)

        # writing first line in file
        if ( self.write_to_log ):
            with open("ApplyMassConservationCheckProcess.log", "w") as logFile:
                logFile.write( "positiveVolume" + "\t" + "negativeVolume" + "\n" )

    def Execute(self):

        # retrieve information if the values were updated
        updated = self.mass_conservation_check_process.GetUpdateStatus()

        if ( updated ):
            posVol = self.mass_conservation_check_process.GetPositiveVolume()
            negVol = self.mass_conservation_check_process.GetNegativeVolume()
            initPosVol = self.mass_conservation_check_process.GetInitialPositiveVolume()
            initNegVol = self.mass_conservation_check_process.GetInitialNegativeVolume()

            # managing the output to the console
            KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess", "Positive Volume = " + str(posVol) + "  ( initially " + str(initPosVol) + ")" )
            KratosMultiphysics.Logger.PrintInfo("ApplyMassConservationCheckProcess", "Negative Volume = " + str(negVol) + "  ( initially " + str(initNegVol) + ")" )
            KratosMultiphysics.Logger.Flush()

            # adds additional lines to the log file
            if ( self.write_to_log ):
                with open("ApplyMassConservationCheckProcess.log", "a+") as logFile:
                    logFile.write( str(posVol) + "\t" + str(negVol) + "\n" )


    def ExecuteInitialize(self):
        self.mass_conservation_check_process.ExecuteInitialize()


    def ExecuteBeforeSolutionLoop(self):
        self.mass_conservation_check_process.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.mass_conservation_check_process.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):

        self.Execute()