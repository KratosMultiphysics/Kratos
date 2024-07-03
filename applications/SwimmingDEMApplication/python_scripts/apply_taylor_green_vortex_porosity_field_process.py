import os
import h5py

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as KratosSDEM

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTaylorGreenVortexPorosityFieldProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyTaylorGreenVortexPorosityFieldProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the fluid model part.
        settings -- Kratos parameters containing process settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {},
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())
        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),"sp_data.hdf5")
        self.fluid_kinetic_energy = []
        self.kinetic_energy_dissipation_rate = []
        self.time = []
        self.group_name = str(1)

        self.TaylorGreenVortexPorosityFieldProcess = KratosSDEM.TaylorGreenVortexPorosityFieldProcess(self.model_part, self.settings)

    def ExecuteBeforeSolutionLoop(self):
        self.TaylorGreenVortexPorosityFieldProcess.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.TaylorGreenVortexPorosityFieldProcess.ExecuteInitializeSolutionStep()
        step = self.model_part.ProcessInfo[KratosMultiphysics.STEP]
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if step == 1:
            data = self.TaylorGreenVortexPorosityFieldProcess.ExecuteInTimeStep()

            self.time.append(time)
            self.fluid_kinetic_energy.append(data[0])
            self.kinetic_energy_dissipation_rate.append(data[1])

            with h5py.File(self.file_path, 'a') as f:
                    self.WriteDataToFile(file_or_group = f,
                                names = ['FLUID_KINETIC_ENERGY', 'DISSIPATION', 'TIME'],
                                data = [self.fluid_kinetic_energy, self.kinetic_energy_dissipation_rate, self.time])

    def ExecuteFinalizeSolutionStep(self):
        data = self.TaylorGreenVortexPorosityFieldProcess.ExecuteInTimeStep()
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        self.time.append(time)
        self.fluid_kinetic_energy.append(data[0])
        self.kinetic_energy_dissipation_rate.append(data[1])

        with h5py.File(self.file_path, 'a') as f:
                self.WriteDataToFile(file_or_group = f,
                            names = ['FLUID_KINETIC_ENERGY', 'DISSIPATION', 'TIME'],
                            data = [self.fluid_kinetic_energy, self.kinetic_energy_dissipation_rate, self.time])

    def WriteDataToFile(self, file_or_group, names, data):
        if self.group_name in file_or_group:
            file_or_group['/'].__delitem__(self.group_name)
            self.sub_group = file_or_group.create_group(self.group_name)
        else:
            self.sub_group = file_or_group.create_group(self.group_name)

        for name, datum in zip(names, data):
            if name in file_or_group:
                file_or_group.__delitem__(name)
        for name, datum in zip(names, data):
            self.sub_group.create_dataset(name = name, data = datum)