# Importing the Kratos Library
import KratosMultiphysics

# Importing the base class
import KratosMultiphysics.ConvectionDiffusionApplication.coupled_structural_thermal_solver as BaseThermoMech

def CreateSolver(main_model_part, custom_settings):
    return CoupledThermoMechanicalSubsteppingSolver(main_model_part, custom_settings)

class CoupledThermoMechanicalSubsteppingSolver(BaseThermoMech.CoupledThermoMechanicalSolver):
    """ A coupled one-way thermo-mechanical solver in which the time step of the two solvers is independent. It is assumed that the dt of the structural part is equal or lower than the thermal one.
    """

    def __init__(self, model, custom_settings):

        super(CoupledThermoMechanicalSubsteppingSolver, self).__init__(model, custom_settings)

        self.var_utils = KratosMultiphysics.VariableUtils()

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.thermal_solver.AdvanceInTime(current_time)
        return new_time


    def SolveSolutionStep(self):

        self.thermal_solver.InitializeSolutionStep()
        self.thermal_solver.Predict()

        KratosMultiphysics.Logger.PrintInfo("\t" + "Solving THERMAL part...")
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()

        current_time = self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]

        dt_solid = self.structural_solver.ComputeDeltaTime()
        dt_thermal = self.thermal_solver.ComputeDeltaTime()

        structural_time = current_time - dt_thermal

        nodes = self.GetComputingModelPart().Nodes

        current_temperature = self.var_utils.GetSolutionStepValuesVector(nodes, KratosMultiphysics.TEMPERATURE, 0)
        old_temperature = self.var_utils.GetSolutionStepValuesVector(nodes, KratosMultiphysics.TEMPERATURE, 1)
        increment_temperature = current_temperature - old_temperature

        # Substepping loop
        while structural_time < current_time:
            structural_time += dt_solid

            if structural_time > current_time:
                structural_time = current_time
            
            KratosMultiphysics.Logger.PrintInfo("\t" + "Solving substepping STRUCTURAL part, Structural TIME: ", structural_time)

            # Now we interpolate the temperature field for this substep
            substep_temperature = old_temperature + increment_temperature * structural_time / current_time

            self.var_utils.SetSolutionStepValuesVector(nodes, KratosMultiphysics.TEMPERATURE, substep_temperature, 0)


            self.structural_solver.InitializeSolutionStep()

            self.structural_solver.Predict()

            solid_is_converged = self.structural_solver.SolveSolutionStep()

            self.RemoveConvectiveVelocity()

        return solid_is_converged and thermal_is_converged
