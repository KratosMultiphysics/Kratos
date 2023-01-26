from KratosMultiphysics import Logger, Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_solver as BaseSolver

class ModifiedSwimmingDEMSolver(BaseSolver.SwimmingDEMSolver):
    def __init__(self, model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager):
        super().__init__(model, project_parameters, field_utility, fluid_solver, dem_solver, variables_manager)

    def _ConstructProjectionModule(self):
        # change PeriodicDomainOption momentarily to trick the projection module

        self.project_parameters["dem_parameters"]["PeriodicDomainOption"].SetBool(False)
        projection_module = super()._ConstructProjectionModule()
        self.project_parameters["dem_parameters"]["PeriodicDomainOption"].SetBool(True)

        return projection_module
