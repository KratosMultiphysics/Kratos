import KratosMultiphysics 
from KratosMultiphysics.assign_master_slave_constraints_to_neighbours_process import AssignMasterSlaveConstraintsToNeighboursProcess

def Factory(settings, model):
    KratosMultiphysics.Logger.PrintWarning("assign_mpcs_to_neighbours_process", "This process is deprecated. Please use 'assign_master_slave_constraints_to_neighbours_process' instead.")
    return AssignMPCsToNeighboursProcess(model, settings["Parameters"])

class AssignMPCsToNeighboursProcess(AssignMasterSlaveConstraintsToNeighboursProcess):
    def __init__(self, model, settings):
        # Check for deprecated settings and update them
        
        # The setting 'assign_mpcs_every_time_step' was renamed to 'reform_constraints_at_each_step'
        # This block ensures backward compatibility by updating the old setting to the new one and issuing a warning.
        if settings.Has("assign_mpcs_every_time_step"):
            KratosMultiphysics.Logger.PrintWarning("AssignMPCsToNeighboursProcess", "'assign_mpcs_every_time_step' is deprecated. Please use 'reform_constraints_at_each_step' instead.")
            settings.AddValue("reform_constraints_at_each_step", settings["assign_mpcs_every_time_step"])
            settings.RemoveValue("assign_mpcs_every_time_step")
        
        super().__init__(model, settings)
