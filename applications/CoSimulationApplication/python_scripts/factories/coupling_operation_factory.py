from KratosMultiphysics.CoSimulationApplication.factories import base_factory

def CreateCouplingOperation(coupling_operation_settings, solvers, parent_coupled_solver_process_info):
    """This function creates and returns the Coupling Operation used for CoSimulation"""
    return base_factory.Create(coupling_operation_settings, [solvers, parent_coupled_solver_process_info], "KratosMultiphysics.CoSimulationApplication.coupling_operations")
