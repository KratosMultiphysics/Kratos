# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_coupled_solver import CoSimulationCoupledSolver
from KratosMultiphysics.CoSimulationApplication.coupling_operations.compute_dem_momentum import ComputeDemMomentum

def Create(settings, models, solver_name):
    return DemFemVolumeCoupledSolver(settings, models, solver_name)

class DemFemVolumeCoupledSolver(CoSimulationCoupledSolver):
    def SolveSolutionStep(self):

    
        for coupling_op in self.coupling_operations_dict.values():
            #print(coupling_op)
            if isinstance(coupling_op,ComputeDemMomentum):
                 coupling_op.InitializeCouplingIteration() # calculate momentum on dem side 
                 for solver_name, solver in self.solver_wrappers.items():
                     if solver_name == "structure":
                         self._SynchronizeInputData(solver_name) # send the momentum and mass  to fem
                         #print("mass and momentum sent to fem")
            else:
                coupling_op.InitializeCouplingIteration() # calculate nodal coupling forces
                #print("else part is executed")


        for solver_name, solver in self.solver_wrappers.items():  # for dem this does nothing, for fem it  maps the forces to the particles
            self._SynchronizeOutputData(solver_name)

        for coupling_op in self.coupling_operations_dict.values(): # multiplies external applied force by mass  and or adds load to bottom particles   
              coupling_op.FinalizeCouplingIteration()

        for solver_name, solver in self.solver_wrappers.items(): # runs both solvers
            solver.SolveSolutionStep()

        return True
