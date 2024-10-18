
import KratosMultiphysics
import KratosMultiphysics.PFEM2Application
import KratosMultiphysics.FluidDynamicsApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

from fluid_dynamics_analysis import FluidDynamicsAnalysis

class PFEM2FluidDynamicsAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers in PFEM2 framework.'''

    def __init__(self,model,parameters):
        super(PFEM2FluidDynamicsAnalysis,self).__init__(model,parameters)

    def _CreateSolver(self):
        solver_module = __import__(self.project_parameters["solver_settings"]["solver_type"].GetString())
        solver = solver_module.CreateSolver(self.model, self.project_parameters["solver_settings"])
        return solver

    def _GetSimulationName(self):
        return "PFEM2 Fluid Dynamics Analysis"

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python fluid_dynamics_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python fluid_dynamics_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = PFEM2FluidDynamicsAnalysis(model,parameters)
    simulation.Run()
