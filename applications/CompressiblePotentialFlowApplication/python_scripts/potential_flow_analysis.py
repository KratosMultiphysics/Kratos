from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as Kratos
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass
from analysis_stage import AnalysisStage

class PotentialFlowAnalysis(AnalysisStage):
        def __init__(self,model,parameters):
            super(PotentialFlowAnalysis,self).__init__(model,parameters)
        def _CreateSolver(self):
            import potential_flow_solver 
            return potential_flow_solver.CreateSolver(self.model,self.project_parameters)
        def _GetSimulationName(self):
            return "Potential Flow Analysis"
if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = PotentialFlowAnalysis(model,parameters)
    simulation.Run()