from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
except ImportError:
    pass

from fluid_dynamics_analysis import FluidDynamicsAnalysis

class ALEFluidAnalysis(FluidDynamicsAnalysis):
    '''Main script for fluid dynamics simulations using the navier_stokes family of python solvers.'''

    def __init__(self, model, project_parameters):

        ale_solver_settings = project_parameters["solver_settings"]["ale_settings"]
        if ale_solver_settings.Has("ale_interface_parts"):
            print("Info, not automatized! ....")
        else:
            # Here the ale-boundary-conditions are extracted from the processes
            # and assigned to the solver such that the VELOCITY can be set to the
            # MESH_VELOCITY during solving
            ale_bc_settings = project_parameters["processes"]["ale_boundary_conditions"]

            ale_interface_parts = KratosMultiphysics.Parameters(""" [ ] """)
            ale_interface_parts = KratosMultiphysics.Parameters(""" [ [] , [] , [] ] """)

            for i_proc in range(ale_bc_settings.size()):
                process_name
                model_part_name = ale_bc_settings[i_proc]["Parameters"]["model_part_name"].GetString()
                for i_comp in range(3):
                    is_constrained = ale_bc_settings[i_proc]["Parameters"] ["constrained"][i_comp].GetBool()
                    if is_constrained:
                        ale_interface_parts[i_comp].Append(model_part_name)

            ale_solver_settings.AddValue("ale_interface_parts", ale_interface_parts)
            print(solver_settings.PrettyPrintJsonString())
            errrrr

        super(ALEFluidAnalysis, self).__init__(model, project_parameters)

    def _CreateSolver(self):
        import ale_fluid_solver
        return ale_fluid_solver.CreateSolver(self.model, self.project_parameters)

    def _GetSimulationName(self):
        return "ALE Fluid Analysis"

if __name__ == '__main__':
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default parameter file (assumed to be called "ProjectParameters.json"):\n'
        err_msg += '    "python ale_fluid_analysis.py"\n'
        err_msg += '- With custom parameter file:\n'
        err_msg += '    "python ale_fluid_analysis.py <my-parameter-file>.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        parameter_file_name = argv[1]
    else: # using default name
        parameter_file_name = "ProjectParameters.json"

    with open(parameter_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = ALEFluidAnalysis(model,parameters)
    simulation.Run()
