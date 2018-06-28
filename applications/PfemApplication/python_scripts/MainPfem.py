from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import MainSolid

class PfemSolution(MainSolid.Solution):

    def __init__(self, file_parameters = "ProjectParameters.json", file_name = None):

        super(PfemSoludion, self).__init__(file_parameters,file_name)

        print("::[---PFEM Solution --]::")

    #### Main internal methods ####

    def _get_processes_parameters(self):
        # get processes parameters from base class
        processes_parameters = MainSolid.Solution._get_processes_parameters(self)
        
        # add process to manage assignation of material properties to particles
        # modify processes_parameters to introduce this process in the problem_process_list
        # particles concept : assign initial material percent and properties vector pointer to nodes
        if(processes_parameters.Has("problem_process_list")):
            problem_processes = processes_parametes["problem_process_list"]
            processes_parameters.__setitem__("problem_process_list", self._set_particle_properties_process(problem_processes))
        
        return processes_parameters

    def _set_particle_properties_process(self, problem_processes):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "python_module" : "assign_properties_to_particles",
            "kratos_module" : "KratosMultiphysics.PfemApplication",
            "Parameters"    : {
                  "fluid_mixture" : true,
                  "solid_mixture" : false
            }        
        }
        """)
        
        problem_processes.Append(default_settings)

        return problem_processes
        
    
if __name__ == "__main__":
    PfemSolution().Run()
