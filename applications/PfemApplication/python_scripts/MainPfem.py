from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.PfemApplication
import MainSolid

class PfemSolution(MainSolid.Solution):

    def __init__(self, file_parameters = "ProjectParameters.json", file_name = None):

        super(PfemSolution, self).__init__(file_parameters,file_name)

    #### Main internal methods ####

    def _get_processes_parameters(self):

        # add fluid processes
        add_fluid_process = True
        if self.ProjectParameters.Has("problem_data"):
            if self.ProjectParameters["problem_data"].Has("domain_type"):
                if(self.ProjectParameters["problem_data"]["domain_type"].GetString() != "Solid"):
                    add_fluid_process = False

        if add_fluid_process is True:
            return self._add_fluid_processes()
        else:
            return MainSolid.Solution._get_processes_parameters(self)

    def _add_fluid_processes(self):

        # get processes parameters from base class
        processes_parameters = MainSolid.Solution._get_processes_parameters(self)

        # add process to manage assignation of material properties to particles
        # modify processes_parameters to introduce this process in the problem_process_list
        # particles concept : assign initial material percent and properties vector pointer to nodes
        if(processes_parameters.Has("problem_process_list")):
            problem_processes = processes_parameters["problem_process_list"]
            #print(" PROBLEM_PROCESSES ", processes_parameters["problem_process_list"].PrettyPrintJsonString())
            extended_problem_processes = self._set_particle_properties_process(problem_processes)
            processes_parameters.AddValue("problem_process_list", extended_problem_processes)
            #extended_problem_processes = self._set_volume_recovery_process(problem_processes)
            #processes_parameters.AddValue("problem_process_list", extended_problem_processes)
            #print(" EXTENDED_PROBLEM_PROCESSES ", processes_parameters["problem_process_list"].PrettyPrintJsonString())

        if(processes_parameters.Has("constraints_process_list")):
            constraints_processes = processes_parameters["constraints_process_list"]
            if(self.echo_level>1):
                print(" CONSTRAINTS_PROCESSES ", processes_parameters["constraints_process_list"].PrettyPrintJsonString())
            extended_constraints_processes = self._set_isolated_nodes_management_process(constraints_processes)
            processes_parameters.AddValue("constraints_process_list", extended_constraints_processes)
            extended_constraints_processes = self._set_selected_elements_management_process(constraints_processes)
            processes_parameters.AddValue("constraints_process_list", extended_constraints_processes)
            if(self.echo_level>1):
                print(" EXTENDED_CONSTRAINTS_PROCESSES ", processes_parameters["constraints_process_list"].PrettyPrintJsonString())

        if(processes_parameters.Has("loads_process_list")):
            loads_processes = processes_parameters["loads_process_list"]
            if(self.echo_level>1):
                print(" LOADS_PROCESSES ", processes_parameters["loads_process_list"].PrettyPrintJsonString())
            extended_loads_processes = self._set_volume_acceleration_process(loads_processes)
            processes_parameters.AddValue("loads_process_list", extended_loads_processes)
            if(self.echo_level>1):
                print(" EXTENDED_LOADS_PROCESSES ", processes_parameters["loads_process_list"].PrettyPrintJsonString())

        return processes_parameters

    def _set_isolated_nodes_management_process(self, constraints_processes):

        default_settings = KratosMultiphysics.Parameters("""
        {
             "python_module" : "manage_isolated_nodes_process",
             "kratos_module" : "KratosMultiphysics.PfemApplication",
             "Parameters"    : {}
        }
        """)

        model_part_name = self.model.GetMainModelPart().Name
        default_settings["Parameters"].AddEmptyValue("model_part_name").SetString(model_part_name)

        constraints_processes.Append(default_settings)

        return constraints_processes

    def _set_selected_elements_management_process(self, constraints_processes):

        default_settings = KratosMultiphysics.Parameters("""
        {
             "python_module" : "manage_selected_elements_process",
             "kratos_module" : "KratosMultiphysics.PfemApplication",
             "Parameters"    : {}
        }
        """)

        model_part_name = self.model.GetMainModelPart().Name
        default_settings["Parameters"].AddEmptyValue("model_part_name").SetString(model_part_name)

        constraints_processes.Append(default_settings)

        return constraints_processes

    def _set_volume_acceleration_process(self, loads_processes):

        default_settings = KratosMultiphysics.Parameters("""
        {
             "python_module" : "assign_modulus_and_direction_to_nodes_process",
             "kratos_module" : "KratosMultiphysics.SolidMechanicsApplication",
             "Parameters"    : {
                  "variable_name"   : "VOLUME_ACCELERATION",
                  "modulus"         : 9.81,
                  "direction"       : [0.0,-1.0,0.0]
            }
        }
        """)

        if(self.ProjectParameters.Has("problem_data")):
            if(self.ProjectParameters["problem_data"].Has("gravity_vector")):
                import math
                #get normalized direction
                direction   = []
                scalar_prod = 0
                for i in range(self.ProjectParameters["problem_data"]["gravity_vector"].size()):
                    direction.append( self.ProjectParameters["problem_data"]["gravity_vector"][i].GetDouble() )
                    scalar_prod = scalar_prod + direction[i]*direction[i]

                norm = math.sqrt(scalar_prod)

                self.value = []
                if( norm != 0.0 ):
                    for j in direction:
                        self.value.append( j/norm )
                else:
                    for j in direction:
                        self.value.append(0.0)

                if(default_settings["Parameters"].Has("modulus")):
                    default_settings["Parameters"]["modulus"].SetDouble(norm)

                if(default_settings["Parameters"].Has("direction")):
                    counter = 0
                    for i in self.value:
                        default_settings["Parameters"]["direction"][counter].SetDouble(i)
                        counter+=1

                model_part_name = self.model.GetMainModelPart().Name
                default_settings["Parameters"].AddEmptyValue("model_part_name").SetString(model_part_name)

                loads_processes.Append(default_settings)

        return loads_processes

    def _set_particle_properties_process(self, problem_processes):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "python_module" : "assign_properties_to_nodes_process",
            "kratos_module" : "KratosMultiphysics.PfemApplication",
            "Parameters"    : {
                  "fluid_mixture" : true,
                  "solid_mixture" : false
            }
        }
        """)

        model_part_name = self.model.GetMainModelPart().Name
        default_settings["Parameters"].AddEmptyValue("model_part_name").SetString(model_part_name)

        problem_processes.Append(default_settings)

        return problem_processes


    def _set_volume_recovery_process(self, problem_processes):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "python_module" : "volume_recovery_process",
            "kratos_module" : "KratosMultiphysics.PfemApplication",
            "Parameters"    : {
            }
        }
        """)

        model_part_name = self.model.GetMainModelPart().Name
        default_settings["Parameters"].AddEmptyValue("model_part_name").SetString(model_part_name)

        problem_processes.Append(default_settings)

        return problem_processes

    @classmethod
    def _class_prefix(self):
        header = "::[--PFEM Simulation--]::"
        return header

if __name__ == "__main__":
    PfemSolution().Run()
