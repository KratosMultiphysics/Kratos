from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
from importlib import import_module

class ProcessHandler(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        self.model = Model

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"                  : 0,
            "constraints_process_list"    : [],
            "loads_process_list"          : [],
            "problem_process_list"        : [],
            "output_process_list"         : [],
            "check_process_list"          : [],
            "processes_sub_model_part_tree_list" : []
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #build and sort the list of processes
        self.list_of_processes = []
        self.Sort()

        #set echo_level
        self.echo_level = self.settings["echo_level"].GetInt()

    #
    def GetVariables(self):
        nodal_variables = []
        for process in self.list_of_processes:
            try:
                nodal_variables = nodal_variables + process.GetVariables()
            except AttributeError:
                # the process does not have GetVariables()
                pass
        return nodal_variables

    #
    def ExecuteInitialize(self):
        for process in self.list_of_processes:
            process.ExecuteInitialize()

    #
    def ExecuteBeforeSolutionLoop(self):
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

    #
    def ExecuteInitializeSolutionStep(self):
        for process in self.list_of_processes:
            process.ExecuteInitializeSolutionStep()

    #
    def ExecuteFinalizeSolutionStep(self):
        for process in self.list_of_processes:
            process.ExecuteFinalizeSolutionStep()

    #
    def ExecuteBeforeOutputStep(self):
        for process in self.list_of_processes:
            process.ExecuteBeforeOutputStep()

    #
    def ExecuteAfterOutputStep(self):
        for process in self.list_of_processes:
            process.ExecuteAfterOutputStep()


    def Sort(self):

        #print("::[--Process_Handler--]:: Create Process List -START-")
        #build sorted list of processes

        # Assumptions: (a) each supplied list is separated and has no intersections
        #              (b) processes_sub_model_part_tree_list sets a hierarchy for process sorting
        #              (c) if no list is supplied the resultant list of processes is the one given by the input

        #sort processes using groups category and order:
        if( self.settings["processes_sub_model_part_tree_list"].size() > 0 ):

            print("::[--Process_Handler--]:: Sorting Loads and Constraints")

            #set group category
            self.Categorize()
            #set group unic order
            self.Order()

            # sorted constraints list
            self.list_of_processes  = self.ConstructSortedList( self.ConstraintsList() )
            # sorted loads list
            self.list_of_processes += self.ConstructSortedList( self.LoadsList() )

        else:

            # sorted constraints list
            self.list_of_processes  = self.ConstructList( self.ConstraintsList() )
            # sorted loads list
            self.list_of_processes += self.ConstructList( self.LoadsList() )


        # sorted problems list
        self.list_of_processes += self.ConstructList( self.ProcessList(self.settings["problem_process_list"]) )
        # sorted restart list
        self.list_of_processes += self.ConstructList( self.ProcessList(self.settings["output_process_list"]) )
        # sorted check list
        self.list_of_processes += self.ConstructList( self.ProcessList(self.settings["check_process_list"]) )

        #print("::[--Process_Handler--]:: Create Process List -END-")
        print("::[--Process_Handler--]:: Process List Ready")


    #
    def ProcessList(self, process_list):

        full_processes_list = []
        for i in range(0,process_list.size()):
            full_processes_list.append(process_list[i])
        return full_processes_list


    #
    def ConstraintsList(self):

        constraints_list = self.settings["constraints_process_list"]
        full_constraints_list = []
        for i in range(0,constraints_list.size()):
            if( constraints_list[i].Has("Parameters") == False ):
                full_constraints_list.append(self.ConstraintsFactory(constraints_list[i]))
            else:
                full_constraints_list.append(constraints_list[i])
            #print(full_constraints_list[i].PrettyPrintJsonString())
        return full_constraints_list


    #
    def LoadsList(self):

        loads_list = self.settings["loads_process_list"]
        full_loads_list = []
        for i in range(0,loads_list.size()):
            if( loads_list[i].Has("Parameters") == False ):
                full_loads_list.append(self.LoadsFactory(loads_list[i]))
            else:
                full_loads_list.append(loads_list[i])
            #print(full_loads_list[i].PrettyPrintJsonString())
        return full_loads_list

    #
    def ConstructProcess(self, process):

        python_module = import_module(process["kratos_module"].GetString() + "." + process["python_module"].GetString())
        return(python_module.Factory(process, self.model))

    #
    def ConstructList(self, process_list):

        #build list
        unsorted_process_list = []
        #add not added sorted processes
        for process in process_list:
            unsorted_process_list.append(self.ConstructProcess(process))

        return unsorted_process_list


    #
    def ConstructSortedList(self, process_list):

        constructed_process = []
        for i in range(0,len(process_list)):
            constructed_process.append(False)

        #build sorted list
        sorted_process_list   = []
        for j in self.list_of_sub_model_parts:
            for process in process_list:
                if( process["Parameters"].Has("model_part_name") ):
                    model_part_name = process["Parameters"]["model_part_name"].GetString()
                    if( j == model_part_name ):
                        sorted_process_list.append(self.ConstructProcess(process))
                        constructed_process[i] = True
                else:
                    print("WARNING: trying to sort a process with no model_part name")

        #add not added sorted processes
        for i in range(0,len(process_list)):
            if( constructed_process[i] == False ):
                sorted_process_list.append(self.ConstructProcess(process_list[i]))

        return sorted_process_list

    #
    def Categorize(self):

        processes_sub_model_parts = self.settings["processes_sub_model_part_tree_list"]

        #set number of categories
        number_of_categories = 1
        for i in range(0,processes_sub_model_parts.size()):
            group = processes_sub_model_parts[i].GetString()
            group_tree = group.split('/')
            if( number_of_categories < len(group_tree) ):
                number_of_categories = len(group_tree)

        self.category_lists = []
        for i in range(0,number_of_categories):
            self.category_lists.append([])

        #set submodelparts in categories
        for i in range(0,processes_sub_model_parts.size()):
            group = processes_sub_model_parts[i].GetString()
            group_tree = group.split('/')
            self.category_lists[len(group_tree)-1].append(group_tree[len(group_tree)-1])

    #
    def Order(self):

        self.list_of_sub_model_parts = []

        processes_sub_model_parts = self.settings["processes_sub_model_part_tree_list"]
        #set unique order
        for category in self.category_lists:
            for i in range(0,len(category)):
                self.list_of_sub_model_parts.append(category[i])


        print(" order", self.list_of_sub_model_parts)


    #
    def ConstraintsFactory(self,custom_settings):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
           "model_part_name"     : "MODEL_PART",
           "variable_name"       : "DISPLACEMENT",
           "value"               : [0.0,0.0,0.0],
           "compound_assignment" : "direct",
           "constrained"         : true,
           "interval"            : [0.0,"End"]
        }
        """)

        ##overwrite the default settings with user-provided parameters
        settings = custom_settings
        settings.ValidateAndAssignDefaults(default_settings)


        ##settings string in json format
        factory_settings = KratosMultiphysics.Parameters("""
        {
           "python_module" : "assign_vector_components_to_nodes_process",
           "kratos_module" : "KratosMultiphysics.PfemFluidDynamicsApplication",
           "process_name"  : "AssignVectorComponentsToNodesProcess",
           "Parameters"    : {
           }
        }
        """)

        factory_settings["Parameters"].AddValue("model_part_name", settings["model_part_name"])
        factory_settings["Parameters"].AddValue("variable_name", settings["variable_name"])
        factory_settings["Parameters"].AddValue("value", settings["value"])
        factory_settings["Parameters"].AddValue("constrained", settings["constrained"])
        factory_settings["Parameters"].AddValue("interval", settings["interval"])
        factory_settings["Parameters"].AddValue("compound_assignment", settings["compound_assignment"])

        return factory_settings

    #
    def LoadsFactory(self,custom_settings):

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
           "python_module"       : "assign_modulus_and_direction_to_conditions_process",
           "model_part_name"     : "MODEL_PART",
           "variable_name"       : "DISPLACEMENT",
           "modulus"             : 0.0,
           "direction"           : [0.0,0.0,0.0],
           "compound_assignment" : "direct",
           "interval"            : [0.0,"End"]
        }
        """)

        #trick to allow "value" to be a string or a double value
        if(custom_settings.Has("modulus")):
            if(custom_settings["modulus"].IsString()):
                default_settings["modulus"].SetString("0.0")


        ##overwrite the default settings with user-provided parameters
        settings = custom_settings
        settings.ValidateAndAssignDefaults(default_settings)


        ##settings string in json format
        factory_settings = KratosMultiphysics.Parameters("""
        {
           "kratos_module" : "KratosMultiphysics.PfemFluidDynamicsApplication",
           "process_name"  : "AssignModulusAndDirectionToConditionsProcess",
           "Parameters"    : {
           }
        }
        """)

        factory_settings.AddValue("python_module", settings["python_module"])
        factory_settings["Parameters"].AddValue("model_part_name", settings["model_part_name"])
        factory_settings["Parameters"].AddValue("variable_name", settings["variable_name"])
        factory_settings["Parameters"].AddValue("modulus", settings["modulus"])
        factory_settings["Parameters"].AddValue("direction", settings["direction"])
        factory_settings["Parameters"].AddValue("interval", settings["interval"])
        factory_settings["Parameters"].AddValue("compound_assignment", settings["compound_assignment"])

        return factory_settings
