from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics 

KratosMultiphysics.CheckForPreviousImport()


class ProcessHandler(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

        self.Model = Model
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "echo_level"                  : 0,
            "constraints_process_list"    : [],
            "loads_process_list"          : [],
            "problem_process_list"        : [],
            "output_process_list"         : [],
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

        #print list of constructed processes
        if(self.echo_level>1):
            for process in self.list_of_processes:
                print(process)
                        
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
            

    #
    def Sort(self):

        print("::[Process_Handler]:: Create Process List -START-")
        #build sorted list of processes

        # Assumptions: (a) each supplied list is separated and has no intersections
        #              (b) processes_sub_model_part_tree_list sets a hierarchy for process sorting
        #              (c) if no list is supplied the resultant list of processes is the one given by the input
        
        #sort processes using groups category and order:
        if( self.settings["processes_sub_model_part_tree_list"].size() > 0 ):

            print("::[Process_Handler]:: Sorting Loads and Constraints")
            
            #set group category
            self.Categorize()
            #set group unic order
            self.Order()

            # sorted constraints list
            self.list_of_processes  = self.ConstructSortedList( self.settings["constraints_process_list"] )
            # sorted loads list
            self.list_of_processes += self.ConstructSortedList( self.settings["loads_process_list"] )

        else:
            
            # sorted constraints list
            self.list_of_processes  = self.ConstructList( self.settings["constraints_process_list"] )
            # sorted loads list
            self.list_of_processes += self.ConstructList( self.settings["loads_process_list"] )

            
        # sorted problems list
        self.list_of_processes += self.ConstructList( self.settings["problem_process_list"] )
        # sorted restart list
        self.list_of_processes += self.ConstructList( self.settings["output_process_list"] )


        print("::[Process_Handler]:: Create Process List -END-")
        
    #
    def ConstructList(self, process_list):

        #build list
        unsorted_process_list = []
        #add not added sorted processes                
        for i in range(0,process_list.size()):
            kratos_module = __import__(process_list[i]["kratos_module"].GetString())
            python_module = __import__(process_list[i]["python_module"].GetString())
            unsorted_process_list.append(python_module.Factory(process_list[i], self.Model))

        return unsorted_process_list
        
    #
    def ConstructSortedList(self, process_list):

        constructed_process = []
        for i in range(0, process_list.size()):
            constructed_process.append(False)

        #build sorted list
        sorted_process_list   = []
        for j in self.list_of_sub_model_parts:
            for i in range(0,process_list.size()):
                if( process_list[i]["Parameters"].Has("model_part_name") ):
                    model_part_name = process_list[i]["Parameters"]["model_part_name"].GetString()
                    
                    if( j == model_part_name ):
                        kratos_module = __import__(process_list[i]["kratos_module"].GetString())
                        python_module = __import__(process_list[i]["python_module"].GetString())
                        sorted_process_list.append(python_module.Factory(process_list[i], self.Model))
                        constructed_process[i] = True
                else:
                    print("WARNING: trying to sort a process with no model_part name")

        #add not added sorted processes                
        for i in range(0,process_list.size()):
            if( constructed_process[i] == False ):
                kratos_module = __import__(process_list[i]["kratos_module"].GetString())
                python_module = __import__(process_list[i]["python_module"].GetString())
                sorted_process_list.append(python_module.Factory(process_list[i], self.Model))

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
