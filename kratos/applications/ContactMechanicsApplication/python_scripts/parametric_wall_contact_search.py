from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateContactSearch(main_model_part, wall_bounding_box, custom_settings):
    return ContactSearch(main_model_part, wall_bounding_box, custom_settings)

class ContactSearch(object):

    #
    def __init__(self, main_model_part, wall_bounding_box, custom_settings):

        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "python_file_name": "parametric_wall_contact_search",
             "search_control_type": "step",
             "search_frequency": 0,
             "contact_parameters":{
                 "contact_condition_type": "PointContactCondition2D1N",
                 "friction_active": false,
                 "friction_law_type": "MorhCoulomb",
                 "variables_of_properties":{
                      "MU_STATIC": 0.3,
                      "MU_DYNAMIC": 0.2,
                      "PENALTY_PARAMETER": 1000
                   }
             }
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)               

        self.search_frequency = self.settings["search_frequency"].GetDouble()

        self.search_control_is_time = False
        serch_control_type = self.settings["search_control_type"].GetString()
        if(search_control_type == "time"):
            self.search_control_is_time = True
        elif(search_control_type == "step"):
            self.search_control_is_time = False

        self.step_count  = 1
        self.next_search = 0.0

        self.search_process = KratosContact.ParametricWallContactSearch( self.main_model_part, wall_bounding_box, CustomSettings["contact_parameters"])

        print("Construction of the Search on Parametric Walls finished")
        
    #
    def Initialize(self):
        pass

  
    #
    def InitializeSearch(self):
        pass

    #
    def FinalizeSearch(self):
        
        # schedule next search
        if(self.search_frequency >= 0.0):
            if(self.search_control_is_time):
                while(self.next_search <= time):
                    self.next_search += self.search_frequency
            else:
                while(self.next_search <= self.step_count):
                    self.next_search += self.search_frequency
                
        
    #
    def ExecuteSearch(self):

        self.serch_process.Execute()

        
