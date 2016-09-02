from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.ContactMechanicsApplication as KratosContact
KratosMultiphysics.CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ParametricWallsProcess(Model, settings["Parameters"])


class ParametricWallsProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)
        
        self.main_model_part = Model[custom_settings["model_part_name"].GetString()]
    
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"         : "Solid Domain",        
            "search_control_type"     : "step",
            "search_frequency"        : 1.0,
            "parametric_walls"        : []
        }
        """)
 
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level        = 1
        self.domain_size       = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        self.search_frequency  = self.settings["search_frequency"].GetDouble()
        
        self.search_control_is_time = False
        search_control_type  = self.settings["search_control_type"].GetString()
        if(search_control_type == "time"):
            self.search_control_is_time = True
        elif(search_control_type == "step"):
            self.search_control_is_time = False

        #construct meshing domains
        self.parametric_walls = []
        walls_list = self.settings["parametric_walls"]
        self.number_of_walls = walls_list.size()
        for i in range(0,self.number_of_walls):
            item = walls_list[i]
            parametric_wall_module = __import__(item["python_module"].GetString())
            wall = parametric_wall_module.CreateParametricWall( self.main_model_part, item)
            self.parametric_walls.append(wall)

        # mesh modeler initial values
        self.search_contact_active = False
        if( self.number_of_walls ):
            self.search_contact_active = True

        self.step_count   = 1
        self.counter      = 1
        self.next_search  = 0.0
                       
    #
    def ExecuteInitialize(self):

        # check restart
        self.restart = False
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] ):
            self.restart = True
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
            
            if self.search_control_is_time:
                self.next_search  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] 
            else:
                self.next_search = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        # execute initialize base class
        if( self.main_model_part.ProcessInfo[KratosPfemBase.INITIALIZED_DOMAINS] == False ):
            self.InitializeDomains()

        for wall in self.parametric_walls:
            wall.Initialize();


    #
    def InitializeDomains(self):

        # initialize the modeler 
        print("::[Walls_Process]:: Initialize Domains ")
            
        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()
        
        # find node neighbours
        domain_utils.SearchNodeNeighbours(self.main_model_part, self.echo_level)
            
        # find element neighbours
        domain_utils.SearchElementNeighbours(self.main_model_part, self.echo_level)
            
        # set neighbour search performed
        neighbour_search_performed = True
               
        self.main_model_part.ProcessInfo.SetValue(KratosPfemBase.INITIALIZED_DOMAINS, True)

        print(self.main_model_part)

    ###

    #
    def ExecuteInitializeSolutionStep(self):

        self.step_count += 1

        #clean all contacts from main_model_part
        for wall in self.parametric_walls:
            wall.InitializeSearch();

        #build all contacts :: when building check if the condition exists and clone it
        if(self.search_contact_active):
            if(self.IsSearchStep()):
                self.SearchContact()


    #
    def ExecuteFinalizeSolutionStep(self):


        pass


    ###

    #
    def SearchContact(self):


        if( self.echo_level > 0 ):
            print("::[Walls_Process]:: CONTACT SEARCH...( call:", self.counter,")")
            
        self.wall_contact_model= KratosContact.ClearPointContactConditions(self.main_model_part, self.echo_level)

        self.wall_contact_model.ExecuteInitialize()

        for wall in self.parametric_walls:
            wall.ExecuteSearch();
            
        self.wall_contact_model.ExecuteFinalize()

        self.counter += 1 


        # schedule next search
        if(self.search_frequency > 0.0): # note: if == 0 always active
            if(self.search_control_is_time):
                while(self.next_search <= time):
                    self.next_search += self.search_frequency
            else:
                while self.next_search <= self.step_count:
                    self.next_search += self.search_frequency



    #
    def GetSearchStep(self):
        return self.counter

    #
    def IsSearchStep(self):

        if(self.search_control_is_time):
            #print( str(self.main_model_part.ProcessInfo[TIME])+">"+ str(self.next_search) )
            return ( self.main_model_part.ProcessInfo[TIME] > self.next_search )
        else:
            return ( self.step_count >= self.next_search )
