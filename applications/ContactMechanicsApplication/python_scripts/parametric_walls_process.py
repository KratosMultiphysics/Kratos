from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact
KratosMultiphysics.CheckForPreviousImport()

from multiprocessing import Pool

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ParametricWallsProcess(Model, settings["Parameters"])


class ParametricWallsProcess(KratosMultiphysics.Process):
    #
    def __init__(self, Model, custom_settings ):

        KratosMultiphysics.Process.__init__(self)

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
        self.search_frequency  = self.settings["search_frequency"].GetDouble()

        self.search_control_is_time = False
        search_control_type  = self.settings["search_control_type"].GetString()
        if(search_control_type == "time"):
            self.search_control_is_time = True
        elif(search_control_type == "step"):
            self.search_control_is_time = False


        self.step_count   = 1
        self.counter      = 1
        self.next_search  = 0.0

        self.Model = Model

    #
    def ExecuteInitialize(self):


        self.main_model_part = self.Model[self.settings["model_part_name"].GetString()]
        self.dimension         = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        #construct parametric wall domains
        self.parametric_walls = []
        walls_list = self.settings["parametric_walls"]
        self.number_of_walls = walls_list.size()
        for i in range(0,self.number_of_walls):
            item = walls_list[i]
            parametric_wall_module = __import__(item["python_module"].GetString())
            wall = parametric_wall_module.CreateParametricWall(self.main_model_part, item)
            self.parametric_walls.append(wall)

        # mesh mesher initial values
        self.search_contact_active = False
        if( self.number_of_walls ):
            self.search_contact_active = True

        # build parametric walls
        for i in range(0,self.number_of_walls):
            self.parametric_walls[i].BuildParametricWall()

        # check restart
        self.restart = False
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] ):
            self.restart = True
            self.step_count = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

            if self.search_control_is_time:
                self.next_search  = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            else:
                self.next_search = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        # initialize wall domains
        print("::[Walls_Process]:: Initialize Domains ")
        import domain_utilities
        domain_utils = domain_utilities.DomainUtilities()
        domain_utils.InitializeDomains(self.main_model_part,self.echo_level)

        for wall in self.parametric_walls:
            wall.Initialize();

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
    def ExecuteSearch(wall):
        wall.ExecuteSearch()

    #
    def SearchContact(self):


        if( self.echo_level > 0 ):
            print("::[Walls_Process]:: CONTACT SEARCH...( call:", self.counter,")")

        self.wall_contact_model= KratosContact.ClearPointContactConditions(self.main_model_part, self.echo_level)

        self.wall_contact_model.ExecuteInitialize()

        #serial
        for wall in self.parametric_walls:
            wall.ExecuteSearch();

        #parallel (not working pickling instances not enabled)
        #walls_number = len(self.parametric_walls)
        #if(walls_number>8):
        #    walls_number = 8

        #pool = Pool(walls_number)
        #pool.map(self.ExecuteSearch,self.parametric_walls)
        #pool.close()
        #pool.joint()


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


    #
    def GetVariables(self):
        nodal_variables = ['RIGID_WALL']
        nodal_variables = nodal_variables + ['CONTACT_FORCE','CONTACT_NORMAL']
        nodal_variables = nodal_variables + ['VOLUME_ACCELERATION']
        nodal_variables = nodal_variables + ['NORMAL', 'NODAL_H']
        return nodal_variables
