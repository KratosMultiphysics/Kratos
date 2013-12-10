#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()

class RigidWallUtility:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part  = model_part
        self.domain_size = domain_size 
        
        # wall parameters
        self.rigid_wall_active = False

        # geometry parameters
        self.tip_radius        = 0
        self.rake_angle        = 0
        self.clearance_angle   = 0
        
        # material parameters
        self.penalty_parameter = 1

        # movement parameters
        self.center            = Vector(3)
        self.velocity          = Vector(3)
        
        # rigid wall bounding box
        self.rigid_wall_bbox   = SpatialBoundingBox(self.center, self.tip_radius, self.velocity) 
        

    #######################################################################
    def Initialize(self,wall_configuration):
        
        self.rigid_wall_active = wall_configuration.rigid_wall

        size_scale = wall_configuration.size_scale

        self.tip_radius         = wall_configuration.tip_radius * size_scale
        self.rake_angle         = wall_configuration.rake_angle
        self.clearance_angle    = wall_configuration.clearance_angle
        self.penalty_parameter  = wall_configuration.penalty_parameter

        for size in range(0,3):
            self.center[size]   = wall_configuration.center[size] * size_scale
            self.velocity[size] = wall_configuration.velocity[size] * size_scale
        
        #print " Center ", self.center
        #print " Velocity ", self.velocity

        if( self.rigid_wall_active ):
            # rigid wall bounding box
            self.rigid_wall_bbox    = RigidToolBoundingBox(self.tip_radius, self.rake_angle, self.clearance_angle, self.center, self.velocity) 
            
            # rigid wall contact search process
            self.contact_search_process  = RigidWallContactSearch( self.rigid_wall_bbox, self.model_part) 


    #######################################################################
    def ExecuteContactSearch(self):
        if( self.rigid_wall_active ):
            # set properties for rigid wall conditions
            self.model_part.Properties[1].SetValue(PENALTY_PARAMETER,self.penalty_parameter)
            self.contact_search_process.ExecuteInitializeSolutionStep()
        
    #######################################################################
    def UpdatePosition(self):
        if( self.rigid_wall_active ):
            self.contact_search_process.ExecuteFinalizeSolutionStep()
                
    #######################################################################   
    def RigidWallVelocity(self):
        return self.velocity  

    #######################################################################   
    def RigidWallActive(self):
        return self.rigid_wall_active  

    #######################################################################   
    def RigidWallCenter(self):
        return self.center  

    #######################################################################   
    def RigidWallTipRadius(self):
        return self.tip_radius  

    #######################################################################   

