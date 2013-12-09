#importing the Kratos Library
from KratosMultiphysics import *
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
        self.young_modulus     = 1
        self.penalty_parameter = 1

        # movement parameters
        self.center            = Array3()
        self.velocity          = Array3()


    #######################################################################
    def Initialize(self,wall_configuration):
        
        self.rigid_wall_active = wall_configuration.rigid_wall

        size_scale = wall_configuration.size_scale

        self.tip_radius         = wall_configuration.tip_radius * size_scale
        self.rake_angle         = wall_configuration.rake_angle
        self.clearance_angle    = wall_configuration.clearance_angle
        self.young_modulus      = wall_configuration.young_modulus / ( size_scale * size_scale )
        self.penalty_parameter  = wall_configuration.penalty_parameter

        for size in range(0,3):
            self.center[size]   = wall_configuration.center[size] * size_scale
            self.velocity[size] = wall_configuration.velocity[size] * size_scale
        
        print " Center ", self.center
        print " Velocity ", self.velocity


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

