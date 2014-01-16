#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()

class RigidWallUtility:
    #######################################################################
    def __init__(self,model_part,domain_size,wall_configuration):

        self.model_part  = model_part
        self.domain_size = domain_size 
        
        # wall parameters
        self.rigid_wall_active = wall_configuration.rigid_wall
        self.number_of_walls   = wall_configuration.number_of_walls

        # material parameters
        self.penalty_parameter = wall_configuration.penalty_parameter

        # geometry parameters
        size_scale = wall_configuration.size_scale
      
        self.wall_labels = wall_configuration.wall_labels
        self.noses_per_wall = Vector(self.number_of_walls)
        for sizek in range(0,self.number_of_walls):
            self.noses_per_wall[sizek] = 0

        for sizei in range(0,self.number_of_walls):
            for sizej in range(0,len(self.wall_labels)):
                if(self.wall_labels[sizej] == sizei+1):
                    self.noses_per_wall[sizei] = self.noses_per_wall[sizei] + 1
                    

        for sizei in range(0,self.number_of_walls):

            self.number_of_noses   = int(self.noses_per_wall[sizei])
            self.tip_radius        = Vector(self.number_of_noses)
            self.rake_angles       = Vector(self.number_of_noses)
            self.clearance_angles  = Vector(self.number_of_noses)
            self.convexities       = Vector(self.number_of_noses)
            self.tip_centers       = Matrix(self.number_of_noses, 3)

            self.wall_velocity         = Vector(3)
            self.wall_angular_velocity = Vector(3)
            self.wall_reference_point  = Vector(3)

            counter = 0
            
            for sizej in range(0,len(wall_configuration.wall_movement_labels)):
                if(wall_configuration.wall_movement_labels[sizej] == sizei+1):
                    for sizek in range(0,3):
                        self.wall_velocity[sizek] = wall_configuration.wall_velocity[sizej][sizek]

            for sizej in range(0,len(wall_configuration.wall_rotation_labels)):
                if(wall_configuration.wall_rotation_labels[sizej] == sizei+1):
                    for sizek in range(0,3):
                        self.wall_angular_velocity[sizek] = wall_configuration.wall_angular_velocity[sizej][sizek]
                        self.wall_reference_point[sizek] = wall_configuration.reference_point[sizej][sizek]


            for sizej in range(0,len(self.wall_labels)):

                if(self.wall_labels[sizej] == sizei+1):

                    self.tip_radius[counter] = wall_configuration.tip_radius[sizej] * size_scale
                    self.rake_angles[counter] = wall_configuration.rake_angles[sizej]
                    self.clearance_angles[counter] = wall_configuration.clearance_angles[sizej]
                    self.convexities[counter] = wall_configuration.nose_convexities[sizej]
                    for sizek in range(0,3):
                        self.tip_centers[counter,sizek] = wall_configuration.tip_centers[sizej][sizek] * size_scale
                    counter = counter + 1
            
            if( self.rigid_wall_active ):
                # rigid wall bounding box
                self.rigid_wall_bbox = []
                self.rigid_wall_bbox.append( RigidWallBoundingBox(self.convexities, self.tip_radius, self.rake_angles, self.clearance_angles, self.tip_centers, self.wall_velocity, self.wall_angular_velocity, self.wall_reference_point) )
            
                if( wall_configuration.contact_condition == "Axisymmetric" ):
                    self.rigid_wall_bbox[sizei].SetAxisymmetric()

                # rigid wall contact search process
                self.contact_search_process  = []
                self.contact_search_process.append( RigidWallContactSearch( self.rigid_wall_bbox[sizei], self.model_part) )



    #######################################################################
    def ExecuteContactSearch(self):
        if( self.rigid_wall_active ):
            # set properties for rigid wall conditions
            self.model_part.Properties[1].SetValue(PENALTY_PARAMETER,self.penalty_parameter)
            for sizei in range(0,self.number_of_walls):
                self.contact_search_process[sizei].ExecuteInitializeSolutionStep()
        
    #######################################################################
    def UpdatePosition(self):
        if( self.rigid_wall_active ):
            for sizei in range(0,self.number_of_walls):
                self.contact_search_process[sizei].ExecuteFinalizeSolutionStep()
                

    #######################################################################   
    def RigidWallActive(self):
        return self.rigid_wall_active  

 
