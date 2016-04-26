from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *
CheckForPreviousImport()


class RigidWallUtility:
    #

    def __init__(self, model_part, domain_size, wall_configuration):

        self.model_part  = model_part
        self.domain_size = domain_size

        # definition of the echo level
        if(hasattr(wall_configuration, "echo_level")):
            self.echo_level = wall_configuration.echo_level

        # wall parameters
        self.rigid_wall_active = wall_configuration.rigid_wall
        self.number_of_walls   = wall_configuration.number_of_walls

         # geometry parameters
        size_scale = wall_configuration.size_scale
       
        # material parameters
        self.penalty_parameters     = []

        # rigid walls
        self.rigid_wall_bbox        = []

        # contact processes
        self.contact_search_process = []

        if(self.rigid_wall_active):

            sizei = 0

            for conditions in wall_configuration.wall_conditions:

                dimension = 2
                axisymmetric = False
                if(conditions["ContactCondition"] == "Axisymmetric"):
                    axisymmetric = True
                if(conditions["ContactCondition"] == "3D"):
                    dimension = 3
                
                if( conditions["WallType"] == "NOSE-WALL" ): 
                    
                    number_of_noses  = int(conditions["NumberOfNoses"])
                    tip_radius       = Vector(number_of_noses)
                    rake_angles      = Vector(number_of_noses)
                    clearance_angles = Vector(number_of_noses)
                    convexities      = Vector(number_of_noses)
                    tip_centers      = Matrix(number_of_noses, 3)
                
                    wall_velocity         = Vector(3)
                    wall_angular_velocity = Vector(3)
                    wall_reference_point  = Vector(3)

                    wall_noses = conditions["WallNoses"]
                    counter = 0
                    for nose in wall_noses:
                        tip_radius[counter]       = nose["TipRadius"] * size_scale
                        rake_angles[counter]      = nose["RakeAngle"]
                        clearance_angles[counter] = nose["ClearanceAngle"]
                        convexities[counter]      = nose["Convexity"]
                        for size in range(0, 3):
                            tip_centers[counter, size] = nose["TipCenter"][size] * size_scale
                        counter +=1
                    
                    for size in range(0, 3):
                        wall_velocity[size]         = conditions["LinearVelocity"][size]
                        wall_angular_velocity[size] = conditions["AngularVelocity"][size]
                        wall_reference_point[size]  = conditions["RotationCenter"][size]

                    self.rigid_wall_bbox.append(RigidNoseWallBoundingBox(int(conditions["Subdomain"]), convexities, tip_radius, rake_angles, clearance_angles, tip_centers, wall_velocity, wall_angular_velocity, wall_reference_point))

                    self.penalty_parameters.append(conditions["PenaltyParameter"])
                    self.rigid_wall_bbox[sizei].SetDimension(dimension)
                    if(axisymmetric):
                        self.rigid_wall_bbox[sizei].SetAxisymmetric()

                    # rigid wall contact search process
                    self.contact_search_process.append(RigidNoseWallContactSearch(self.rigid_wall_bbox[sizei], self.model_part, self.echo_level))

                elif(conditions["WallType"] == "PLANE"):
                    
                    wall_point  = Vector(3)
                    wall_normal = Vector(3)
                    
                    wall_velocity         = Vector(3)
                    wall_angular_velocity = Vector(3)
                    wall_reference_point  = Vector(3)
                    
                    wall_plane = conditions["WallPlane"]                        
                    convexity = wall_plane["Convexity"]
                    
                    for size in range(0, 3):
                        wall_point[size]            = wall_plane["WallPoint"][size]
                        wall_normal[size]           = wall_plane["WallNormal"][size]
                        wall_velocity[size]         = conditions["LinearVelocity"][size]
                        wall_angular_velocity[size] = conditions["AngularVelocity"][size]
                        wall_reference_point[size]  = conditions["RotationCenter"][size]
                        
                    self.rigid_wall_bbox.append(RigidPlaneWallBoundingBox(int(conditions["Subdomain"]), convexity, wall_point, wall_normal, wall_velocity, wall_angular_velocity, wall_reference_point))

                    self.penalty_parameters.append(conditions["PenaltyParameter"])
                    self.rigid_wall_bbox[sizei].SetDimension(dimension)
                    if(axisymmetric):
                        self.rigid_wall_bbox[sizei].SetAxisymmetric()

                    # rigid wall contact search process
                    self.contact_search_process.append(RigidPlaneWallContactSearch(self.rigid_wall_bbox[sizei], self.model_part, 8))

                        
                elif(conditions["WallType"] == "CIRCLE"):
                    
                    wall_center  = Vector(3)
                    
                    wall_velocity         = Vector(3)
                    wall_angular_velocity = Vector(3)
                    wall_reference_point  = Vector(3)
    
                    wall_circle = conditions["WallCircle"]                        
                    radius = wall_circle["Radius"]
                    convexity = wall_circle["Convexity"]

                    for size in range(0, 3):
                        wall_center[size]           = wall_circle["Center"][size]
                        wall_velocity[size]         = conditions["LinearVelocity"][size]
                        wall_angular_velocity[size] = conditions["AngularVelocity"][size]
                        wall_reference_point[size]  = conditions["RotationCenter"][size]

                    self.rigid_wall_bbox.append(RigidCircleWallBoundingBox(int(conditions["Subdomain"]), convexity, radius, wall_center, wall_velocity, wall_angular_velocity, wall_reference_point))

                    self.penalty_parameters.append(conditions["PenaltyParameter"])
                    self.rigid_wall_bbox[sizei].SetDimension(dimension)
                    if(axisymmetric):
                        self.rigid_wall_bbox[sizei].SetAxisymmetric()

                    # rigid wall contact search process
                    self.contact_search_process.append(RigidCircleWallContactSearch(self.rigid_wall_bbox[sizei], self.model_part, self.echo_level))

                sizei += 1

    #
    def GetPenaltyParameter(self):
        
        penalty_parameter = 0
        for size in range(0,len(self.penalty_parameters)):
            if((penalty_parameter < self.penalty_parameters[size]) or penalty_parameter == 0):
                penalty_parameter = self.penalty_parameters[size]
        
        return penalty_parameter
    #
    def ExecuteContactSearch(self):
        if(self.rigid_wall_active):
            # set properties for rigid wall conditions
            penalty_parameter = self.GetPenaltyParameter()
            size = self.model_part.NumberOfProperties(0)-1 #take in account that penalty will be set to this property:: other option size=0
            self.model_part.Properties[size].SetValue(PENALTY_PARAMETER, penalty_parameter)
            for size in range(0, self.number_of_walls):
                self.contact_search_process[size].ExecuteInitializeSolutionStep()

    #
    def UpdatePosition(self):
        if(self.rigid_wall_active):
            for size in range(0, self.number_of_walls):
                self.contact_search_process[size].ExecuteFinalizeSolutionStep()

    #
    def RigidWallActive(self):
        return self.rigid_wall_active

    #
    def RigidWallBoundingBoxes(self):
        return self.rigid_wall_bbox
