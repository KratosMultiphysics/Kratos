from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemSolidMechanicsApplication import *
CheckForPreviousImport()


class ModelerUtility:
    #

    def __init__(self, model_part, domain_size, remesh_domains, contact_search, rigid_wall_contact_search):

        self.echo_level = 1        
        self.model_part = model_part
        self.domain_size = domain_size

        # set remesh flags
        self.modeler_active = False
        self.remesh_domains = remesh_domains
        self.contact_search = contact_search
        self.rigid_wall_contact_search = rigid_wall_contact_search
        self.neighbours_set = False

        # set mesh modeler
        self.counter = 1
        if(domain_size >= 2):
            self.mesh_modeler = TriangularMesh2DModeler()
        # else:
            # self.mesh_modeler = TetrahedronMesh3DModeler()

        # set contact modeler
        if(domain_size >= 2):
            self.contact_modeler = ContactDomain2DModeler()
        # else:
            # self.contact_modeler = ContactDomain3DModeler()
            

        # mesh modeler parameters
        self.alpha_shape        = 2.4
        self.h_factor           = 0.5
        self.avoid_tip_elements = False
        self.offset_factor      = 0
        self.remesh_frequencies = []

        # contact modeler parameters
        self.mu_static  = 0.3
        self.mu_dynamic = 0.2

        self.initial_transfer      = True
        self.contact_alpha_shape   = 1.4
        self.contact_constrained   = False
        self.penalty_contact       = False
        self.friction_active       = False
        self.penalty_parameter     = 1
        self.stability_parameter   = 1
        self.contact_offset_factor = 0

        self.contact_condition = "ContactDomainLM2DCondition"

        # time step meshing control parameters

        self.remesh_executed = False
        self.contact_transfer_done = False             

    #
    def Initialize(self):

        self.remesh_executed = False
        self.contact_transfer_done = False
        
    #
    def InitializeDomains(self, ReloadFile = False):

        if( self.modeler_active ):        
            print("::[Modeler_Utility]:: Initialize Domains ")
            
            # set active search
            self.search_active = False

            if(self.remesh_domains or self.contact_search or self.rigid_wall_contact_search):
                self.search_active = True

            self.neighbours_set = False
            if(self.search_active):
                # find neighbours
                self.SearchNeighbours()
                # find skin and boundary normals
                if(ReloadFile == False):
                    self.BuildBoundarySkin()
                self.neighbours_set = True

    #
    def SearchNeighbours(self):

        self.SearchNodeNeighbours()
        self.SearchElementNeighbours()

    #
    def SearchNodeNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        mesh_id = 0

        # define search utility
        nodal_neighbour_search = NodalNeighboursSearch(self.model_part, self.echo_level, number_of_avg_elems, number_of_avg_nodes, mesh_id)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Nodal Search executed ")

    #
    def SearchElementNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
        mesh_id = 0
         
        # define search utility
        elemental_neighbour_search = ElementalNeighboursSearch(self.model_part, self.domain_size, self.echo_level, number_of_avg_elems, mesh_id)

        # execute search:
        elemental_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Elemental Search executed ")

    #
    def ComputeBoundaryNormals(self):

        # define calculation utility
        normals_calculation = BoundaryNormalsCalculation()

        # execute calculation:
        #(scaled normals)
        normals_calculation.CalculateBoundaryNormals(self.model_part, self.echo_level)
        #(unit normals)
        # normals_calculation.CalculateBoundaryUnitNormals(model_part, self.echo_level)

        print("::[Modeler_Utility]:: Boundary Normals computed ")

    #
    def BuildBoundarySkin(self):

        print("::[Modeler_Utility]:: Build Boundary Skin ")
        # set building options:
        mesh_id = 0

        # define building utility
        skin_build = BuildBoundarySkin(self.model_part, self.domain_size, self.echo_level, mesh_id)

        # execute building:
        skin_build.Execute()

        print("::[Modeler_Utility]:: Boundary Skin Build executed ")

    #
    def SearchNodalH(self):

        if(self.neighbours_set):
            # define search utility
            nodal_h_search = FindNodalHProcess(self.model_part)
            # execute search:
            nodal_h_search.Execute()

            # for node in self.model_part.Nodes:
                # nodal_h  = node.GetSolutionStepValue(NODAL_H);
                # print "nodal_h:",nodal_h

            print("::[Modeler_Utility]:: Nodal H Search executed ")

    #
    def BuildMeshModeler(self, configuration):

        # definition of the echo level
        if(hasattr(configuration, "echo_level")):
            self.echo_level = configuration.echo_level

        self.mesh_modeler.SetEchoLevel(self.echo_level)
            
         # check domain consistency
        if(configuration.number_domains != len(configuration.mesh_conditions)):
            print("::[Modeler_Utility]:: Number of Domain Meshing Conditions do not match ")

        # check mesh consistency
        # if(configuration.number_domains != self.model_part.NumberOfMeshes):
            # print("::[Modeler_Utility]:: Number of Domain Meshing Conditions and Meshes in model_part do not match " )

        # set the domains number to mesh modeler   
        number_of_domains = self.model_part.NumberOfMeshes();
        if( number_of_domains > configuration.number_domains ): # mesh and constraint meshes
            number_of_domains = configuration.number_domains
            if( number_of_domains == 1 ): # mesh 0 and mesh 1
                number_of_domains += 1
            

        self.mesh_modeler.SetInitialMeshData(number_of_domains)

        # set modeler utilities
        self.modeler_utils = ModelerUtilities()

        # set the domain labels to mesh modeler
        self.modeler_utils.SetDomainLabels(self.model_part)

        mesh_id = 0
        for conditions in configuration.mesh_conditions:

            mesh_id = int(conditions["Subdomain"])
            
            if(conditions["StructuralType"] == "Solid"):
                self.modeler_active = True

            # set remesh-refine conditions to mesh modeler
            critical_mesh_size = conditions["CriticalMeshSize"]

            # set mesh refinement based on wall tip discretization size
            if(conditions["TipRadiusRefine"]):
                # tip arch opening (in degrees = 5-7.5-10)
                tool_arch_opening = 12
                # tip surface length
                tool_arch_length = tool_arch_opening * (3.1416 / 180.0)
                # critical mesh size based on wall tip
                critical_mesh_size = tool_arch_length * conditions["CriticalTipRadius"]

            critical_mesh_size = critical_mesh_size * configuration.size_scale

            print("::[Modeler_Utility]:: Domain Mesh [",conditions["Subdomain"],"] [ Remesh:",conditions["Remesh"],"] [ Refine:",conditions["Refine"],"]" )
            if( conditions["Remesh"] ):
                print("(Type:",conditions["MeshElement"],")")

            #remesh data
            self.mesh_modeler.SetRemeshData(conditions["MeshElement"], "CompositeCondition2D2N", conditions["Remesh"], conditions["Constrained"], conditions["MeshSmoothing"], conditions["JacobiSmoothing"], self.avoid_tip_elements, self.alpha_shape, self.offset_factor, mesh_id)

            #refine data
            self.mesh_modeler.SetRefineData(conditions["Refine"], self.h_factor, critical_mesh_size, conditions["DissipationVariable"], conditions["CriticalDissipation"], conditions["ErrorVariable"], conditions["CriticalError"], mesh_id)

            box_refinement_only = conditions["RefineOnBoxOnly"]

            if(box_refinement_only):

                radius_box = conditions["BoxRadius"] * configuration.size_scale
                center_box = Vector(self.domain_size)
                velocity_box = Vector(self.domain_size)

                for size in range(0, self.domain_size):
                    center_box[size] = conditions["BoxCenter"][size] * configuration.size_scale
                    velocity_box[size] = conditions["BoxVelocity"][size] * configuration.size_scale

                self.mesh_modeler.SetMeshRefiningBox(radius_box, center_box, velocity_box, mesh_id)

            # set remesh frequency
            self.remesh_frequencies.append(conditions["RemeshFrequency"])

    #
    def GetRemeshFrequency(self):
        
        remesh_frequency = 0
        for size in range(0,len(self.remesh_frequencies)):
            if((remesh_frequency > self.remesh_frequencies[size]) or remesh_frequency == 0):
                remesh_frequency = self.remesh_frequencies[size]
        
        return remesh_frequency
            

    #
    def BuildContactModeler(self, contact_config):

        # if restart file is not loaded geometric searches are needed previously
        # find neighbours,find model skin, find nodal_h

        # definition of the echo level
        if(hasattr(contact_config, "echo_level")):
            self.echo_level = contact_config.echo_level

        self.contact_modeler.SetEchoLevel(self.echo_level)

        self.contact_condition = contact_config.contact_condition

        self.constrained_contact = contact_config.constrained_contact

        self.friction_active = contact_config.friction_active

        self.contact_offset_factor = contact_config.offset_factor
        self.penalty_parameter = contact_config.penalty_parameter
        self.stability_parameter = contact_config.stability_parameter
        self.mu_static = contact_config.mu_static
        self.mu_dynamic = contact_config.mu_dynamic

        self.initial_transfer = True

        # set contact search frequency
        self.contact_search_frequency = contact_config.contact_search_frequency

    #
    def InitialContactSearch(self):

        if(self.contact_search):
            print("::[Modeler_Utility]:: CONTACT SEARCH START: ", self.contact_condition)
            self.ContactTransfer()
            self.contact_modeler.GenerateContactMesh(self.model_part, "Element2D", self.contact_condition, self.constrained_contact, self.alpha_shape, self.h_factor, self.contact_offset_factor, self.penalty_parameter, self.stability_parameter, self.friction_active, self.mu_static, self.mu_dynamic);

    #
    def InitializeStep(self):

        self.remesh_executed = False
        self.contact_transfer_done = False

        if(self.initial_transfer):
            self.initial_transfer = False

    #
    def ContactTransfer(self):

        if(self.contact_transfer_done == False):
            print("::[Modeler_Utility]:: TRANSFER CONTACTS")
            self.contact_modeler.TransferContactData(self.model_part, self.initial_transfer);
            self.contact_transfer_done = True

    #
    def ContactSearch(self):

        if(self.contact_search):
            print("::[Modeler_Utility]:: CONTACT SEARCH : ", self.contact_condition)

            self.ContactTransfer()
            self.contact_modeler.GenerateContactMesh(self.model_part, "Element2D", self.contact_condition, self.constrained_contact, self.alpha_shape, self.h_factor, self.contact_offset_factor, self.penalty_parameter, self.stability_parameter, self.friction_active, self.mu_static, self.mu_dynamic);
            
    #
    def RemeshDomains(self):

        if(self.remesh_domains):
            if(self.contact_search):
                self.ContactTransfer()

            if( self.echo_level > 0 ):
                print("::[Modeler_Utility]:: MESH DOMAIN...", self.counter)

            self.mesh_modeler.GenerateMesh(self.model_part);
            self.remesh_executed = True
            self.counter += 1 

    #
    def SetRigidWall(self, rigid_wall):

        if( rigid_wall.RigidWallActive() ):
            rigid_wall_bbox = rigid_wall.RigidWallBoundingBoxes()
            for sizei in range(0, len(rigid_wall_bbox)):
                self.mesh_modeler.SetRigidWall( rigid_wall_bbox[sizei] )

    #
