from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()


class ModelerUtility:
    #

    def __init__(self, model_part, domain_size, remesh_domains, contact_search):

        self.model_part = model_part
        self.domain_size = domain_size

        # set remesh flags
        self.remesh_domains = False
        if(remesh_domains == "True"):
            self.remesh_domains = True

        self.contact_search = False
        if(contact_search == "True"):
            self.contact_search = True

        self.neighbours_set = False

        # set mesh modeler
        if(domain_size == 2):
            self.mesh_modeler = TriangleMesh2DModeler()
        # else:
            # self.mesh_modeler = TetrahedronMesh3DModeler()

        # set contact modeler
        if(domain_size == 2):
            self.contact_modeler = ContactDomain2DModeler()
        # else:
            # self.contact_modeler = ContactDomain3DModeler()

        # mesh modeler parameters
        self.alpha_shape = 2.4
        self.h_factor = 0.5

        self.remesh = False
        self.refine = False
        self.constrained = False
        self.laplacian_smoothing = False
        self.jacobi_smoothing = False
        self.avoid_tip_elements = False
        self.offset_factor = 0

        # contact modeler parameters
        self.mu_static = 0.3
        self.mu_dynamic = 0.2

        self.initial_transfer = True
        self.contact_alpha_shape = 1.4
        self.contact_constrained = False
        self.penalty_contact = False
        self.friction_active = False
        self.penalty_factor = 1
        self.contact_offset_factor = 0

        self.contact_condition = "ContactDomain2DCondition"

        # time step meshing control parameters
        self.remesh_frequency = 1
        self.contact_search_frequency = 1

        self.remesh_executed = False
        self.remesh_step = 0
        self.contact_transfer_done = False
        self.contact_search_step = 0

    #
    def Initialize(self, remesh_step, contact_search_step):

        self.remesh_step = remesh_step
        self.contact_search_step = contact_search_step

    #
    def InitializeDomains(self):

        # set active search
        self.search_active = False

        if(self.remesh_domains or self.contact_search):
            self.search_active = True

        self.neighbours_set = False
        if(self.search_active):
            # find neighbours
            self.SearchNeighbours()
            # find skin and boundary normals
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
        method = 0

        # define search utility
        nodal_neighbour_search = NodalNeighboursSearch(self.model_part, number_of_avg_elems, number_of_avg_nodes, method)

        # execute search:
        nodal_neighbour_search.Execute()

        # print " Nodal Search executed "

    #
    def SearchElementNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
        method = 0

        # define search utility
        elemental_neighbour_search = ElementalNeighboursSearch(self.model_part, self.domain_size, number_of_avg_elems, method)

        # execute search:
        elemental_neighbour_search.Execute()

        # print " Element Search executed "

    #
    def ComputeBoundaryNormals(self):

        # define calculation utility
        normals_calculation = BoundaryNormalsCalculation()

        # execute calculation:
        normals_calculation.CalculateBoundaryNormals(self.model_part, self.domain_size)
        # normals_calculation.CalculateBoundaryUnitNormals(model_part)

    #
    def BuildBoundarySkin(self):

        # set building options:
        preserve = 1

        # define building utility
        skin_build = BuildBoundarySkin(self.model_part, self.domain_size, preserve)

        # execute building:
        skin_build.Execute()

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

            # print " Nodal H search executed "

    #
    def BuildMeshModeler(self, configuration):

         # check domain consistency
        if(configuration.number_domains != len(configuration.mesh_conditions)):
            print(" Number of Domain Meshing Conditions do not match ")

        # check mesh consistency
        # if(configuration.number_domains != self.model_part.NumberOfMeshes):
            # print " Number of Domain Meshing Conditions and Meshes in model_part do not match "

        # set the domains number to mesh modeler
        self.mesh_modeler.SetInitialMeshData(configuration.number_domains)

        # set the domain labels to mesh modeler
        self.modeler_utils = ModelerUtilities()
        self.modeler_utils.SetDomainLabels(self.model_part)

        # set remesh-refine conditions to mesh modeler
        radius_critical = configuration.critical_radius * configuration.size_scale

        if hasattr(configuration, 'box_refinement_only'):
            box_refinement_only = configuration.box_refinement_only

        for conditions in configuration.mesh_conditions:

            if(conditions["Remesh"] == 1):
                self.remesh = True
            else:
                self.remesh = False

            if(conditions["Refine"] == 1):
                self.refine = True
            else:
                self.refine = False

            if(conditions["Constrained"] == 1):
                self.constrained = True
            else:
                self.constrained = False

            if(conditions["MeshSmoothing"] == 1):
                self.laplacian_smoothing = True
            else:
                self.laplacian_smoothing = False

            if(conditions["JacobiSmoothing"] == 1):
                self.jacobi_smoothing = True
            else:
                self.jacobi_smoothing = False

            domain = int(conditions["Subdomain"])

            print("SET MESH DOMAIN DATA")
            print(" --> Domain [", conditions["Subdomain"], "] ", conditions["MeshElement"], " remesh: ", conditions["Remesh"])
            # self.mesh_modeler.SetRemeshData(conditions["MeshElement"],"Condition2D",,self.remesh,,self.constrained,,self.laplacian_smoothing,,self.jacobi_smoothing,,self.avoid_tip_elements,,self.alpha_shape,domain);
            self.mesh_modeler.SetRemeshData(conditions["MeshElement"], "SkinCondition2D", self.remesh, self.constrained, self.laplacian_smoothing, self.jacobi_smoothing, self.avoid_tip_elements, self.alpha_shape, self.offset_factor, domain)
            self.mesh_modeler.SetRefineData(self.refine, self.h_factor, configuration.critical_dissipation, radius_critical, configuration.reference_error, domain)

            if(box_refinement_only == "True"):

                center_box = Vector(self.domain_size)
                velocity_box = Vector(self.domain_size)
                for size in range(0, self.domain_size):
                    center_box[size] = configuration.box_center[size] * configuration.size_scale
                    velocity_box[size] = configuration.box_velocity[size] * configuration.size_scale

                radius_box = configuration.box_radius * configuration.size_scale

                self.mesh_modeler.SetRefiningBox(radius_box, center_box, velocity_box)

        # set remesh frequency
        self.remesh_frequency = configuration.remesh_frequency

    #
    def BuildContactModeler(self, contact_config):

        # if restart file is not loaded geometric searches are needed previously
        # find neighbours,find model skin, find nodal_h

        self.contact_condition = contact_config.contact_condition

        if(contact_config.constrained_contact == "True"):
            self.constrained_contact = True
        else:
            self.constrained_contact = False

        if(contact_config.friction_active == "True"):
            self.friction_active = True
        else:
            self.friction_active = False

        if(contact_config.penalty_contact == "True"):
            self.penalty_contact = True
        else:
            self.penalty_contact = False

        self.contact_offset_factor = contact_config.offset_factor
        self.penalty_factor = contact_config.penalty_factor
        self.mu_static = contact_config.mu_static
        self.mu_dynamic = contact_config.mu_dynamic

        self.initial_transfer = True

        # set contact search frequency
        self.contact_search_frequency = contact_config.contact_search_frequency

    #
    def InitialContactSearch(self):

        if(self.contact_search):
            print(" CONTACT SEARCH START: ", self.contact_condition)
            self.ContactTransfer()
            self.contact_modeler.GenerateContactMesh(self.model_part, "Element2D", self.contact_condition, self.constrained_contact, self.alpha_shape, self.h_factor, self.contact_offset_factor, self.penalty_factor, self.friction_active, self.mu_static, self.mu_dynamic, self.penalty_contact)

    #
    def InitializeStep(self):

        self.remesh_executed = False
        self.contact_transfer_done = False

        if(self.initial_transfer):
            self.initial_transfer = False

    #
    def ContactTransfer(self):

        if(self.contact_transfer_done == False):
            print("TRANSFER CONTACTS")
            self.contact_modeler.TransferContactData(self.model_part, self.initial_transfer)
            self.contact_transfer_done = True

    #
    def ContactSearch(self, current_step):

        if(self.contact_search):
            print(" CONTACT SEARCH : ", self.contact_condition)

            if(self.remesh_executed or current_step == self.contact_search_step):
                self.ContactTransfer()
                self.contact_modeler.GenerateContactMesh(self.model_part, "Element2D", self.contact_condition, self.constrained_contact, self.alpha_shape, self.h_factor, self.contact_offset_factor, self.penalty_factor, self.friction_active, self.mu_static, self.mu_dynamic, self.penalty_contact);

                if(current_step == self.contact_search_step):
                    self.contact_search_step = self.contact_search_step + self.contact_search_frequency

    #
    def RemeshDomains(self, current_step):

        if(self.remesh_domains):
            if(current_step == self.remesh_step):
                if(self.contact_search):
                    self.ContactTransfer()

                print("MESH DOMAIN")
                self.mesh_modeler.GenerateMesh(self.model_part);
                self.remesh_executed = True

                self.remesh_step = self.remesh_step + self.remesh_frequency

    #
    def RemeshDomains(self, current_step, rigid_wall):

        if(self.remesh_domains):
            if(current_step == self.remesh_step):
                if(self.contact_search):
                    self.ContactTransfer()

                self.SetRigidWall(rigid_wall);
                print("MESH DOMAIN")
                self.mesh_modeler.GenerateMesh(self.model_part);
                self.remesh_executed = True

                self.remesh_step = self.remesh_step + self.remesh_frequency

    #
    def SetRigidWall(self, rigid_wall):

        rigid_wall_active = rigid_wall.RigidWallActive()

        # if rigid walls are present
        if(rigid_wall_active):
            center = rigid_wall.RigidWallCenter()
            tip_radius = rigid_wall.RigidWallTipRadius()
            self.mesh_modeler.SetToolTip(tip_radius, center);

    #
