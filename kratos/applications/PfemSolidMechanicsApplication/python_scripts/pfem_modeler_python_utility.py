from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PfemBaseApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
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

        # mesh modeler vector
        self.counter = 1
        self.mesh_ids = []
        self.mesh_modelers = []

        # mesh modeler parameters
        self.alpha_shape        = 2.4
        self.h_factor           = 0.5
        self.avoid_tip_elements = False
        self.offset_factor      = 0.0
        self.remesh_frequencies = []

        # contact modeler parameters
        self.mu_static  = 0.3
        self.mu_dynamic = 0.2

        self.initial_transfer      = True
        self.contact_alpha_shape   = 1.4
        self.contact_constrained   = False
        self.penalty_contact       = False
        self.friction_active       = False
        self.penalty_parameter     = 1.0
        self.stability_parameter   = 1.0
        self.contact_offset_factor = 0.0

        self.reference_element = "Element2D"
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
                    self.BuildMeshBoundary()
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
    def BuildMeshBoundary(self):

        print("::[Modeler_Utility]:: Build Mesh Boundary ")
        # set building options:
        mesh_id = 0

        # define building utility
        skin_build = BuildMeshBoundary(self.model_part, self.domain_size, self.echo_level, mesh_id)

        # execute building:
        skin_build.Execute()

        print("::[Modeler_Utility]:: Mesh Boundary Build executed ")

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
    def BuildMeshModelers(self, configuration, rigid_wall):

        # definition of the echo level
        if(hasattr(configuration, "echo_level")):
            self.echo_level = configuration.echo_level

            
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
                            
        # set mesh modeler array
        self.mesh_modeler_vector = []

        # set modeler utilities
        self.modeler_utils = ModelerUtilities()

        # set transfer utilities
        self.transfer_utils = MeshDataTransferUtilities()
                
        # set the domain labels to mesh modeler
        self.modeler_utils.SetDomainLabels(self.model_part)

        mesh_id = 0

        for parameters in configuration.mesh_conditions:

            # set mesh modeler
            mesh_modeler = TriangularMesh2DModeler()
            #if(self.domain_size == 3):
            #    mesh_modeler = TetrahedronMesh3DModeler()
                
            
            mesh_modeler.Initialize()

            mesh_modeler.SetEchoLevel(self.echo_level)

            mesh_id = int(parameters["Subdomain"])
            
            if(parameters["StructuralType"] == "Solid"):
                self.modeler_active = True

            # create info parameters
            self.InfoParameters   = InfoParameters()
 
            # set refine parameters to mesh modeler
            self.RefiningParameters = RefiningParameters()
            
            self.RefiningParameters.Initialize()

            refining_options = Flags()
            refining_options.Set(ModelerUtilities.REFINE, parameters["Refine"])
            refining_options.Set(ModelerUtilities.REFINE_ADD_NODES, True)
            #refine elements
            refining_options.Set(ModelerUtilities.REFINE_ELEMENTS, True)
            refining_options.Set(ModelerUtilities.REFINE_ELEMENTS_ON_DISTANCE, True)
            refining_options.Set(ModelerUtilities.REFINE_ELEMENTS_ON_THRESHOLD, True)  
            #refine boundary
            refining_options.Set(ModelerUtilities.REFINE_BOUNDARY, True)
            #refining_options.Set(ModelerUtilities.REFINE_BOUNDARY_ON_DISTANCE, True)
            refining_options.Set(ModelerUtilities.REFINE_BOUNDARY_ON_THRESHOLD, True)  

            self.RefiningParameters.SetRefiningOptions(refining_options)
            
            removing_options = Flags()
            if( parameters["Refine"] ):
                removing_options.Set(ModelerUtilities.REMOVE_NODES, True)
                removing_options.Set(ModelerUtilities.REMOVE_NODES_ON_DISTANCE, True)
                removing_options.Set(ModelerUtilities.REMOVE_NODES_ON_ERROR, True)
                removing_options.Set(ModelerUtilities.REMOVE_NODES_ON_THRESHOLD, True)
                
                #removing_options.Set(ModelerUtilities.REMOVE_BOUNDARY_NODES, True)
                #removing_options.Set(ModelerUtilities.REMOVE_BOUNDARY_NODES_ON_DISTANCE, True)
                #removing_options.Set(ModelerUtilities.REMOVE_BOUNDARY_NODES_ON_ERROR, True)

            self.RefiningParameters.SetRemovingOptions(removing_options)

            self.RefiningParameters.SetAlphaParameter(self.alpha_shape)

            critical_mesh_size = parameters["CriticalMeshSize"]

            # set mesh refinement based on wall tip discretization size
            if(parameters["TipRadiusRefine"]):
                # tip arch opening (in degrees = 5-7.5-10)
                tool_arch_opening = 12
                # tip surface length
                tool_arch_length = tool_arch_opening * (3.1416 / 180.0)
                # critical mesh size based on wall tip
                critical_mesh_size = tool_arch_length * parameters["CriticalTipRadius"]

            critical_mesh_size = critical_mesh_size * configuration.size_scale
            critical_mesh_side = critical_mesh_size * 3

            self.RefiningParameters.SetCriticalRadius(critical_mesh_size)                       
            self.RefiningParameters.SetCriticalSide(critical_mesh_side)

            # set mesh refinement in box
            box_refinement_only = parameters["RefineOnBoxOnly"]

            if(box_refinement_only):

                radius_box = parameters["BoxRadius"] * configuration.size_scale
                center_box = Vector(self.domain_size)
                velocity_box = Vector(self.domain_size)

                for size in range(0, self.domain_size):
                    center_box[size] = parameters["BoxCenter"][size] * configuration.size_scale
                    velocity_box[size] = parameters["BoxVelocity"][size] * configuration.size_scale

                refining_box = SpatialBoundingBox(center_box, radius_box, velocity_box)
                
                self.RefiningParameters.SetRefiningBox(refining_box)
                

            self.RefiningParameters.SetThresholdVariable(globals()[parameters["DissipationVariable"]])
            self.RefiningParameters.SetReferenceThreshold(parameters["CriticalDissipation"])

            self.RefiningParameters.SetErrorVariable(globals()[parameters["ErrorVariable"]])
            self.RefiningParameters.SetReferenceError(parameters["CriticalError"])


            # set transfer parameters
            self.TransferParameters = TransferParameters()
            cauchy_stress = "CAUCHY_STRESS_VECTOR"
            deformation_gradient = "DEFORMATION_GRADIENT"
            self.TransferParameters.SetVariable(globals()[cauchy_stress])
            self.TransferParameters.SetVariable(globals()[deformation_gradient])

            # set meshing parameters to mesh modeler
            self.MeshingParameters = MeshingParameters()
            
            self.MeshingParameters.Initialize()

            self.MeshingParameters.SetInfoParameters(self.InfoParameters)
            self.MeshingParameters.SetRefiningParameters(self.RefiningParameters)
            self.MeshingParameters.SetTransferParameters(self.TransferParameters)
                        
            meshing_options = Flags()

            meshing_options.Set(ModelerUtilities.REMESH, parameters["Remesh"])
            meshing_options.Set(ModelerUtilities.CONSTRAINED, parameters["Constrained"])
            meshing_options.Set(ModelerUtilities.REFINE, parameters["Refine"])
            meshing_options.Set(ModelerUtilities.MESH_SMOOTHING, parameters["MeshSmoothing"])
            meshing_options.Set(ModelerUtilities.VARIABLES_SMOOTHING, parameters["JacobiSmoothing"])
            
            self.MeshingParameters.SetOptions(meshing_options)

            self.MeshingParameters.SetOffsetFactor(self.offset_factor)
            self.MeshingParameters.SetAlphaParameter(self.alpha_shape)
                
            self.MeshingParameters.SetReferenceElement(parameters["MeshElement"])
            self.MeshingParameters.SetReferenceCondition("CompositeCondition2D2N")
                

            #Pre Meshing Processes

            # A:
            refine_mesh_elements = RefineMeshElements(self.model_part, self.RefiningParameters, mesh_id, self.echo_level)
            mesh_modeler.SetPreMeshingProcess(refine_mesh_elements)

            # B:
            #refine_mesh_boundary = RefineMeshBoundary(self.model_part, self.RefiningParameters, self.InfoParameters, mesh_id, self.echo_level)
            
            rigid_walls_container = BoundingBoxContainer()
            if( rigid_wall.RigidWallActive() ):
                rigid_wall_bbox = rigid_wall.RigidWallBoundingBoxes()
                for sizei in range(0, len(rigid_wall_bbox)):
                    rigid_walls_container.PushBack( rigid_wall_bbox[sizei] )

            
            refine_mesh_boundary = ContactRefineMeshBoundary(self.model_part, rigid_walls_container, self.MeshingParameters, mesh_id, self.echo_level)

            mesh_modeler.SetPreMeshingProcess(refine_mesh_boundary)

            # C:
            remove_mesh_nodes = RemoveMeshNodes(self.model_part, self.MeshingParameters, mesh_id, self.echo_level)
            mesh_modeler.SetPreMeshingProcess(remove_mesh_nodes)


            #Post Meshing Processes
            rebuild_mesh_boundary = ReconstructMeshBoundary(self.model_part, self.MeshingParameters, mesh_id, self.echo_level)

            
            
            mesh_modeler.SetPostMeshingProcess(rebuild_mesh_boundary)
        
            print("::[Modeler_Utility]:: Domain Mesh [",parameters["Subdomain"],"] [ Remesh:",parameters["Remesh"],"] [ Refine:",parameters["Refine"],"]" )
            if( parameters["Remesh"] ):
                print("(Type:",parameters["MeshElement"],")")

            mesh_modeler.SetMeshingParameters(self.MeshingParameters, mesh_id)

            self.mesh_modelers.append(mesh_modeler)

            self.mesh_ids.append(mesh_id)

            ###############

            # set remesh frequency
            self.remesh_frequencies.append(parameters["RemeshFrequency"])

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

        # set contact modeler
        if(self.domain_size >= 2):
            self.contact_modeler = ContactDomain2DModeler()
            # else:
            # self.contact_modeler = ContactDomain3DModeler()
            

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

        #set contact modeler parameters
        self.contact_modeler.Initialize()

        self.contact_modeler.SetEchoLevel(self.echo_level)
            
        # create info parameters
        self.ContactInfoParameters   = InfoParameters()
 
        # set refine parameters to mesh modeler
        self.ContactRefiningParameters = RefiningParameters()
            
        self.ContactRefiningParameters.Initialize()
        self.ContactRefiningParameters.SetAlphaParameter(self.alpha_shape)

        # set transfer parameters
        self.ContactTransferParameters = TransferParameters()
        cauchy_stress = "CAUCHY_STRESS_VECTOR"
        deformation_gradient = "DEFORMATION_GRADIENT"
        self.ContactTransferParameters.SetVariable(globals()[cauchy_stress])
        self.ContactTransferParameters.SetVariable(globals()[deformation_gradient])

        # set meshing parameters to mesh modeler
        self.ContactMeshingParameters = MeshingParameters()
            
        self.ContactMeshingParameters.Initialize()
        
        self.ContactMeshingParameters.SetInfoParameters(self.ContactInfoParameters)
        self.ContactMeshingParameters.SetRefiningParameters(self.ContactRefiningParameters)
        self.ContactMeshingParameters.SetTransferParameters(self.ContactTransferParameters)
                        
        contact_meshing_options = Flags()

        contact_meshing_options.Set(ModelerUtilities.REMESH, True)
        contact_meshing_options.Set(ModelerUtilities.CONSTRAINED, contact_config.constrained_contact)
            
        self.ContactMeshingParameters.SetOptions(contact_meshing_options)

        self.ContactMeshingParameters.SetOffsetFactor(self.offset_factor)
        self.ContactMeshingParameters.SetAlphaParameter(self.alpha_shape)
                
        self.ContactMeshingParameters.SetReferenceElement(self.reference_element)
        self.ContactMeshingParameters.SetReferenceCondition(self.contact_condition)
                
        mesh_id = 0
        self.contact_modeler.SetMeshingParameters(self.ContactMeshingParameters, mesh_id)

        # set contact search frequency
        self.contact_search_frequency = contact_config.contact_search_frequency

    #
    def InitialContactSearch(self):

        if(self.contact_search):
            print("::[Modeler_Utility]:: CONTACT SEARCH START: ", self.contact_condition)
            self.ContactTransfer()
            self.contact_modeler.GenerateContactMesh(self.model_part, self.penalty_parameter, self.stability_parameter, self.friction_active, self.mu_static, self.mu_dynamic);

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
            # set transfer parameters
            ContactTransferParameters = TransferParameters()
            cauchy_stress = "CAUCHY_STRESS_VECTOR"
            deformation_gradient = "DEFORMATION_GRADIENT"
            ContactTransferParameters.SetVariable(globals()[cauchy_stress])
            ContactTransferParameters.SetVariable(globals()[deformation_gradient])

            self.contact_modeler.TransferContactData(self.model_part, ContactTransferParameters, self.initial_transfer);
            self.contact_transfer_done = True

    #
    def ContactSearch(self):

        if(self.contact_search):
            print("::[Modeler_Utility]:: CONTACT SEARCH : ", self.contact_condition)

            self.ContactTransfer()
            self.contact_modeler.GenerateContactMesh(self.model_part, self.penalty_parameter, self.stability_parameter, self.friction_active, self.mu_static, self.mu_dynamic);
            
    #
    def RemeshDomains(self):

        if(self.remesh_domains):
            if(self.contact_search):
                self.ContactTransfer()

            if( self.echo_level >= 0 ):
                print("::[Modeler_Utility]:: MESH DOMAIN...", self.counter)

            meshing_options = Flags()
            self.model_meshing = ModelMeshing(self.model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()

            id = 0
            for mesher in self.mesh_modelers:

                mesh_id = self.mesh_ids[id]
                
                mesher.InitializeMeshModeler(self.model_part,mesh_id)

                mesher.GenerateMesh(self.model_part,mesh_id);

                mesher.FinalizeMeshModeler(self.model_part,mesh_id)

                self.remesh_executed = True

                id+=1

            self.model_meshing.ExecuteFinalize()

            self.counter += 1 


    #
