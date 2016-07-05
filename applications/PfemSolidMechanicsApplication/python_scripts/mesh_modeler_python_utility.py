from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PfemSolidMechanicsApplication as KratosPfemSolid
KratosMultiphysics.CheckForPreviousImport()


class ModelerUtility:
    #

    def __init__(self, model_part, domain_size, remesh_domains, contact_search, wall_contact_search):

        self.echo_level = 1        
        self.model_part  = model_part
        self.domain_size = domain_size

        # set remesh flags
        self.modeler_active = False
        self.remesh_domains = remesh_domains
        self.contact_search = contact_search
        self.wall_contact_search = wall_contact_search
        self.neighbours_set = False

        # mesh modeler vector
        self.counter = 1
        self.meshing_domains = [] 
        self.remesh_frequencies = []

        # contact modeler parameters
        self.mu_static  = 0.3
        self.mu_dynamic = 0.2

        self.reference_element = "Element2D3N"
        
        self.initial_transfer      = True
        self.contact_alpha_shape   = 1.4
        self.contact_constrained   = False
        self.penalty_contact       = False
        self.friction_active       = False
        self.penalty_parameter     = 1.0
        self.stability_parameter   = 1.0
        self.contact_offset_factor = 0.0

        # time step meshing control parameters

        self.remesh_executed = False
        self.contact_transfer_done = False             

    #
    def GetMeshingStep(self):
        return self.counter
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

            if(self.remesh_domains or self.contact_search or self.wall_contact_search):
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
        nodal_neighbour_search = KratosPfemBase.NodalNeighboursSearch(self.model_part, self.echo_level, number_of_avg_elems, number_of_avg_nodes, mesh_id)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Nodal Search executed ")

    #
    def SearchElementNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
        mesh_id = 0
         
        # define search utility
        elemental_neighbour_search = KratosPfemBase.ElementalNeighboursSearch(self.model_part, self.domain_size, self.echo_level, number_of_avg_elems, mesh_id)

        # execute search:
        elemental_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Elemental Search executed ")

    #
    def ComputeBoundaryNormals(self):

        # define calculation utility
        normals_calculation = KratosPfemBase.BoundaryNormalsCalculation()

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
        skin_build = KratosPfemBase.BuildMeshBoundary(self.model_part, self.domain_size, self.echo_level, mesh_id)

        # execute building:
        skin_build.Execute()

        print("::[Modeler_Utility]:: Mesh Boundary Build executed ")

    #
    def SearchNodalH(self):

        if(self.neighbours_set):
            # define search utility
            nodal_h_search = KratosMultiphysics.FindNodalHProcess(self.model_part)
            # execute search:
            nodal_h_search.Execute()

            # for node in self.model_part.Nodes:
                # nodal_h  = node.GetSolutionStepValue(NODAL_H);
                # print "nodal_h:",nodal_h

            print("::[Modeler_Utility]:: Nodal H Search executed ")

    #
    def BuildMeshModelers(self, meshing_domains ):

        if(self.remesh_domains):
            self.modeler_active = True

        # set mesing domains
        self.meshing_domains = meshing_domains

        # set modeler utilities
        self.modeler_utils = KratosPfemBase.ModelerUtilities()

        # set transfer utilities
        self.transfer_utils = KratosPfemBase.MeshDataTransferUtilities()
                
        # set the domain labels to mesh modeler
        self.modeler_utils.SetDomainLabels(self.model_part)

        # set remesh frequency vector
        #for domain in self.meshing_domains:
        #    self.remesh_frequencies.append(domain.GetMeshingFrequency())
        
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
            self.contact_modeler =  KratosPfemSolid.ContactDomain2DModeler()
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
        self.ContactInfoParameters =  KratosPfemBase.InfoParameters()
 
        # set refine parameters to mesh modeler
        self.ContactRefiningParameters =  KratosPfemBase.RefiningParameters()
            
        self.ContactRefiningParameters.Initialize()
        self.ContactRefiningParameters.SetAlphaParameter(self.contact_alpha_shape)

        # set transfer parameters
        self.ContactTransferParameters =  KratosPfemBase.TransferParameters()
        cauchy_stress = "CAUCHY_STRESS_VECTOR"
        deformation_gradient = "DEFORMATION_GRADIENT"
        self.ContactTransferParameters.SetVariable(KratosMultiphysics.KratosGlobals.GetVariable(cauchy_stress))
        self.ContactTransferParameters.SetVariable(KratosMultiphysics.KratosGlobals.GetVariable(deformation_gradient))

        # set meshing parameters to mesh modeler
        self.ContactMeshingParameters =  KratosPfemBase.MeshingParameters()
            
        self.ContactMeshingParameters.Initialize()
        
        self.ContactMeshingParameters.SetInfoParameters(self.ContactInfoParameters)
        self.ContactMeshingParameters.SetRefiningParameters(self.ContactRefiningParameters)
        self.ContactMeshingParameters.SetTransferParameters(self.ContactTransferParameters)
                        
        contact_meshing_options = KratosMultiphysics.Flags()

        contact_meshing_options.Set(KratosPfemBase.ModelerUtilities.REMESH, True)
        contact_meshing_options.Set(KratosPfemBase.ModelerUtilities.CONSTRAINED, contact_config.constrained_contact)
            
        self.ContactMeshingParameters.SetOptions(contact_meshing_options)

        self.ContactMeshingParameters.SetOffsetFactor(self.contact_offset_factor)
        self.ContactMeshingParameters.SetAlphaParameter(self.contact_alpha_shape)
                

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
            ContactTransferParameters = KratosPfemBase.TransferParameters()
            cauchy_stress = "CAUCHY_STRESS_VECTOR"
            deformation_gradient = "DEFORMATION_GRADIENT"
            ContactTransferParameters.SetVariable(KratosMultiphysics.KratosGlobals.GetVariable(cauchy_stress))
            ContactTransferParameters.SetVariable(KratosMultiphysics.KratosGlobals.GetVariable(deformation_gradient))

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

            if( self.echo_level > 0 ):
                print("::[Modeler_Utility]:: MESH DOMAIN...", self.counter)

            meshing_options = KratosMultiphysics.Flags()
            self.model_meshing =  KratosPfemBase.ModelMeshing(self.model_part, meshing_options, self.echo_level)

            self.model_meshing.ExecuteInitialize()

            id = 0
            for domain in self.meshing_domains:

                domain.ExecuteMeshing();

                self.remesh_executed = True

                id+=1

            self.model_meshing.ExecuteFinalize()

            self.counter += 1 


    #
