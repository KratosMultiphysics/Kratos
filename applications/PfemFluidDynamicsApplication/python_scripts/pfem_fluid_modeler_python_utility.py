from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay


class ModelerUtility:
    #

    def __init__(self, model_part, dimension, remesh_domains):

        self.echo_level = 1        
        self.model_part = model_part
        self.dimension = dimension

        # set remesh flags
        self.modeler_active = False
        self.remesh_domains = remesh_domains
        self.neighbours_set = False

        # mesh modeler vector
        self.counter = 1
        self.mesh_modelers = []

        # mesh modeler parameters
        self.alpha_shape        = 1.4
        self.h_factor           = 0.5
        self.avoid_tip_elements = False
        self.offset_factor      = 0
        self.remesh_frequencies = []

        # time step meshing control parameters
        self.remesh_executed = False

    #
    def Initialize(self):

        self.remesh_executed = False
        
    #
    def InitializeDomains(self, ReloadFile = False):

        if( self.modeler_active ):        
            print("::[Modeler_Utility]:: Initialize Domains ")
            
            # set active search
            self.search_active = False

            if(self.remesh_domains):
                self.search_active = True

            self.neighbours_set = False
            if(self.search_active):
                # find neighbours
                self.SearchNeighbours()
                # find skin and boundary normals
                if(ReloadFile == False):
                    self.BuildMeshBoundaryForFluids()
                    #self.BuildMeshBoundary()
                self.neighbours_set = True

    #
    def SearchNeighbours(self):

        print("::[Modeler_Utility]:: Search Neighbours ")

        self.SearchNodeNeighbours()
        self.SearchElementNeighbours()

    #
    def SearchNodeNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
        number_of_avg_nodes = 10

        # define search utility
        nodal_neighbour_search = KratosDelaunay.NodalNeighboursSearch(self.model_part, self.echo_level, number_of_avg_elems, number_of_avg_nodes)

        # execute search:
        nodal_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Nodal Search executed ")

    #
    def SearchElementNeighbours(self):

        # set search options:
        number_of_avg_elems = 10
         
        # define search utility
        elemental_neighbour_search = KratosDelaunay.ElementalNeighboursSearch(self.model_part, self.dimension, self.echo_level, number_of_avg_elems)

        # execute search:
        elemental_neighbour_search.Execute()

        print("::[Modeler_Utility]:: Elemental Search executed ")

    #
    def ComputeBoundaryNormals(self):

        # define calculation utility
        # normals_calculation = BoundaryNormalsCalculation()
        normals_calculation = KratosDelaunay.BoundaryNormalsCalculation()

        # execute calculation:
        #(scaled normals)
        normals_calculation.CalculateWeightedBoundaryNormals(self.model_part, self.echo_level)
        #(unit normals)
        # normals_calculation.CalculateUnitBoundaryNormals(model_part, self.echo_level)

        print("::[Modeler_Utility]:: Boundary Normals computed ")

    #
    def BuildMeshBoundary(self):

        print("::[Modeler_Utility]:: Build Mesh Boundary ")
        # set building options:

        # define building utility
        # skin_build = BuildMeshBoundary(self.model_part, self.dimension, self.echo_level)
        skin_build = KratosDelaunay.BuildMeshBoundary(self.model_part, self.echo_level)

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
    def ComputeAverageMeshParameters(self):
     
        for domain in self.meshing_domains:
            if(domain.Active()):
                domain.ComputeAverageMeshParameters()       
#
    def ComputeInitialAverageMeshParameters(self):
     

        for domain in self.meshing_domains:
            if(domain.Active()):
                domain.ComputeInitialAverageMeshParameters()       
#

    def BuildMeshModelers(self, configuration):

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

        ## set transfer utilities
        self.transfer_utils = MeshDataTransferUtilities()
                
        # set the domain labels to mesh modeler
        self.modeler_utils.SetDomainLabels(self.model_part)


        for parameters in configuration.mesh_conditions:

            # set mesh modeler
            # if(self.dimension == 2):
            mesh_modeler = TriangularMesh2DModeler()
            # else:
            # mesh_modeler = TetrahedronMesh3DModeler()
                
            
            mesh_modeler.Initialize()

            mesh_modeler.SetEchoLevel(self.echo_level)

            mesh_id = int(parameters["Subdomain"])
            
            if(parameters["StructuralType"] == "Solid"):
                self.modeler_active = True

            # create info parameters
            self.InfoParameters   = MeshingInfoParameters()
 
            # set refine parameters to mesh modeler
            self.RefiningParameters = RefiningParameters()
            
            self.RefiningParameters.Initialize()

            refining_options = Flags()
            refining_options.Set(ModelerUtilities.REFINE, parameters["Refine"])
            refining_options.Set(ModelerUtilities.REFINE_ADD_NODES, True)
            ##refine elements
            refining_options.Set(ModelerUtilities.REFINE_ELEMENTS, True)
            refining_options.Set(ModelerUtilities.REFINE_ELEMENTS_ON_DISTANCE, True)
            refining_options.Set(ModelerUtilities.REFINE_ELEMENTS_ON_THRESHOLD, True)  
            ##refine boundary
            #refining_options.Set(ModelerUtilities.REFINE_BOUNDARY, True)
            ##refining_options.Set(ModelerUtilities.REFINE_BOUNDARY_ON_DISTANCE, True)
            #refining_options.Set(ModelerUtilities.REFINE_BOUNDARY_ON_THRESHOLD, True)  

            self.RefiningParameters.SetRefiningOptions(refining_options)


            removing_options = Flags()
            if( parameters["Refine"] ):
                removing_options.Set(ModelerUtilities.REMOVE_NODES, True)
                removing_options.Set(ModelerUtilities.REMOVE_NODES_ON_DISTANCE, True)
                removing_options.Set(ModelerUtilities.REMOVE_NODES_ON_ERROR, False)
                removing_options.Set(ModelerUtilities.REMOVE_NODES_ON_THRESHOLD, False)

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
                center_box = Vector(self.dimension)
                velocity_box = Vector(self.dimension)

                for size in range(0, self.dimension):
                    center_box[size] = parameters["BoxCenter"][size] * configuration.size_scale
                    velocity_box[size] = parameters["BoxVelocity"][size] * configuration.size_scale

                refining_box = SpatialBoundingBox(center_box, radius_box, velocity_box)
                
                self.RefiningParameters.SetRefiningBox(refining_box)
                

            self.RefiningParameters.SetThresholdVariable(globals()[parameters["DissipationVariable"]])
            self.RefiningParameters.SetReferenceThreshold(parameters["CriticalDissipation"])

            self.RefiningParameters.SetErrorVariable(globals()[parameters["ErrorVariable"]])
            self.RefiningParameters.SetReferenceError(parameters["CriticalError"])

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
            #meshing_options.Set(ModelerUtilities.MESH_SMOOTHING, parameters["MeshSmoothing"])
            #meshing_options.Set(ModelerUtilities.VARIABLES_SMOOTHING, parameters["JacobiSmoothing"])
            
            self.MeshingParameters.SetOptions(meshing_options)

            self.MeshingParameters.SetOffsetFactor(self.offset_factor)
            self.MeshingParameters.SetAlphaParameter(self.alpha_shape)
                
            self.MeshingParameters.SetReferenceElement(parameters["MeshElement"])
            self.MeshingParameters.SetReferenceCondition("CompositeCondition2D2N")
            #self.MeshingParameters.SetReferenceCondition("WallCondition2D")
                

            #Pre Meshing Processes
            #remove_mesh_nodes = RemoveMeshNodes(self.model_part, self.MeshingParameters, self.echo_level)
            remove_mesh_nodes = RemoveMeshNodesForFluids(self.model_part, self.MeshingParameters, self.echo_level)

            mesh_modeler.SetPreMeshingProcess(remove_mesh_nodes)

            ##refine_mesh_boundary = RefineMeshBoundary(self.model_part, self.RefiningParameters, self.InfoParameters, self.echo_level) commented the 26 04 2016
            #mesh_modeler.SetPreMeshingProcess(refine_mesh_boundary)

            #Post Meshing Processes
            rebuild_mesh_boundary = ReconstructMeshBoundary(self.model_part, self.MeshingParameters, self.echo_level)
            #rebuild_mesh_boundary = ReconstructMeshBoundaryForFluids(self.model_part, self.MeshingParameters, self.echo_level)
           
            
            mesh_modeler.SetPostMeshingProcess(rebuild_mesh_boundary)
        
            print("::[Modeler_Utility]:: Domain Mesh [",parameters["Subdomain"],"] [ Remesh:",parameters["Remesh"],"] [ Refine:",parameters["Refine"],"]" )
            if( parameters["Remesh"] ):
                print("(Type:",parameters["MeshElement"],")")

            mesh_modeler.SetMeshingParameters(self.MeshingParameters)

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
    def InitializeStep(self):

        self.remesh_executed = False

    #
