from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMeshModeler(main_model_part, custom_settings):
    return MeshModeler(main_model_part, meshing_parameters, mesh_id)

class MeshModeler:
    
    #
    def __init__(self, main_model_part, meshing_parameters, mesh_id): 
        
        self.echo_level        = 0
        self.mesh_id           = mesh_id
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters

    #
    def Initialize(self, imposed_walls, domain_size):
        
        self.domain_size   =  domain_size
        self.imposed_walls =  imposed_walls
        # set mesh modeler
        if(self.domain_size == 2):
            self.mesher = TriangularMesh2DModeler()
        elif(self.domain_size == 3):
            self.mesher = TetrahedralMesh3DModeler()

        self.mesher.SetMeshingParameters(self.MeshingParameters,self.mesh_id)

        self.SetPreMeshingProcess()
        self.SetPostMeshingProcess()    

        self.mesher.Initialize()
        self.mesher.SetEchoLevel(self.echo_level)

    #
    def SetPreMeshingProcesses(self):
        
        # The order set is the order of execution:

        # process to remove nodes / remove boundary
        remove_mesh_nodes = RemoveMeshNodes(self.main_model_part, self.MeshingParameters,  self.mesh_id, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)

        # process to refine elements /refine boundary
        refine_mesh_elements  = RefineMeshElements(self.main_model_part, self.RefiningParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_elements)

        #refine_mesh_boundary = RefineMeshBoundary(self.model_part, self.RefiningParameters, self.mesh_id, self.echo_level)            
        #self.mesher.SetPreMeshingProcess(refine_mesh_boundary)
                
        #set imposed walls (rigid walls)
        rigid_walls_container = BoundingBoxContainer()

        if( self.imposed_walls.RigidWallActive() ):
            rigid_wall_bbox = self.imposed_walls.RigidWallBoundingBoxes()
            for sizei in range(0, len(rigid_wall_bbox)):
                rigid_walls_container.PushBack( rigid_wall_bbox[sizei] )

        #refine_mesh_boundary
        refine_mesh_boundary = ContactRefineMeshBoundary(self.model_part, rigid_walls_container, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_boundary)

        
    #
    def SetPostMeshingProcesses(self):

        #select mesh elements
        generate_particles  = GenerateNewParticles(self.model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPostMeshingProcess(generate_particles)

        #select mesh elements
        select_mesh_elements  = SelectMeshElements(self.model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = BuildMeshElements(self.model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = ReconstructMeshBoundary(self.model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)

    #
    def InitializeMeshing(self):
        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = Flags()

        execution_options.Set(ModelerUtilities.SET_NODES, False)
        execution_options.Set(ModelerUtilities.SET_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SET_FACES, False)  
        execution_options.Set(ModelerUtilities.SELECT_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SELECT_NODES, False)
        execution_options.Set(ModelerUtilities.PASS_ALPHA_SHAPE, False)
        execution_options.Set(ModelerUtilities.ENGAGED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        if( self.domain_size = 2 ):

        #REFINE
        #ADD NODES
        #to add_nodes automatically and refine the mesh ("q"-quality mesh and "a"-area constraint switches)
        # "YYJaqrn" "YJq1.4arn" "Jq1.4arn"
        #refine
        modeler_flags = "YJq1.4arnQ" 
        #refine constrained
        modeler_flags = "pYJq1.4arnCQ"

        #INSERT NODES
        #to insert a set of given points and refine the mesh
        # "rinYYJQ" "rinYYJQ" "rinJQ" "rinQ"
        #refine
        modeler_flags = "rinJQ" 
        #refine constrained
        modeler_flags = "rinYYJQ"

        #refine without adding nodes
        modeler_flags = "YJrnQ" 

        #RECONNECT
        #to reconnect a set of points only
        modeler_flags = "nQP"
        #constrained
        modeler_flags = "pnBYYQ"

        #BOUNDARY SEARCH
        #to get conectivities, boundaries and neighbours only
        modeler_flags = "ncEBQ" 
        
        if( self.domain_size = 3 ):


    #
    def FinalizeMeshing(self):
        
        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = Flags()
        execution_options.Set(ModelerUtilities.SET_NODES, False)
        execution_options.Set(ModelerUtilities.SET_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SET_FACES, False)  
        execution_options.Set(ModelerUtilities.SELECT_ELEMENTS, False)
        execution_options.Set(ModelerUtilities.SELECT_NODES, False)
        execution_options.Set(ModelerUtilities.PASS_ALPHA_SHAPE, False)
        execution_options.Set(ModelerUtilities.ENGAGED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
    #
    def ExecuteMeshing(self):

    
        self.mesher.InitializeMeshing()  #set execution flags and modeler flags
        
        self.mesher.ExecuteMeshing()
        
        self.mesher.FinalizeMeshing()    #set execution flags and modeler flags
