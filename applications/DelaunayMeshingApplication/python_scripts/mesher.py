from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMesher(main_model_part, meshing_parameters):
    return Mesher(main_model_part, meshing_parameters)

class Mesher(object):

    #
    def __init__(self, main_model_part, meshing_parameters):

        self.echo_level             = 1
        self.main_model_part        = main_model_part
        self.MeshingParameters      = meshing_parameters

        self.model_part = self.main_model_part
        if( self.main_model_part.Name != self.MeshingParameters.GetSubModelPartName() ):
            self.model_part = self.main_model_part.GetSubModelPart(self.MeshingParameters.GetSubModelPartName())

        print(self._class_prefix()+" Ready")

    #
    def Initialize(self, dimension):

        self.dimension   =  dimension

        # set mesher
        if(self.dimension == 2):
            self.mesher = KratosDelaunay.TriangularMesh2DMesher()
        elif(self.dimension == 3):
            self.mesher = KratosDelaunay.TetrahedralMesh3DMesher()

        self.mesher.SetEchoLevel(self.echo_level)
        self.mesher.SetMeshingParameters(self.MeshingParameters)

        self.SetPreMeshingProcesses()
        self.SetPostMeshingProcesses()

        self.mesher.Initialize()


    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, False)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)

        # set mesher flags: to set options for the mesher (triangle 2D, tetgen 3D)
        if( self.dimension == 2 ):
            pass
            #REFINE
            #ADD NODES
            #to add_nodes automatically and refine the mesh ("q"-quality mesh and "a"-area constraint switches)
            # "YYJaqrn" "YJq1.4arn" "Jq1.4arn"
            #refine
            #mesher_flags = "YJq1.4arnQ"
            #refine constrained
            #mesher_flags = "pYJq1.4arnCQ"

            #INSERT NODES
            #to insert a set of given points and refine the mesh
            # "rinYYJQ" "rinYYJQ" "rinJQ" "rinQ"
            #refine
            #mesher_flags = "rinJQ"
            #refine constrained
            #mesher_flags = "rinYYJQ"

            #refine without adding nodes
            #mesher_flags = "YJrnQ"

            #RECONNECT
            #to reconnect a set of points only
            #mesher_flags = "nQP"
            #constrained
            #mesher_flags = "pnBYYQ"

            #BOUNDARY SEARCH
            #to get conectivities, boundaries and neighbours only
            #mesher_flags = "ncEBQ"

        if( self.dimension == 3 ):
            #other flags
            pass

    #
    def SetPreMeshingProcesses(self):

        # The order set is the order of execution:


        # process to refine elements / refine boundary
        refine_mesh_elements  = KratosDelaunay.RefineElementsOnThreshold(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_elements)

        # process to refine boundary / contact boundary
        refine_mesh_boundary = RefineConditions(self.model_part, self.RefiningParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(refine_mesh_boundary)


        # process to remove nodes / remove boundary
        remove_mesh_nodes = KratosDelaunay.RemoveNodes(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPreMeshingProcess(remove_mesh_nodes)


    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:

        #generate new particles
        generate_particles  = KratosDelaunay.GenerateNewNodes(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(generate_particles)

        #select mesh elements
        select_mesh_elements  = KratosDelaunay.SelectElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosDelaunay.GenerateNewElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosDelaunay.GenerateNewConditions(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)

    #
    def FinalizeMeshing(self):

        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        # all flags
        execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT, False)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, False)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)

        self.MeshingParameters.FinalizeMeshing()

    #
    def ExecuteMeshing(self):

        self.InitializeMeshing()  #set execution flags and mesher flags

        self.mesher.ExecuteMeshing(self.model_part)

        self.FinalizeMeshing()    #set execution flags and mesher flags

    #
    def SetEchoLevel(self, echo_level):
        self.echo_level = echo_level

    #
    @classmethod
    def _class_prefix(self):
        header = "::[-------Mesher------]::"
        return header
