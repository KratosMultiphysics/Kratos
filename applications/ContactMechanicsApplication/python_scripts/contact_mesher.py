from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh mesher (the base class for the mesher derivation)
import mesher

def CreateMesher(main_model_part, meshing_parameters):
    return ContactMesher(main_model_part, meshing_parameters)

class ContactMesher(mesher.Mesher):

    #
    def __init__(self, main_model_part, meshing_parameters):

        mesher.Mesher.__init__(self, main_model_part, meshing_parameters)
        self.echo_level = 0

    #
    def Initialize(self, dimension):

        self.dimension = dimension

        # set mesh mesher
        if(self.dimension == 2):
            self.mesher = KratosContact.ContactDomain2DMesher()
        elif(self.dimension == 3):
            self.mesher = KratosContact.ContactDomain3DMesher()

        self.mesher.SetEchoLevel(self.echo_level)
        self.mesher.SetMeshingParameters(self.MeshingParameters)

        self.SetPreMeshingProcesses()
        self.SetPostMeshingProcesses()

        self.mesher.Initialize()

    #
    def InitializeMeshing(self):

        self.MeshingParameters.InitializeMeshing()

        meshing_options = self.MeshingParameters.GetOptions()

        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosDelaunay.MesherUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
            execution_options.Set(KratosDelaunay.MesherUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)

        execution_options.Set(KratosDelaunay.MesherUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosDelaunay.MesherUtilities.KEEP_ISOLATED_NODES, True)


        self.MeshingParameters.SetExecutionOptions(execution_options)

        # set mesher flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT

        mesher_flags = ""
        mesher_info  = "Reconnect a cloud of points"
        if( self.dimension == 2 ):

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                mesher_flags = "pBYYQ"
            else:
                mesher_flags = "QNP"


        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosDelaunay.MesherUtilities.CONSTRAINED) ):
                mesher_flags = "pMYYCJFQ"     #tetgen 1.5.0
                #mesher_flags = "pJFBMYYCCQu0"  #tetgen 1.4.3
                #mesher_flags = "pJFBMYYCCQ"  #tetgen 1.5.0
            else:
                mesher_flags = "JFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(mesher_flags)
        self.MeshingParameters.SetTessellationInfo(mesher_info)

    #
    def SetPreMeshingProcesses(self):

        # The order set is the order of execution:

        # clear contact conditions
        clear_contact_conditions= KratosContact.ClearContactConditions(self.model_part, self.echo_level)
        self.mesher.SetPreMeshingProcess(clear_contact_conditions)

        #print GiD mesh output for checking purposes
        #print_output_mesh = KratosDelaunay.PrintMeshOutput(self.model_part, self.MeshingParameters, "input", self.echo_level)
        #self.mesher.SetPreMeshingProcess(print_output_mesh)
    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:

        #print GiD mesh output for checking purposes (current print)
        #print_output_mesh = KratosDelaunay.PrintMeshOutput(self.model_part, self.MeshingParameters, "output", self.echo_level)
        #self.mesher.SetPostMeshingProcess(print_output_mesh)

        #select mesh elements
        select_mesh_elements  = KratosDelaunay.SelectElements(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(select_mesh_elements)

        # build contact conditions
        build_contact_conditions= KratosContact.BuildContactConditions(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(build_contact_conditions)

    #
    def FinalizeMeshing(self):

        # reset execution flags: to unset the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()
        # all false
        self.MeshingParameters.SetExecutionOptions(execution_options)

        self.MeshingParameters.FinalizeMeshing()

    #
    @classmethod
    def _class_prefix(self):
        header = "::[---Contact Mesher--]::"
        return header
