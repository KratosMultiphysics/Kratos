from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemApplication as KratosPfem
import KratosMultiphysics.ContactMechanicsApplication as KratosContact

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh modeler (the base class for the modeler derivation)
import mesh_modeler

def CreateMeshModeler(main_model_part, meshing_parameters):
    return ContactModeler(main_model_part, meshing_parameters)

class ContactModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters): 
        
        mesh_modeler.MeshModeler.__init__(self, main_model_part, meshing_parameters)        

        print("::[Contact_Mesh_Modeler]:: -BUILT-")
  
    #
    def Initialize(self, dimension):
        
        self.dimension   =  dimension

        # set mesh modeler
        if(self.dimension == 2):
            self.mesher = KratosContact.ContactDomain2DModeler()
        elif(self.dimension == 3):
            self.mesher = KratosContact.ContactDomain3DModeler()

        self.mesher.SetEchoLevel(self.echo_level)
        self.mesher.SetMeshingParameters(self.MeshingParameters)

        self.SetPreMeshingProcesses()
        self.SetPostMeshingProcesses()    

        self.mesher.Initialize()

    #
    def InitializeMeshing(self):

        meshing_options = self.MeshingParameters.GetOptions()
        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosPfem.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosPfem.ModelerUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, False)
        execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, False)

        if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
            execution_options.Set(KratosPfem.ModelerUtilities.TRANSFER_KRATOS_FACES_TO_MESHER, True)
                              
        execution_options.Set(KratosPfem.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosPfem.ModelerUtilities.KEEP_ISOLATED_NODES, True)


        self.MeshingParameters.SetExecutionOptions(execution_options)
        
        # set modeler flags: to set options for the mesher (triangle 2D, tetgen 3D)
        # RECONNECT
            
        modeler_flags = ""
        modeler_info  = "Reconnect a cloud of points"
        if( self.dimension == 2 ):
           
            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                modeler_flags = "pBYYQ"  
            else:
                modeler_flags = "QNP"

            
        elif( self.dimension == 3 ):

            if( meshing_options.Is(KratosPfem.ModelerUtilities.CONSTRAINED) ):
                #modeler_flags = "pMYYCJFu0"     #tetgen 1.5.0
                modeler_flags = "pJFBMYYCCQu0"  #tetgen 1.4.3
                #modeler_flags = "pJFBMYYCCQ"  #tetgen 1.5.0
            else:
                modeler_flags = "JFMQO4/4"

        self.MeshingParameters.SetTessellationFlags(modeler_flags)
        self.MeshingParameters.SetTessellationInfo(modeler_info)

    #
    def SetPreMeshingProcesses(self):
 
        # The order set is the order of execution:

        # clear contact conditions
        clear_contact_conditions= KratosContact.ClearContactConditions(self.model_part, self.echo_level)
        self.mesher.SetPreMeshingProcess(clear_contact_conditions)

    #
    def SetPostMeshingProcesses(self):

        # The order set is the order of execution:

        #print GiD mesh output for checking purposes
        print_output_mesh = KratosPfem.PrintOutputMeshProcess(self.model_part, self.MeshingParameters, self.echo_level)
        self.mesher.SetPostMeshingProcess(print_output_mesh)        
        
        #select mesh elements
        select_mesh_elements  = KratosPfem.SelectMeshElements(self.model_part, self.MeshingParameters, self.echo_level)
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
