from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemBaseApplication as KratosPfemBase

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mesh modeler (the base class for the modeler derivation)
import mesh_modeler

def CreateMeshModeler(main_model_part, meshing_parameters, mesh_id):
    return TransferModeler(main_model_part, meshing_parameters, mesh_id)

class TransferModeler(mesh_modeler.MeshModeler):
    
    #
    def __init__(self, main_model_part, meshing_parameters, mesh_id): 
        
        self.echo_level        = 0
        self.mesh_id           = mesh_id
        self.main_model_part   = main_model_part 
        self.MeshingParameters = meshing_parameters

        print("Construction of RecoonectModeler finished")
  
    #
    def InitializeMeshing(self):
        
        # set execution flags: to set the options to be executed in methods and processes
        execution_options = KratosMultiphysics.Flags()

        execution_options.Set(KratosPfemBase.ModelerUtilities.INITIALIZE_MESHER_INPUT, True)
        execution_options.Set(KratosPfemBase.ModelerUtilities.FINALIZE_MESHER_INPUT, True)

        execution_options.Set(KratosPfemBase.ModelerUtilities.TRANSFER_KRATOS_NODES_TO_MESHER, True)
        execution_options.Set(KratosPfemBase.ModelerUtilities.TRANSFER_KRATOS_ELEMENTS_TO_MESHER, True)
        execution_options.Set(KratosPfemBase.ModelerUtilities.TRANSFER_KRATOS_NEIGHBOURS_TO_MESHER, True)
                              
        execution_options.Set(KratosPfemBase.ModelerUtilities.SELECT_TESSELLATION_ELEMENTS, True)
        execution_options.Set(KratosPfemBase.ModelerUtilities.KEEP_ISOLATED_NODES, False)

        self.MeshingParameters.SetExecutionOptions(execution_options)
        

    #
    def SetPreMeshingProcesses(self):
        #nothing to do: only transfer
        pass


   #
    def SetPostMeshingProcesses(self):
        
        #nothing to do: only reconnection
        #select mesh elements
        #select_mesh_elements  = KratosPfemBase.SelectMeshElements(self.main_model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        #self.mesher.SetPostMeshingProcess(select_mesh_elements)

        #rebuild elements
        rebuild_mesh_elements = KratosPfemBase.BuildMeshElements(self.main_model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_elements)

        #rebuild boundary
        rebuild_mesh_boundary = KratosPfemBase.ReconstructMeshBoundary(self.main_model_part, self.MeshingParameters, self.mesh_id, self.echo_level)
        self.mesher.SetPostMeshingProcess(rebuild_mesh_boundary)

