
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid

from KratosMultiphysics.DelaunayMeshingApplication import meshing_strategy

from importlib import import_module

def CreateMeshingStrategy(main_model_part, custom_settings):
    return FluidMeshingStrategy(main_model_part, custom_settings)

class FluidMeshingStrategy(meshing_strategy.MeshingStrategy):

    #
    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "python_module": "meshing_strategy",
             "meshing_frequency": 0.0,
             "refine": false,
             "refine_module" : "KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_complete_mesher",
             "remesh": false,
             "remesh_module" : "KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_keeping_nodes_mesher",
             "transfer" : false,
             "transfer_module" : "KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_none_mesher",
             "reference_element_type": "Element2D3N",
             "reference_condition_type": "CompositeCondition2D3N"
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.echo_level = 0

    #
    def Initialize(self,meshing_parameters,dimension):

        #meshing parameters
        self.MeshingParameters = meshing_parameters

        meshing_options = KratosMultiphysics.Flags()

        meshing_options.Set(KratosDelaunay.MesherUtilities.REMESH, self.settings["remesh"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.REFINE, self.settings["refine"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.TRANSFER, self.settings["transfer"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.RECONNECT, False)
        meshing_options.Set(KratosDelaunay.MesherUtilities.CONSTRAINED, False)
        meshing_options.Set(KratosDelaunay.MesherUtilities.MESH_SMOOTHING, False)
        meshing_options.Set(KratosDelaunay.MesherUtilities.VARIABLES_SMOOTHING, False)

        self.MeshingParameters.SetOptions(meshing_options)
        self.MeshingParameters.SetReferenceElement(self.settings["reference_element_type"].GetString())
        self.MeshingParameters.SetReferenceCondition(self.settings["reference_condition_type"].GetString())

        #set variables to global transfer
        self.MeshDataTransfer   = KratosDelaunay.MeshDataTransferUtilities()
        self.TransferParameters = KratosDelaunay.TransferParameters()
        self.global_transfer    = False

        #mesh meshers for the current strategy
        self.meshers = []

        #configure meshers:
        self.SetMeshers();

        self.model_part = self.main_model_part
        if( self.main_model_part.Name != self.MeshingParameters.GetSubModelPartName() ):
            self.model_part = self.main_model_part.GetSubModelPart(self.MeshingParameters.GetSubModelPartName())

        for mesher in self.meshers:
            mesher.SetEchoLevel(self.echo_level)
            mesher.Initialize(dimension)

        self.number_of_nodes      = 0
        self.number_of_elements   = 0
        self.number_of_conditions = 0

    #
    def SetMeshers(self):

        meshers_list = []
        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):
            meshers_list.append(self.settings["refine_module"].GetString())
        elif( self.settings["remesh"].GetBool() ):
            meshers_list.append(self.settings["remesh_module"].GetString())
        elif( self.settings["transfer"].GetBool() ):
            meshers_list.append(self.settings["transfer_module"].GetString())

        for mesher in meshers_list:
            full_module_name = mesher
            meshing_module = import_module(full_module_name)
            new_mesher = meshing_module.CreateMesher(self.main_model_part,self.MeshingParameters)
            self.meshers.append(new_mesher)

    #
    def GetMeshers(self):

        meshers_list = []

        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):

            meshers_list.append("pre_refining_mesher")
            meshers_list.append("post_refining_mesher")

        elif( self.settings["remesh"].GetBool() ):

            meshers_list.append("reconnect_mesher")

        return meshers_list
    #
