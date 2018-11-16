from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateMeshingStrategy(main_model_part, custom_settings):
    return MeshingStrategy(main_model_part, custom_settings)

class MeshingStrategy(object):

    #
    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
             "python_module": "meshing_strategy",
             "meshing_frequency": 0.0,
             "remesh": false,
             "refine": false,
             "reconnect": false,
             "transfer": false,
             "constrained": false,
             "mesh_smoothing": false,
             "variables_smoothing": false,
             "elemental_variables_to_smooth":[ "DETERMINANT_F" ],
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
        meshing_options.Set(KratosDelaunay.MesherUtilities.RECONNECT, self.settings["reconnect"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.TRANSFER, self.settings["transfer"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.CONSTRAINED, self.settings["constrained"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.MESH_SMOOTHING, self.settings["mesh_smoothing"].GetBool())
        meshing_options.Set(KratosDelaunay.MesherUtilities.VARIABLES_SMOOTHING, self.settings["variables_smoothing"].GetBool())

        self.MeshingParameters.SetOptions(meshing_options)
        self.MeshingParameters.SetReferenceElement(self.settings["reference_element_type"].GetString())
        self.MeshingParameters.SetReferenceCondition(self.settings["reference_condition_type"].GetString())

        #set variables to global transfer
        self.MeshDataTransfer   = KratosDelaunay.MeshDataTransferUtilities()
        self.TransferParameters = KratosDelaunay.TransferParameters()
        self.global_transfer    = False
        if( self.settings["variables_smoothing"].GetBool() == True ):
            self.global_transfer = True
            transfer_variables = self.settings["elemental_variables_to_smooth"]
            #for variable in transfer_variables:
            #    self.TransferParameters.SetVariable( KratosMultiphysics.KratosGlobals.GetVariable( variable.GetString() ) )
            for i in range(0, transfer_variables.size() ):
                self.TransferParameters.SetVariable(KratosMultiphysics.KratosGlobals.GetVariable(transfer_variables[i].GetString()))


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

        print(self._class_prefix()+" Ready")

    #
    def GetMeshers(self):

        meshers_list = []

        if( self.settings["remesh"].GetBool() and self.settings["refine"].GetBool() ):

            meshers_list.append("pre_refining_mesher")
            meshers_list.append("post_refining_mesher")

        elif( self.settings["remesh"].GetBool() ):

            meshers_list.append("reconnect_mesher")

        elif( self.settings["transfer"].GetBool() ):

            meshers_list.append("transfer_mesher")

        return meshers_list

    #
    def SetMeshers(self):

        meshers_list = self.GetMeshers()

        if( self.echo_level > 0 ):
            print(self._class_prefix()+"  ["+self.MeshingParameters.GetSubModelPartName()+" model part ] (REMESH:",self.settings["remesh"].GetBool(),"/ REFINE:",self.settings["refine"].GetBool(),"/ TRANSFER:",self.settings["transfer"].GetBool(),")")

        for mesher in meshers_list:
            meshing_module =__import__(mesher)
            new_mesher = meshing_module.CreateMesher(self.main_model_part,self.MeshingParameters)
            self.meshers.append(new_mesher)

    #
    def SetMeshInfo(self):

        info_parameters = self.MeshingParameters.GetInfoParameters()

        number_of_new_nodes = self.main_model_part.NumberOfNodes() - info_parameters.GetNumberOfNodes()
        number_of_new_elements = self.main_model_part.NumberOfElements() - info_parameters.GetNumberOfElements()
        number_of_new_conditions = self.main_model_part.NumberOfConditions() - info_parameters.GetNumberOfConditions()

        if( number_of_new_nodes > 0 ):
            info_parameters.SetNumberOfNewNodes(number_of_new_nodes)
        else:
            info_parameters.SetNumberOfNewNodes(0)

        if( number_of_new_elements > 0 ):
            info_parameters.SetNumberOfNewElements(number_of_new_elements)
        else:
            info_parameters.SetNumberOfNewElements(0)

        if( number_of_new_conditions > 0 ):
            info_parameters.SetNumberOfNewConditions(number_of_new_conditions)
        else:
            info_parameters.SetNumberOfNewConditions(0)


        info_parameters.SetNumberOfNodes(self.main_model_part.NumberOfNodes())
        info_parameters.SetNumberOfElements(self.main_model_part.NumberOfElements())
        info_parameters.SetNumberOfConditions(self.main_model_part.NumberOfConditions())



    #
    def InitializeMeshGeneration(self):

        info_parameters = self.MeshingParameters.GetInfoParameters()
        info_parameters.Initialize()

        self.SetMeshInfo()

        if( self.global_transfer == True ):
            if( self.echo_level > 0 ):
                print(self._class_prefix()+" Elements To Nodes transfer ")
            self.MeshDataTransfer.TransferElementalValuesToNodes(self.TransferParameters,self.model_part)

    #
    def FinalizeMeshGeneration(self):

        self.SetMeshInfo()

        info_parameters    = self.MeshingParameters.GetInfoParameters()
        smoothing_required = info_parameters.CheckMechanicalSmoothing()

        refining_parameters = self.MeshingParameters.GetRefiningParameters()

        if( self.global_transfer == True ):
            self.MeshDataTransfer.TransferNodalValuesToElements(self.TransferParameters,self.model_part)

    #        if( self.global_transfer == True ):
    #            if(smoothing_required):
    #                #smooth only on selected part based on a threshold variable
    #                print(" smooth only on threshold ")
    #                self.MeshDataTransfer.TransferNodalValuesToElementsOnThreshold(self.TransferParameters,refining_parameters,self.model_part)
    #            else:
    #                #smooth all domain
    #                print(" smooth all domain ")
    #                self.MeshDataTransfer.TransferNodalValuesToElements(self.TransferParameters,self.model_part)



    #
    def GenerateMesh(self):

        self.InitializeMeshGeneration()

        for mesher in self.meshers:
            mesher.ExecuteMeshing()

        self.FinalizeMeshGeneration()

    #
    def SetEchoLevel(self, echo_level):
        self.echo_level = echo_level

    #
    @classmethod
    def _class_prefix(self):
        header = "::[--Meshing Strategy-]::"
        return header
