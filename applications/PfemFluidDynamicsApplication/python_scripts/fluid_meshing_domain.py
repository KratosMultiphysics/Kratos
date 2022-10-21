#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.DelaunayMeshingApplication as KratosDelaunay
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
from importlib import import_module


def CreateMeshingDomain(main_model_part, custom_settings):
    return FluidMeshingDomain(main_model_part, custom_settings)

class FluidMeshingDomain(object):

    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part.
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the mesher is already filled
    def __init__(self, main_model_part, custom_settings):

        self.echo_level = 1
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
	    "python_module": "meshing_domain",
            "model_part_name": "model_part_name",
            "alpha_shape": 1.25,
            "meshing_strategy":{
               "python_module": "meshing_strategy",
               "remesh": false,
               "refine": false,
               "transfer": false,
               "reference_element_type": "Element2D3N",
               "reference_condition_type": "CompositeCondition2D2N"
            },
            "spatial_bounding_box":{
                "use_bounding_box" : true,
                "initial_time"     : 0.0,
                "final_time"       : 1000.0,
                "upper_point"      : [10,10,10],
                "lower_point"      : [-10,-10,-10]
            },
            "spatial_refining_box"            : {
                    "use_refining_box" : false,
                    "mesh_size"        : 0.1,
                    "initial_time"     : 0.0,
                    "final_time"       : 1,
                    "upper_point"      : [10,10,10],
                    "lower_point"      : [-10,-10,-10]
            },
            "spatial_refining_box_list"            : [{
                    "use_refining_box"    : false,
                    "transition_elements" : 4,
                    "initial_time"        : 0.0,
                    "final_time"          : 1,
                    "upper_point"         : [10,10,10],
                    "lower_point"         : [-10,-10,-10]
            }]
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the meshing strategy
        python_module_name = "KratosMultiphysics.PfemFluidDynamicsApplication"
        full_module_name = python_module_name + "." + self.settings["meshing_strategy"]["python_module"].GetString()
        meshing_module = import_module(full_module_name)
        #meshing_module = __import__(self.settings["meshing_strategy"]["python_module"].GetString())
        self.MeshingStrategy = meshing_module.CreateMeshingStrategy(self.main_model_part, self.settings["meshing_strategy"])

        self.active_remeshing = False
        if( self.settings["meshing_strategy"]["remesh"].GetBool() or self.settings["meshing_strategy"]["transfer"].GetBool() ):
            self.active_remeshing = True

        print("::[Meshing_Domain]:: (",self.settings["model_part_name"].GetString()," ) -BUILT-")


    ####

    def Initialize(self):

        self.dimension = self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION]

        # Set MeshingParameters
        self.SetMeshingParameters()

        # Meshing Stratety
        self.MeshingStrategy.SetEchoLevel(self.echo_level)
        self.MeshingStrategy.Initialize(self.MeshingParameters, self.dimension)

    ####

    #
    def SetInfoParameters(self):

        # Create InfoParameters
        self.InfoParameters  = KratosDelaunay.MeshingInfoParameters()
        self.InfoParameters.Initialize()

    #
    def SetRefiningParameters(self):

        # Create RefiningParameters
        self.RefiningParameters = KratosDelaunay.RefiningParameters()
        self.RefiningParameters.Initialize()

        # parameters
        self.RefiningParameters.SetAlphaParameter(self.settings["alpha_shape"].GetDouble())
        number_of_refining_boxes=self.settings["spatial_refining_box_list"].size()
        self.MeshingParameters.InitializeRefiningBoxParameters(number_of_refining_boxes)       
        # set mesh refinement in box
        index=0    
        for item in self.settings["spatial_refining_box_list"]:
            refining_box_list = item
            self.MeshingParameters.SetUseRefiningBox(index,refining_box_list["use_refining_box"].GetBool()) 
            self.MeshingParameters.SetRefiningBoxElementsInTransitionZone(index,refining_box_list["transition_elements"].GetInt())
            self.MeshingParameters.SetRefiningBoxTimeInterval(index,refining_box_list["initial_time"].GetDouble(),refining_box_list["final_time"].GetDouble())
            self.MeshingParameters.SetRefiningBoxMinimumPoint(index,refining_box_list["lower_point"][0].GetDouble(),refining_box_list["lower_point"][1].GetDouble(),refining_box_list["lower_point"][2].GetDouble()) 
            self.MeshingParameters.SetRefiningBoxShiftedMinimumPoint(index,refining_box_list["lower_point"][0].GetDouble(),refining_box_list["lower_point"][1].GetDouble(),refining_box_list["lower_point"][2].GetDouble()) 
            self.MeshingParameters.SetRefiningBoxMaximumPoint(index,refining_box_list["upper_point"][0].GetDouble(),refining_box_list["upper_point"][1].GetDouble(),refining_box_list["upper_point"][2].GetDouble()) 
            self.MeshingParameters.SetRefiningBoxShiftedMaximumPoint(index,refining_box_list["upper_point"][0].GetDouble(),refining_box_list["upper_point"][1].GetDouble(),refining_box_list["upper_point"][2].GetDouble()) 
            index+=1   

        removing_options = KratosMultiphysics.Flags()

        #remove nodes
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_NODES, True)
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_NODES_ON_DISTANCE, True)
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_NODES_ON_ERROR, False)
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_NODES_ON_THRESHOLD, False)

        #remove boundary
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_BOUNDARY_NODES, False)
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_BOUNDARY_NODES_ON_DISTANCE, False)
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_BOUNDARY_NODES_ON_ERROR, False)
        removing_options.Set(KratosDelaunay.MesherUtilities.REMOVE_BOUNDARY_NODES_ON_THRESHOLD, False)

        refining_options = KratosMultiphysics.Flags()
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE, self.settings["meshing_strategy"]["refine"].GetBool())

        #refine elements
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_ELEMENTS, False)
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_ELEMENTS_ON_DISTANCE, False)
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_ELEMENTS_ON_ERROR, False)
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_ELEMENTS_ON_THRESHOLD, False)

        #refine boundary
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_BOUNDARY, False)
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_BOUNDARY_ON_DISTANCE, False)
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_BOUNDARY_ON_ERROR, False)
        refining_options.Set(KratosDelaunay.MesherUtilities.REFINE_BOUNDARY_ON_THRESHOLD, False)

        self.RefiningParameters.SetRefiningOptions(refining_options)
        self.RefiningParameters.SetRemovingOptions(removing_options)

    #
    def SetMeshingParameters(self):

        self.MeshingParameters = KratosDelaunay.MeshingParameters()
        self.MeshingParameters.Initialize()

        self.MeshingParameters.SetSubModelPartName(self.settings["model_part_name"].GetString())

        if(self.active_remeshing):

            self.MeshingParameters.SetAlphaParameter(self.settings["alpha_shape"].GetDouble())

            self.SetInfoParameters()
            self.SetRefiningParameters()

            self.MeshingParameters.SetInfoParameters(self.InfoParameters)
            self.MeshingParameters.SetRefiningParameters(self.RefiningParameters)


            bounding_box = self.settings["spatial_bounding_box"]
            if(bounding_box["use_bounding_box"].GetBool()):
              self.MeshingParameters.SetUseBoundingBox(True)
              self.MeshingParameters.SetBoundingBoxLowerPoint(bounding_box["lower_point"][0].GetDouble(),bounding_box["lower_point"][1].GetDouble(),bounding_box["lower_point"][2].GetDouble())
              self.MeshingParameters.SetBoundingBoxUpperPoint(bounding_box["upper_point"][0].GetDouble(),bounding_box["upper_point"][1].GetDouble(),bounding_box["upper_point"][2].GetDouble())
              self.MeshingParameters.SetBoundingBoxTimeInterval(bounding_box["initial_time"].GetDouble(),bounding_box["final_time"].GetDouble())

    #
    def ExecuteMeshing(self):

        if( self.active_remeshing ):
            self.MeshingStrategy.GenerateMesh()

    #
    def Check(self):

        # set mesher utilities
        self.mesher_utils = KratosDelaunay.MesherUtilities()

        critical_radius = self.mesher_utils.CheckCriticalRadius(self.main_model_part,critical_mesh_size)
        print(" CriticalRadius ", critical_radius)

    #
    def Active(self):
        return self.active_remeshing

    #
    def SetEchoLevel(self, echo_level):
        self.echo_level = echo_level

    #
    def GetVariables(self):

        nodal_variables = []
        transfer_variables = self.settings["elemental_variables_to_transfer"]
        for i in range(0, transfer_variables.size() ):
            nodal_variables.append(transfer_variables[i].GetString())

        return nodal_variables

    #
    def ComputeAverageMeshParameters(self):

        MesherUtils = KratosDelaunay.MesherUtilities();
        self.domain_volume =  MesherUtils.ComputeModelPartVolume(self.main_model_part)
        self.element_mean_volume = 0

        number_of_elements =  self.main_model_part.NumberOfElements()
        nodes_for_element  =  self.main_model_part.ProcessInfo[KratosMultiphysics.SPACE_DIMENSION] + 1

        if(number_of_elements != 0):
            self.element_mean_volume = self.domain_volume/float(number_of_elements*nodes_for_element)

        self.RefiningParameters.SetMeanVolume(self.element_mean_volume)

    #
    def GetMeanVolume(self):

        return self.element_mean_volume

    #
    def GetTotalVolume(self):

        return self.domain_volume

    #
    def ComputeInitialAverageMeshParameters(self):

        self.mesh_parameters =  KratosPfemFluid.ComputeAveragePfemMeshParameters(self.main_model_part, self.MeshingParameters,self.echo_level)
        self.mesh_parameters.Execute()

        # numFluid=0
        # mean_nodal_h=0
        # for node in self.main_model_part.Nodes:
        #     if (node.Is(KratosMultiphysics.FLUID)):
        #         numFluid+=1
        #         nodal_h=node.GetSolutionStepValue(KratosMultiphysics.NODAL_H)
        #         mean_nodal_h+=nodal_h

        # mean_nodal_h*=1.0/numFluid;

        # self.RefiningParameters.SetCriticalRadius(mean_nodal_h)
        # self.RefiningParameters.SetInitialRadius(mean_nodal_h)

        # delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        # self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.INITIAL_DELTA_TIME,delta_time)
        # self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.CURRENT_DELTA_TIME,delta_time)
        # self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.PREVIOUS_DELTA_TIME,delta_time)
        # self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.TIME_INTERVAL_CHANGED,False)

    def SetTimeDataOnProcessInfo(self):

        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.INITIAL_DELTA_TIME,delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.CURRENT_DELTA_TIME,delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.PREVIOUS_DELTA_TIME,delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosPfemFluid.TIME_INTERVAL_CHANGED,False)


    #
