import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow
import KratosMultiphysics.ShapeOptimizationApplication as KSO
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
import time
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LevelSetRemeshingProcess(Model, settings["Parameters"])

def RotateModelPart(origin, angle, model_part):
    ox,oy,oz=origin
    for node in model_part.Nodes:
        node.X = ox+math.cos(angle)*(node.X - ox)-math.sin(angle)*(node.Y - oy)
        node.Y = oy+math.sin(angle)*(node.X - ox)+math.cos(angle)*(node.Y - oy)

## All the processes python should be derived from "Process"
class LevelSetRemeshingProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name": "insert_model_part",
                "skin_model_part_name": "insert_skin_model_part",
                "maximum_iterations": 1,
                "problem_name": "",
                "input_type": "calculate",
                "update_coefficient": 0.5,
                "remeshing_flag": false,
                "perform_moving": true,
                "ray_casting_tolerance": 1e-9,
                "initial_angle_of_attack" : 0.0,
                "moving_parameters":    {
                    "origin"                        : [0.0,0.0,0.0],
                    "rotation_angle"                : 0.0,
                    "sizing_multiplier"             : 1.0
                },
                "metric_parameters":  {
                    "minimal_size"                         : 5e-3,
                    "maximal_size"                         : 1.0,
                    "sizing_parameters": {
                        "reference_variable_name"               : "DISTANCE",
                        "boundary_layer_max_distance"           : 1.0,
                        "interpolation"                         : "constant"
                    },
                    "enforce_current"                      : true,
                    "anisotropy_remeshing"                 : false,
                    "anisotropy_parameters": {
                        "reference_variable_name"          : "DISTANCE",
                        "hmin_over_hmax_anisotropic_ratio"      : 0.5,
                        "boundary_layer_max_distance"           : 1,
                        "interpolation"                         : "Linear"
                    }
                },
                "distance_modification_parameters":{
                    "distance_threshold"                          : 0.001,
                    "check_at_each_time_step"                     : true,
                    "avoid_almost_empty_elements"                 : true,
                    "deactivate_full_negative_elements"           : true,
                    "full_negative_elements_fixed_variables_list" : ["VELOCITY_POTENTIAL","AUXILIARY_VELOCITY_POTENTIAL"]
                }
            }  """ );
        settings.ValidateAndAssignDefaults(default_parameters)
        self.model=Model
        self.main_model_part = Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.skin_model_part_name=settings["skin_model_part_name"].GetString()

        '''Remeshing loop parameters'''
        self.do_remeshing = settings["remeshing_flag"].GetBool()
        self.step = 0
        self.max_iter = settings["maximum_iterations"].GetInt()
        self.update_coefficient = settings["update_coefficient"].GetDouble()
        self.problem_name = settings["problem_name"].GetString()
        self.input_type = settings["input_type"].GetString()

        self.moving_parameters = settings["moving_parameters"]
        self.perform_moving = settings["perform_moving"].GetBool()
        # Synchronizing parameters for the wake process
        if self.moving_parameters.Has("rotation_point"):
            self.main_model_part.ProcessInfo.SetValue(CompressiblePotentialFlow.WAKE_ORIGIN, self.moving_parameters["rotation_point"].GetVector())
        else:
            self.main_model_part.ProcessInfo.SetValue(CompressiblePotentialFlow.WAKE_ORIGIN, self.moving_parameters["origin"].GetVector())
        self.main_model_part.ProcessInfo.SetValue(CompressiblePotentialFlow.ROTATION_ANGLE, self.moving_parameters["rotation_angle"].GetDouble()+settings["initial_angle_of_attack"].GetDouble())

        self.metric_parameters = settings["metric_parameters"]
        self.distance_modification_parameters = settings["distance_modification_parameters"]
        # self.ray_casting_tolerance = settings["ray_casting_tolerance"].GetDouble()
        self.ray_casting_tolerance = 1e-16
        KratosMultiphysics.Logger.PrintInfo("WARNING", "Hard coding ray casting  tolerance to 1e-16")
    def ExecuteInitialize(self):
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Executing Initialize Geometry')
        self._InitializeSkinModelPart()
        # self.wake_model_part = self.model.CreateModelPart("aux_wake")
        # self.target_h_wake = 0.2
        # self._DefineWakeModelPart()

        ini_time=time.time()
        # self._CalculateDiscontinuousDistanceAndComputeWakeMetric()
        self._CalculateDistance()

        ini_time=time.time()
        if self.do_remeshing:
            print("EXECUTING LEVEL SET REMESHING. Source nnodes:", self.main_model_part.NumberOfNodes())
        while self.step < self.max_iter and self.do_remeshing:
            self.step += 1
            KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','##### Executing refinement #', self.step, ' #####')
            self._ExtendDistance()
            self._RefineMesh()
            # if self.step < self.max_iter:
                # self._CalculateDiscontinuousDistanceAndComputeWakeMetric()
            self._CalculateDistance()
            self._UpdateParameters()

        self._ModifyFinalDistance()
        self._CopyAndDeleteDefaultDistance()


        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Elapsed time: ',time.time()-ini_time)

    def _InitializeSkinModelPart(self):
        ''' This function loads and moves the skin_model_part in the main_model_part to the desired initial point (origin).
            It also rotates the skin model part around the origin point according to the rotation_angle'''

        ini_time=time.time()
        if not self.model.HasModelPart("skin"):
            self.skin_model_part=self.model.CreateModelPart("skin")
            self.skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            self.skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
            self.skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_SENSITIVITY)
            self.skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
            self.skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
            self.skin_model_part.AddNodalSolutionStepVariable(KSO.DF1DX_MAPPED)
            # Reading skin model part
            KratosMultiphysics.ModelPartIO(self.skin_model_part_name).ReadModelPart(self.skin_model_part)
        else:
            self.skin_model_part = self.model.GetModelPart("skin")
            print("GETTING SKIN MODEL PART")
        if self.perform_moving:
            # Moving and rotating the skin model part
            angle=math.radians(-self.moving_parameters["rotation_angle"].GetDouble())
            self.moving_parameters["rotation_angle"].SetDouble(angle)

            CompressiblePotentialFlow.MoveModelPartProcess(self.skin_model_part, self.moving_parameters).Execute()

        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','InitializeSkin time: ',time.time()-ini_time)

    def _CalculateDistance(self):
        ''' This function calculate the distance to skin for every node in the main_model_part.'''
        ini_time=time.time()
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.main_model_part, self.skin_model_part,self.ray_casting_tolerance).Execute()
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','CalculateDistance time: ',time.time()-ini_time)

    def _ExtendDistance(self):
        ''' This function extends the distance field to all the nodes of the main_model_part in order to
            remesh the background mesh.'''
        ini_time=time.time()
        # Construct the variational distance calculation process
        maximum_iterations = 2 #TODO: Make this user-definable

        ###Defining linear solver to be used by the variational distance process###
        from KratosMultiphysics import python_linear_solver_factory #Linear solver for variational distance process
        linear_solver_settings=KratosMultiphysics.Parameters("""
        {
            "solver_type": "amgcl",
            "max_iteration": 200,
            "gmres_krylov_space_dimension": 100,
            "smoother_type":"ilu0",
            "coarsening_type":"ruge_stuben",
            "coarse_enough" : 5000,
            "krylov_type": "lgmres",
            "tolerance": 1e-8,
            "verbosity": 0,
            "scaling": false
        }""")

        linear_solver = python_linear_solver_factory.ConstructSolver(linear_solver_settings)
        variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
            self.main_model_part,
            linear_solver,
            maximum_iterations)
        variational_distance_process.Execute()

        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Variational distance process time: ',time.time()-ini_time)


    def _RefineMesh(self):
        ''' This function remeshes the main_model_part according to the distance, using the MMG process from the MeshingApplication.
            In order to perform the refinement, it is needed to calculate the distance gradient, the initial nodal_h and the level_set metric.
        '''
        ini_time=time.time()
        # for elem in self.main_model_part.Elements:
        #     is_negative = True
        #     elem.Set(KratosMultiphysics.SOLID, False)
        #     for node in elem.GetNodes():
        #         distance=node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
        #         if distance > 0.0:
        #             is_negative = False
        #     if not is_negative:
        #         elem.Set(KratosMultiphysics.SOLID, True)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()


        # gid_output = GiDOutputProcess(self.main_model_part,
        #                             "distance_"+str(self.step),
        #                             KratosMultiphysics.Parameters("""
        #                                 {
        #                                     "result_file_configuration" : {
        #                                         "gidpost_flags": {
        #                                             "GiDPostMode": "GiD_PostBinary",
        #                                             "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                             "WriteConditionsFlag": "WriteConditions",
        #                                             "MultiFileFlag": "SingleFile"
        #                                         },
        #                                         "nodal_results" : ["DISTANCE", "DISTANCE_GRADIENT"],
        #                                         "nodal_nonhistorical_results": ["METRIC_TENSOR_2D","TEMPERATURE"]
        #                                     }
        #                                 }
        #                                 """)
        #                             )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,self.main_model_part.Nodes)
        # for node in self.main_model_part.Nodes:
        #     final_size = node.GetValue(KratosMultiphysics.NODAL_H)
        #     metric_tensor_2d = KratosMultiphysics.Vector(3, 0.0)
        #     metric_tensor_2d[0] =1/final_size/final_size
        #     metric_tensor_2d[1] =1/final_size/final_size
        #     node.SetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D, metric_tensor_2d)


        metric_parameters_hessian = KratosMultiphysics.Parameters("""
                {
                    "minimal_size"                        : 0.0000001,
                    "maximal_size"                        : 1000000.0,
                    "enforce_current"                     : false,
                    "hessian_strategy_parameters":
                    {
                        "non_historical_metric_variable"  : true,
                        "estimate_interpolation_error"    : false,
                        "interpolation_error"             : 5e-3
                    },
                    "anisotropy_remeshing"                : false
                }
                """)

        # for node in self.main_model_part.Nodes:
        #     gradient = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE_GRADIENT)
        #     node.SetValue(KratosMultiphysics.TEMPERATURE, gradient.norm_2())


        # metric_x = MeshingApplication.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.TEMPERATURE, metric_parameters_hessian)
        # metric_x.Execute()

        # metric_x = MeshingApplication.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.DISTANCE_GRADIENT_X, metric_parameters_hessian)
        # metric_x.Execute()
        # metric_y = MeshingApplication.ComputeHessianSolMetricProcess(self.main_model_part, KratosMultiphysics.DISTANCE_GRADIENT_Y, metric_parameters_hessian)
        # metric_y.Execute()
        # minimal_size=self.metric_parameters["minimal_size"].GetDouble()
        # self.metric_parameters["sizing_parameters"]["boundary_layer_max_distance"].SetDouble(minimal_size*15.0)

            ##COMPUTE LEVEL SET METRIC



        # gid_output = GiDOutputProcess(self.main_model_part,
        #                             "extendend_distance_"+str(self.step),
        #                             KratosMultiphysics.Parameters("""
        #                                 {
        #                                     "result_file_configuration" : {
        #                                         "gidpost_flags": {
        #                                             "GiDPostMode": "GiD_PostBinary",
        #                                             "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                             "WriteConditionsFlag": "WriteConditions",
        #                                             "MultiFileFlag": "SingleFile"
        #                                         },
        #                                         "nodal_results" : ["DISTANCE", "DISTANCE_GRADIENT"],
        #                                         "nodal_nonhistorical_results": ["METRIC_TENSOR_2D","TEMPERATURE","AUXILIAR_HESSIAN","AUXILIAR_GRADIENT","WATER_PRESSURE"]
        #                                     }
        #                                 }
        #                                 """)
        #                             )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()


        # KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SELECTED, False, self.main_model_part.Elements)
        # KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.TO_SPLIT, False, self.main_model_part.Elements)

        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.metric_parameters)
        metric_process.Execute()


        # bound_layer = 0.1
        # distance_to_te = 1000.0
        # distance_to_le = 1000.0
        # minimal_size=self.metric_parameters["minimal_size"].GetDouble()
        # te_size = 1e-5

        # # base_multiplier = 25
        # for node in self.skin_model_part.GetSubModelPart("TrailingEdgeNode").Nodes:
        #     te_node = node
        #     break
        # for node in self.skin_model_part.GetSubModelPart("LeadingEdgeNode").Nodes:
        #     le_node = node
        #     break
        # for node in self.main_model_part.Nodes:
        #     distance_to_te = math.sqrt((node.X-te_node.X)**2+(node.Y-te_node.Y)**2)
        #     distance_to_le = math.sqrt((node.X-le_node.X)**2+(node.Y-le_node.Y)**2)
        #     final_distance = min(distance_to_le, distance_to_te)
        #     if final_distance <= bound_layer:
        #         # multiplier =  base_multiplier - final_distance * base_multiplier / bound_layer + 1.0
        #         # metric_tensor_2d = node.GetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D)
        #         # node.SetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D, multiplier*metric_tensor_2d)
        #         # final_size = te_size + (minimal_size-te_size)*final_distance/bound_layer
        #         final_size = te_size + (minimal_size-te_size)*(1-(math.cos(final_distance/bound_layer*math.pi)+1)/2)
        #         metric_tensor_2d = KratosMultiphysics.Vector(3, 0.0)
        #         metric_tensor_2d[0] =1/final_size/final_size
        #         metric_tensor_2d[1] =1/final_size/final_size
        #         node.SetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D, metric_tensor_2d)

        # gid_output = GiDOutputProcess(self.main_model_part,
        #                             "metric_tensor_"+str(self.step),
        #                             KratosMultiphysics.Parameters("""
        #                                 {
        #                                     "result_file_configuration" : {
        #                                         "gidpost_flags": {
        #                                             "GiDPostMode": "GiD_PostBinary",
        #                                             "WriteDeformedMeshFlag": "WriteUndeformed",
        #                                             "WriteConditionsFlag": "WriteConditions",
        #                                             "MultiFileFlag": "SingleFile"
        #                                         },
        #                                         "nodal_results" : ["DISTANCE", "DISTANCE_GRADIENT"],
        #                                         "nodal_nonhistorical_results": ["METRIC_TENSOR_2D","TEMPERATURE","WATER_PRESSURE","PRESSURE"]
        #                                     }
        #                                 }
        #                                 """)
        #                             )

        # gid_output.ExecuteInitialize()
        # gid_output.ExecuteBeforeSolutionLoop()
        # gid_output.ExecuteInitializeSolutionStep()
        # gid_output.PrintOutput()
        # gid_output.ExecuteFinalizeSolutionStep()
        # gid_output.ExecuteFinalize()
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Preremesh time: ',time.time()-ini_time)
        ini_time = time.time()
        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"                  : "STANDARD",
            "save_external_files"              : false,
            "initialize_entities"              : false,
            "preserve_flags"              : false,
            "interpolate_nodal_values"              : false,
            "echo_level"                       : 0,
            "advanced_parameters"                  :
            {
                "force_hausdorff_value"              : false,
                "hausdorff_value"                    : 0.1,
                "force_gradation_value"              : false,
                "gradation_value"                    : 0.1
            }
        }
        """)

        mmg_process = MeshingApplication.MmgProcess2D(self.main_model_part, mmg_parameters)
        mmg_process.Execute()

        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Remesh time: ',time.time()-ini_time)

    def _UpdateParameters(self):
        ''' This process updates remeshing parameters in case more than one iteration is performed'''
        previous_size=self.metric_parameters["minimal_size"].GetDouble()
        self.metric_parameters["minimal_size"].SetDouble(previous_size*self.update_coefficient)

    def _ModifyFinalDistance(self):
        ''' This function modifies the distance field to avoid ill defined cuts.
        '''
        ini_time = time.time()
        KratosMultiphysics.FluidDynamicsApplication.DistanceModificationProcess(self.main_model_part,self.distance_modification_parameters).Execute()
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Modify distance time: ',time.time()-ini_time)

    def _CopyAndDeleteDefaultDistance(self):
        ''' This function copies the distance field to an auxiliary distance variable and sets
        to zero the default one.
        '''
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.GEOMETRY_DISTANCE, self.main_model_part.Nodes)
        # KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.DISTANCE, self.main_model_part.Nodes)

    def _DefineWakeModelPart(self):
        ''' This function generates the modelpart of the wake. TODO: make end of the domain user-definable.
        '''

        skin_model_part = self.model["skin"]
        te_node = -1
        for node in skin_model_part.GetSubModelPart("TrailingEdgeNode").Nodes:
            te_node = node
        le_node = -1
        for node in skin_model_part.GetSubModelPart("LeadingEdgeNode").Nodes:
            le_node = node
        print(te_node.Id)
        print(le_node.Id)
        te_weight=0.95
        self.wake_model_part.CreateNewNode(1, (1-te_weight)*le_node.X+te_weight*te_node.X, (1-te_weight)*le_node.Y+te_weight*te_node.Y, 0.0)
        print("WAKE_ORIGIN", (1-te_weight)*le_node.X+te_weight*te_node.X, (1-te_weight)*le_node.Y+te_weight*te_node.Y)
        self.wake_model_part.CreateNewNode(2, te_node.X, te_node.Y, 0.0)
        self.wake_model_part.CreateNewNode(3, 200.0, te_node.Y, 0.0)
        self.wake_model_part.CreateNewElement("Element2D2N", 1, [1,2], KratosMultiphysics.Properties(0))
        self.wake_model_part.CreateNewElement("Element2D2N", 2, [2,3], KratosMultiphysics.Properties(0))

    def _CalculateDiscontinuousDistanceAndComputeWakeMetric(self):
        # KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,self.main_model_part.Nodes)


        ini_time=time.time()
        # KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(self.main_model_part, self.wake_model_part).Execute()
        # for elem in self.main_model_part.Elements:
        #     if elem.Is(KratosMultiphysics.TO_SPLIT):
        #         for node in elem.GetNodes():
        #             tensor = KratosMultiphysics.Vector(3, 0.0)
        #             tensor[0] = 1/self.target_h_wake/self.target_h_wake
        #             tensor[1] = 1/self.target_h_wake/self.target_h_wake
        #             node.SetValue(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D, tensor)

        # KratosMultiphysics.FindGlobalNodalElementalNeighboursProcess(self.main_model_part).Execute()
        # CompressiblePotentialFlow.PotentialFlowUtilities.ComputeWakeMetrics(self.main_model_part, self.target_h_wake)
        # KratosMultiphysics.VariableUtils().SaveNonHistoricalVariable(CompressiblePotentialFlow.POTENTIAL_METRIC_TENSOR_3D,KratosMultiphysics.MeshingApplication.METRIC_TENSOR_3D, self.main_model_part.Nodes)

        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Wake metric time ',time.time()-ini_time)
    # @staticmethod
    # def _HasSplitNeighbour(elem):
    #     for node in elem.GetNodes():
    #         neigh_elems = node.GetValue(KratosMultiphysics.NEIGHBOUR_ELEMENTS)
    #         for neigh_elem in neigh_elems:
    #             if neigh_elem.Is(KratosMultiphysics.TO_SPLIT):
    #                 return True
    #     return False