import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
import time

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LevelSetRemeshingProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class LevelSetRemeshingProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name": "insert_model_part",
                "skin_model_part_name": "insert_skin_model_part",
                "maximum_iterations": 1,
                "update_coefficient": 0.5,
                "remeshing_flag": false,
                "ray_casting_tolerance": 1e-9,
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

        self.moving_parameters = settings["moving_parameters"]
        # Synchronizing parameters for the wake process
        if self.moving_parameters.Has("rotation_point"):
            self.main_model_part.ProcessInfo.SetValue(CompressiblePotentialFlow.WAKE_ORIGIN, self.moving_parameters["rotation_point"].GetVector())
        else:
            self.main_model_part.ProcessInfo.SetValue(CompressiblePotentialFlow.WAKE_ORIGIN, self.moving_parameters["origin"].GetVector())
        self.main_model_part.ProcessInfo.SetValue(CompressiblePotentialFlow.ROTATION_ANGLE, self.moving_parameters["rotation_angle"].GetDouble())

        self.metric_parameters = settings["metric_parameters"]
        self.distance_modification_parameters = settings["distance_modification_parameters"]
        self.ray_casting_tolerance = settings["ray_casting_tolerance"].GetDouble()

    def ExecuteInitialize(self):
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Executing Initialize Geometry')
        self._InitializeSkinModelPart()
        self._CalculateDistance()

        ini_time=time.time()
        while self.step < self.max_iter and self.do_remeshing:
            self.step += 1
            KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','##### Executing refinement #', self.step, ' #####')
            self._ExtendDistance()
            self._RefineMesh()
            self._CalculateDistance()
            self._UpdateParameters()
        self._ModifyFinalDistance()
        self._CopyAndDeleteDefaultDistance()
        KratosMultiphysics.Logger.PrintInfo('LevelSetRemeshing','Elapsed time: ',time.time()-ini_time)

    def _InitializeSkinModelPart(self):
        ''' This function loads and moves the skin_model_part in the main_model_part to the desired initial point (origin).
            It also rotates the skin model part around the origin point according to the rotation_angle'''
        self.skin_model_part=self.model.CreateModelPart("skin")

        ini_time=time.time()
        # Reading skin model part
        KratosMultiphysics.ModelPartIO(self.skin_model_part_name).ReadModelPart(self.skin_model_part)
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
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariableToZero(KratosMultiphysics.MeshingApplication.METRIC_TENSOR_2D,self.main_model_part.Nodes)

        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.metric_parameters)
        metric_process.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"                  : "STANDARD",
            "save_external_files"              : false,
            "initialize_entities"              : false,
            "echo_level"                       : 0
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
        KratosMultiphysics.VariableUtils().SetHistoricalVariableToZero(KratosMultiphysics.DISTANCE, self.main_model_part.Nodes)