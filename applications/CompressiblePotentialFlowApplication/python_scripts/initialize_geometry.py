import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
import sys
import pprint


def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InitializeGeometryProcess(Model, settings["Parameters"])
def RotateModelPart(origin, angle, model_part):
    ox,oy=origin        
    for node in model_part.Nodes:
        node.X = ox+math.cos(angle)*(node.X - ox)-math.sin(angle)*(node.Y - oy)
        node.Y = oy+math.sin(angle)*(node.X - ox)+math.cos(angle)*(node.Y - oy)   
## All the processes python should be derived from "Process"
class InitializeGeometryProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        print("Initialize Geometry process")
        default_parameters = KratosMultiphysics.Parameters( """
            {   
                "model_part_name": "insert_model_part",
                "skin_model_part_name": "insert_skin_model_part",
                "geometry_parameter": 1.0,
                "remeshing_flag": false,
                "initial_point": [0,0],
                "metric_parameters":  {
                    "minimal_size"                         : 1e-2, 
                    "enforce_current"                      : true, 
                    "anisotropy_remeshing"                 : true, 
                    "anisotropy_parameters": {   
                        "reference_variable_name"          : "DISTANCE",
                        "hmin_over_hmax_anisotropic_ratio"      : 0.5, 
                        "boundary_layer_max_distance"           : 0.5, 
                        "interpolation"                         : "Linear"
                    }
                }
            }  """ );

        settings.ValidateAndAssignDefaults(default_parameters)
        self.model=Model
        self.main_model_part = Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.geometry_parameter = settings["geometry_parameter"].GetDouble();       
        self.do_remeshing = settings["remeshing_flag"].GetBool();   
        self.initial_point = KratosMultiphysics.Vector(2)
        self.initial_point[0] = settings["initial_point"][0].GetDouble()
        self.initial_point[1] = settings["initial_point"][1].GetDouble()
        self.skin_model_part=self.model.CreateModelPart("skin")
        self.skin_model_part_name=settings["skin_model_part_name"].GetString()
        KratosMultiphysics.ModelPartIO(self.skin_model_part_name).ReadModelPart(self.skin_model_part)
        
        self.MetricParameters = settings["metric_parameters"]
        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(3)
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        for node in self.main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_2D, ZeroVector)

        import linear_solver_factory #Linear solver for variational distance process
        linear_solver_settings=KratosMultiphysics.Parameters("""
        {
            "solver_type": "AMGCL",
            "max_iteration": 400,
            "gmres_krylov_space_dimension": 100,
            "smoother_type":"ilu0",
            "coarsening_type":"ruge_stuben",
            "coarse_enough" : 5000,
            "krylov_type": "lgmres",
            "tolerance": 1e-9,
            "verbosity": 0,
            "scaling": false
        }""")
        
        self.linear_solver = linear_solver_factory.ConstructSolver(linear_solver_settings)

    def Execute(self):
        print("Executing Initialize Geometry")
        self.InitializeSkinModelPart()
        self.CalculateDistance()
        if (self.do_remeshing):
            self.RefineMesh()
            self.CalculateDistance()

        self.ApplyFlags()
        print("Level Set geometry initialized")
        
        
    def ExecuteInitialize(self):
        self.Execute()

    def InitializeSkinModelPart(self):
        if self.skin_model_part_name=='naca0012':
            angle=math.radians(-self.geometry_parameter)
            origin=[0.25+self.initial_point[0],0+self.initial_point[1]] #to be defined from skin model part
            RotateModelPart(origin,angle,self.skin_model_part)
            for node in self.skin_model_part.Nodes:
                node.X=node.X+1e-5
                node.Y=node.Y+1e-5
        elif self.skin_model_part_name=='circle':
            for node in self.skin_model_part.Nodes:
                node.X=self.geometry_parameter*node.X+1e-5
                node.Y=self.geometry_parameter*node.Y+1e-5


    def CalculateDistance(self):
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.main_model_part, self.skin_model_part).Execute()
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.LEVEL_SET_DISTANCE, self.main_model_part.Nodes)

    def RefineMesh(self):
        KratosMultiphysics.VariableUtils().CopyScalarVar(CompressiblePotentialFlow.LEVEL_SET_DISTANCE,KratosMultiphysics.DISTANCE, self.main_model_part.Nodes)
        
        # Construct the variational distance calculation process
        maximum_iterations = 2 #TODO: Make this user-definable
        if self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                self.main_model_part,
                self.linear_solver,
                maximum_iterations)
        else:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                self.main_model_part,
                self.linear_solver,
                maximum_iterations)
        variational_distance_process.Execute()     

        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.main_model_part)
        find_nodal_h.Execute()

        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.MetricParameters)
        metric_process.Execute()

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "save_external_files"              : true,
            "initialize_entities"              : false,
            "echo_level"                       : 0
        }
        """)

        # We create the remeshing utility
        # mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess2D(self.main_model_part, mmg_parameters)

        # We remesh
        # print(self.main_model_part)
        
        mmg_process.Execute()
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.LEVEL_SET_DISTANCE, self.main_model_part.Nodes)

    def ApplyFlags(self):

        for element in self.main_model_part.Elements:
            IsPositive=False
            IsNegative=False
            for node in element.GetNodes():
                distance=node.GetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE)
                if distance>0:
                    IsPositive=True
                else:
                    IsNegative=True
            if IsPositive and IsNegative:
                element.Set(KratosMultiphysics.BOUNDARY,True)
                for node in element.GetNodes():
                    node.Set(KratosMultiphysics.VISITED,True)
            elif IsPositive:
                element.Set(KratosMultiphysics.FLUID,True)
                for node in element.GetNodes():
                    node.Set(KratosMultiphysics.VISITED,True)
            else:
                element.Set(KratosMultiphysics.FLUID,False)
                element.Set(KratosMultiphysics.ACTIVE,False)
        # for node in self.main_model_part.Nodes:                    
        #     if node.IsNot(KratosMultiphysics.VISITED):
        #         node.Set(KratosMultiphysics.VISITED,True)
        #         # node.SetSolutionStepValue(CompressiblePotentialFlow.POSITIVE_POTENTIAL)
        #         # node.SetSolutionStepValue(CompressiblePotentialFlow.NEGATIVE_POTENTIAL,0)
        #         node.Fix(CompressiblePotentialFlow.POSITIVE_POTENTIAL)
        #         node.Fix(CompressiblePotentialFlow.NEGATIVE_POTENTIAL)

        
        
