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
                "geometry_parameter": 0.0,
                "remeshing_flag": false,
                "isosurface_flag": false,
                "initial_point": [0,0],
                "node_id": -1,
                "epsilon": 1e-9,
                "metric_parameters":  {
                    "minimal_size"                         : 5e-3, 
                    "maximal_size"                         : 1.0, 
                    "custom_settings":
                    {
                        "ratio"                     : 1.0,
                        "max_x"                     : 1.0,
                        "min_x"                     : 0.0
                    },
                    "sizing_parameters": {   
                        "reference_variable_name"               : "DISTANCE",
                        "boundary_layer_max_distance"           : 1.0, 
                        "interpolation"                         : "Constant"
                    },
                    "enforce_current"                      : true, 
                    "anisotropy_remeshing"                 : false, 
                    "anisotropy_parameters": {   
                        "reference_variable_name"          : "DISTANCE",
                        "hmin_over_hmax_anisotropic_ratio"      : 0.5, 
                        "boundary_layer_max_distance"           : 1, 
                        "interpolation"                         : "Linear"
                    }
                }
            }  """ );
        settings.ValidateAndAssignDefaults(default_parameters)
        self.model=Model
        self.main_model_part = Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.geometry_parameter = settings["geometry_parameter"].GetDouble();       
        self.do_remeshing = settings["remeshing_flag"].GetBool();   
        self.isosurface_flag = settings["isosurface_flag"].GetBool();   
        self.initial_point = KratosMultiphysics.Vector(2)
        self.initial_point[0] = settings["initial_point"][0].GetDouble()
        self.initial_point[1] = settings["initial_point"][1].GetDouble()
        self.skin_model_part=self.model.CreateModelPart("skin")
        self.skin_model_part_name=settings["skin_model_part_name"].GetString()
        self.boundary_model_part = self.main_model_part.CreateSubModelPart("boundary_model_part")
        KratosMultiphysics.ModelPartIO(self.skin_model_part_name).ReadModelPart(self.skin_model_part)

        self.node_id=settings["node_id"].GetInt()
        self.epsilon=settings["epsilon"].GetDouble()

        self.MetricParameters = settings["metric_parameters"]
        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(3)
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        for node in self.main_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_2D, ZeroVector)

        import python_linear_solver_factory #Linear solver for variational distance process
        linear_solver_settings=KratosMultiphysics.Parameters("""
        {
            "solver_type": "amgcl",
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
        
        self.linear_solver = python_linear_solver_factory.ConstructSolver(linear_solver_settings)

    def Execute(self):
        print("Executing Initialize Geometry")
        self.InitializeSkinModelPart()
        self.CalculateDistance()
        # self.PerturbateDistanceNumericalGradient()

        if (self.do_remeshing):
            if self.isosurface_flag:
                # self.RefineMesh()
                self.FindIsosurface()
                self.CalculateDistance()
            else:
                # print("Executing first refinement ")
                # import time
                # ini_time=time.time()
                # initial_min=self.MetricParameters["minimal_size"].GetDouble()
                # initial_distance = self.MetricParameters["sizing_parameters"]["boundary_layer_max_distance"].GetDouble()
                # self.MetricParameters["minimal_size"].SetDouble(5e-3)
                # self.MetricParameters["sizing_parameters"]["boundary_layer_max_distance"].SetDouble(1.0)
                # self.RefineMesh()
                # print("Elapsed time: ",time.time()-ini_time)
                # self.CalculateDistance()

                # print("Executing second refinement ")
                # self.MetricParameters["minimal_size"].SetDouble(initial_min)
                # self.MetricParameters["sizing_parameters"]["boundary_layer_max_distance"].SetDouble(initial_distance)
                # print("Executing second refinement ")
                # ini_time=time.time()
                # self.RefineMesh()
                # print("Elapsed time: ",time.time()-ini_time)
                # self.CalculateDistance()
                for i in range (0,3):
                    print("Executing refinement #",i+1)
                    import time
                    ini_time=time.time()
                    self.RefineMesh()
                    self.CalculateDistance()
                    print("Elapsed time: ",time.time()-ini_time)




        self.ApplyFlags()
        self.ComputeKuttaNodeAndElement()
        print("Level Set geometry initialized")
        
        
    def ExecuteInitialize(self):
        self.Execute()

    def InitializeSkinModelPart(self):
        if self.skin_model_part_name=='naca0012':
            angle=math.radians(-self.geometry_parameter)
            self.origin=[0.25+self.initial_point[0],0+self.initial_point[1]] #to be defined from skin model part 
            for node in self.skin_model_part.Nodes:
                node.X=self.initial_point[0]+node.X+1e-5
                node.Y=self.initial_point[1]+node.Y+1e-5
            RotateModelPart(self.origin,angle,self.skin_model_part)
        elif self.skin_model_part_name=='circle':
            for node in self.skin_model_part.Nodes:
                node.X=node.X+1e-5
                node.Y=node.Y+1e-5
        elif self.skin_model_part_name=='ellipse' or self.skin_model_part_name=='ellipse_fine':
            angle=math.radians(-self.geometry_parameter)   
            self.origin=[self.initial_point[0],self.initial_point[1]] #to be defined from skin model part  
            for node in self.skin_model_part.Nodes:
                node.X=self.initial_point[0]+node.X+1e-5
                node.Y=self.initial_point[1]+node.Y+1e-5
            RotateModelPart(self.origin,angle,self.skin_model_part)
            # a=1
            # b=a/4
            # for node in self.main_model_part.Nodes:
            #     distance = ((node.X*math.cos(angle)+node.Y*math.sin(angle))/a)**2 + ((node.X*math.sin(angle)+node.Y*math.cos(angle))/b)**2 -1
            #     if distance==0:
            #         distance=1e-6
            #     node.SetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE,distance)


    def CalculateDistance(self):
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.main_model_part, self.skin_model_part).Execute()
        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.LEVEL_SET_DISTANCE, self.main_model_part.Nodes)
        return

    def FindIsosurface(self):
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
        import mmg_process
        mmg_parameters = KratosMultiphysics.Parameters("""
        {  
            "Parameters":{
                "model_part_name"                  : "Parts_Fluid",
                "save_external_files"              : false,
                "initialize_entities"              : false,
                "automatic_remesh"                 : false,
                "echo_level"                       : 1  ,      
                "discretization_type"                  : "Standard",
                "minimal_size"                         : 0.001, 
                "maximal_size"                         : 1.0,
                "sizing_parameters":
                {
                    "reference_variable_name"          : "DISTANCE",
                    "boundary_layer_max_distance"      : 1.0,
                    "interpolation"                    : "Linear"
                },
                "enforce_current"                      : false, 
                "anisotropy_remeshing"                 : false, 
                "anisotropy_parameters": {   
                    "reference_variable_name"          : "DISTANCE",
                    "hmin_over_hmax_anisotropic_ratio"      : 0.5, 
                    "boundary_layer_max_distance"           : 1.0, 
                    "interpolation"                         : "Linear"
                },
                "isosurface_parameters"                    : {
                        "remove_regions"                   : true
                }
            }
            
        }""")

        mmg=mmg_process.Factory(mmg_parameters,self.model)
        mmg.ExecuteInitialize()
        mmg.ExecuteInitializeSolutionStep()

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

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()
        self.FindXMaxAndMinSkin()
        self.MetricParameters["custom_settings"]["max_x"].SetDouble(self.max_x)
        self.MetricParameters["custom_settings"]["max_x"].SetDouble(self.min_x)

        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.MetricParameters)
        metric_process.Execute()

        self.PrintOutput('metric_output')

        mmg_parameters = KratosMultiphysics.Parameters("""
        {
            "discretization_type"                  : "STANDARD",
            "save_external_files"              : false,
            "initialize_entities"              : false,
            "echo_level"                       : 0,
            "isosurface_parameters"                    : {
                    "remove_regions"                   : false
                }
        }
        """)

        # We create the remeshing utility
        # mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
        mmg_process = MeshingApplication.MmgProcess2D(self.main_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()

        self.PrintOutput('remeshed_output')

        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.LEVEL_SET_DISTANCE, self.main_model_part.Nodes)

    def ApplyFlags(self):
        max_x=-1e10
        for element in self.main_model_part.Elements:
            IsPositive=False
            IsNegative=False
            element.Set(KratosMultiphysics.ACTIVE,True)
            for node in element.GetNodes():
                distance=node.GetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE)
                if distance>0:
                    IsPositive=True
                else:
                    IsNegative=True
            if IsPositive and IsNegative:
                element.Set(KratosMultiphysics.BOUNDARY,True)
                self.boundary_model_part.Elements.append(element)
            elif IsPositive:
                element.Set(KratosMultiphysics.FLUID,True)
            else:
                element.Set(KratosMultiphysics.FLUID,False)
                element.Set(KratosMultiphysics.ACTIVE,False)
                element.Set(KratosMultiphysics.BOUNDARY,False)
                for node in element.GetNodes():
                    if node.X>max_x:
                        max_x=node.X
                        self.max_node=node
    def ComputeKuttaNodeAndElement(self):    
        x0 =  self.origin[0]
        y0 =  self.origin[1]
        direction = KratosMultiphysics.Vector(3)
        angle=math.radians(self.geometry_parameter)
        direction[0]=-math.cos(angle)
        direction[1]=math.sin(angle)
        n = KratosMultiphysics.Vector(3)
        xn = KratosMultiphysics.Vector(3)
        n[0] = -direction[1]
        n[1] = direction[0]
        n[2] = 0.0
        max_x_center = -1e10
        for elem in self.boundary_model_part.Elements:         
            distances = KratosMultiphysics.Vector(len(elem.GetNodes()))
            '''with positive epsilon'''
            npos = 0
            nneg = 0
            for elnode in elem.GetNodes():
                xn[0] = elnode.X - x0
                xn[1] = elnode.Y - y0
                xn[2] = 0.0
                d =  xn[0]*n[0] + xn[1]*n[1]
                if(d < 0):
                    nneg += 1
                else:
                    npos += 1
                elnode.SetValue(KratosMultiphysics.TEMPERATURE,d)
            if(nneg>0 and npos>0):
                center_X=elem.GetGeometry().Center().X                
                if (center_X>self.origin[0]):
                    elem.Set(KratosMultiphysics.THERMAL,True)
                    elem.Set(KratosMultiphysics.ACTIVE,False)
                    elem.Set(KratosMultiphysics.BOUNDARY,False)
                    if (center_X > max_x_center):
                        max_x_center = center_X
                        max_elem = elem

        max_elem.Set(KratosMultiphysics.ACTIVE,False)
        max_elem.Set(KratosMultiphysics.BOUNDARY,False)
        max_x_node=-1e10
        for elnode in max_elem.GetNodes():
            if elnode.X>max_x_node:
                max_node=elnode
                max_x_node=elnode.X

        self.main_model_part.CreateSubModelPart('KuttaLS').AddNode(max_node,0)
        for node in self.main_model_part.GetSubModelPart('KuttaLS').Nodes:
            node.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.TRAILING_EDGE,True)
    def PrintOutput(self,filename):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(self.main_model_part,
                                    filename,
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results" : ["LEVEL_SET_DISTANCE"],
                                                "nodal_nonhistorical_results": ["METRIC_TENSOR_2D"]

                                            }
                                        }
                                        """)
                                    )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()
        
    def FindXMaxAndMinSkin(self):
        self.max_x=-1e10
        self.min_x=1e10
        for node in self.skin_model_part.Nodes:
            if node.X > self.max_x:
                self.max_x=node.X
            if node.X < self.min_x:
                self.min_x=node.X
        print("MIN AND MAX X", self.min_x, self.max_x)

