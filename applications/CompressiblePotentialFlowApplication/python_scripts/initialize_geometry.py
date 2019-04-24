import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow
import KratosMultiphysics.MeshingApplication as MeshingApplication
import math
import time


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
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','Initialize Geometry process')
        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name": "insert_model_part",
                "skin_model_part_name": "insert_skin_model_part",
                "geometry_parameter": 0.0,
                "maximum_iterations": 3,
                "remeshing_flag": false,
                "initial_point": [0,0],
                "metric_parameters":  {
                    "minimal_size"                         : 5e-3,
                    "maximal_size"                         : 1.0,
                    "custom_settings":
                    {
                        "ratio"                         : 0.5,
                        "chord_ratio"                   : 4.0,
                        "max_node"                      : [0.0,0.0,0.0],
                        "min_node"                      : [0.0,0.0,0.0]
                    },
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
                }
            }  """ );
        settings.ValidateAndAssignDefaults(default_parameters)
        self.model=Model
        self.main_model_part = Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.geometry_parameter = settings["geometry_parameter"].GetDouble()
        self.do_remeshing = settings["remeshing_flag"].GetBool()
        self.initial_point = KratosMultiphysics.Vector(2)
        self.initial_point[0] = settings["initial_point"][0].GetDouble()
        self.initial_point[1] = settings["initial_point"][1].GetDouble()
        self.skin_model_part=self.model.CreateModelPart("skin")
        self.skin_model_part_name=settings["skin_model_part_name"].GetString()
        self.boundary_model_part = self.main_model_part.CreateSubModelPart("boundary_model_part")
        KratosMultiphysics.ModelPartIO(self.skin_model_part_name).ReadModelPart(self.skin_model_part)

        '''Loop parameters'''
        self.step = 0
        self.max_iter = settings["maximum_iterations"].GetInt()

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
        KratosMultiphysics.Logger.PrintInfo('Executing Initialize Geometry')
        self.InitializeSkinModelPart()
        self.CalculateDistance()

        ini_time=time.time()
        while self.step < self.max_iter and self.do_remeshing:

            self.step += 1
            previous_ratio=self.MetricParameters["custom_settings"]["ratio"].GetDouble()
            self.MetricParameters["custom_settings"]["ratio"].SetDouble(previous_ratio*0.5)

            KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','##### Executing refinement #', self.step, ' #####')
            self.ExtendDistance()
            self.RefineMesh()
            self.CalculateDistance()
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','Elapsed time: ',time.time()-ini_time)

        KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.LEVEL_SET_DISTANCE, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DISTANCE, 0.0, self.main_model_part.Nodes)
        self.ApplyFlags()
        self.ComputeKuttaNodeAndElement()
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','Level Set geometry initialized')


    def ExecuteInitialize(self):
        self.Execute()

    def InitializeSkinModelPart(self):
        ini_time=time.time()
        if self.skin_model_part_name=='naca0012':
            angle=math.radians(-self.geometry_parameter)
            self.origin=[0.25+self.initial_point[0],0+self.initial_point[1]] #to be defined from skin model part
            for node in self.skin_model_part.Nodes:
                node.X=self.initial_point[0]+node.X
                node.Y=self.initial_point[1]+node.Y
            RotateModelPart(self.origin,angle,self.skin_model_part)
        elif self.skin_model_part_name=='circle':
            for node in self.skin_model_part.Nodes:
                node.X=node.X+1e-5
                node.Y=node.Y+1e-5
        elif self.skin_model_part_name=='ellipse' or self.skin_model_part_name=='ellipse_fine':
            angle=math.radians(-self.geometry_parameter)
            self.origin=[self.initial_point[0],self.initial_point[1]] #to be defined from skin model part
            for node in self.skin_model_part.Nodes:
                node.X=self.initial_point[0]+node.X
                node.Y=self.initial_point[1]+node.Y
            RotateModelPart(self.origin,angle,self.skin_model_part)
        self.FindXMaxAndMinSkin()
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','InitializeSkin Time: ',time.time()-ini_time)

    def CalculateDistance(self):
        ini_time=time.time()
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.main_model_part, self.skin_model_part).Execute()
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','CalculateDistance Time: ',time.time()-ini_time)

    def ExtendDistance(self):
        ini_time=time.time()
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
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','Variational distance process Time: ',time.time()-ini_time)


    def RefineMesh(self):
        ini_time=time.time()
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(self.main_model_part)
        find_nodal_h.Execute()


        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.MetricParameters)
        metric_process.Execute()

        self.PrintOutput('metric_output'+str(self.step))

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

        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','Remesh Time: ',time.time()-ini_time)

        self.PrintOutput('remeshed_output'+str(self.step))

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
                # for node in element.GetNodes():
                #     if node.X>max_x:
                #         max_x=node.X
                #         self.max_node=node
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
        # x0 =  self.origin[0]
        # y0 =  self.origin[1]
        # direction = KratosMultiphysics.Vector(3)
        # angle=math.radians(self.geometry_parameter)
        # direction[0]=-math.cos(angle)
        # direction[1]=math.sin(angle)
        # n = KratosMultiphysics.Vector(3)
        # xn = KratosMultiphysics.Vector(3)
        # n[0] = -direction[1]
        # n[1] = direction[0]
        # n[2] = 0.0
        # max_x_center = -1e10
        # for elem in self.boundary_model_part.Elements:
        #     distances = KratosMultiphysics.Vector(len(elem.GetNodes()))
        #     '''with positive epsilon'''
        #     npos = 0
        #     nneg = 0
        #     for elnode in elem.GetNodes():
        #         xn[0] = elnode.X - x0
        #         xn[1] = elnode.Y - y0
        #         xn[2] = 0.0
        #         d =  xn[0]*n[0] + xn[1]*n[1]
        #         if(d < 0):
        #             nneg += 1
        #         else:
        #             npos += 1
        #         elnode.SetValue(KratosMultiphysics.TEMPERATURE,d)
        #     if(nneg>0 and npos>0):
        #         center_X=elem.GetGeometry().Center().X
        #         if (center_X>self.origin[0]):
        #             elem.Set(KratosMultiphysics.THERMAL,True)
        #             elem.Set(KratosMultiphysics.ACTIVE,False)
        #             elem.Set(KratosMultiphysics.BOUNDARY,False)
        #             if (center_X > max_x_center):
        #                 max_x_center = center_X
        #                 max_elem = elem

        # max_elem.Set(KratosMultiphysics.ACTIVE,False)
        # max_elem.Set(KratosMultiphysics.BOUNDARY,False)
        # max_x_node=-1e10
        # for elnode in max_elem.GetNodes():
        #     if elnode.X>max_x_node:
        #         max_node=elnode
        #         max_x_node=elnode.X
        # self.main_model_part.CreateSubModelPart('KuttaLS').AddNode(max_node,0)
        # self.main_model_part.CreateSubModelPart('KuttaLS').AddNode(self.max_node,0)



        max_x_node=-1e10
        deactivated_model_part=self.main_model_part.CreateSubModelPart('deactivated')
        for element in self.boundary_model_part.Elements:
            for elnode in element.GetNodes():
                if elnode.Id == self.max_node.Id:
                    self.DeactivateActive(element)
                    deactivated_model_part.Elements.append(element)
        max_x_center=-1e10
        for element in deactivated_model_part.Elements:
            n_center = 0
            center_X=element.GetGeometry().Center().X
            if center_X>max_x_center:
                max_x_center=center_X
                max_elem=element

        kutta_model_part=self.main_model_part.CreateSubModelPart('KuttaLS')
        for elnode in max_elem.GetNodes():
            if elnode.X>max_elem.GetGeometry().Center().X:
                n_center+=1
                kutta_model_part.AddNode(elnode,0)
        for elnode in max_elem.GetNodes():
            if elnode.X>max_x_node:
                max_x_node=elnode.X
                max_node=elnode
                # max_elem=element
                # max_n_center=n_center

        kutta_model_part.AddNode(max_node,0)

        for element in self.boundary_model_part.Elements:
            center_X=element.GetGeometry().Center().X
            if center_X>self.max_node.X:
                self.DeactivateBoundary(element)

        for node in self.main_model_part.GetSubModelPart('KuttaLS').Nodes:
            node.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.TRAILING_EDGE,True)
    def DeactivateActive(self,elem):
        elem.Set(KratosMultiphysics.ACTIVE,False)
        elem.Set(KratosMultiphysics.BOUNDARY,False)
    def DeactivateBoundary(self,elem):
        elem.Set(KratosMultiphysics.BOUNDARY,False)

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
                                                "nodal_results" : ["DISTANCE"],
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
        max_x=-1e10
        min_x=1e10
        for node in self.skin_model_part.Nodes:
            if node.X > max_x:
                max_x=node.X
                max_node=node
            if node.X < min_x:
                min_x=node.X
                min_node=node
        min_node_vec=KratosMultiphysics.Vector(3)
        min_node_vec[0]=min_node.X
        min_node_vec[1]=min_node.Y
        min_node_vec[2]=min_node.Z
        max_node_vec=KratosMultiphysics.Vector(3)
        max_node_vec[0]=max_node.X
        max_node_vec[1]=max_node.Y
        max_node_vec[2]=max_node.Z
        self.MetricParameters["custom_settings"]["max_node"].SetVector(max_node_vec)
        self.MetricParameters["custom_settings"]["min_node"].SetVector(min_node_vec)
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','MIN SKIN NODE', min_node.X, min_node.Y)
        KratosMultiphysics.Logger.PrintInfo('InitializeGeometry','MAX SKIN NODE', max_node.X, max_node.Y)


