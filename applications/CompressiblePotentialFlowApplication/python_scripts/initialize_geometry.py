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
                    "enforce_current"                      : true, 
                    "anisotropy_remeshing"                 : true, 
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
        self.PerturbateDistanceNumericalGradient()

        if (self.do_remeshing):
            if self.isosurface_flag:
                # self.RefineMesh()
                self.FindIsosurface()
                self.CalculateDistance()
            else:
                print("Executing first refinement ")
                import time
                ini_time=time.time()
                self.RefineMesh()
                print("Elapsed time: ",time.time()-ini_time)
                self.CalculateDistance()
                # self.MetricParameters["enforce_current"].SetBool(True)
                # self.MetricParameters["minimal_size"].SetDouble(5e-4)
                # print("Executing second refinement ")
                # self.RefineMesh()
                # self.CalculateDistance()

        self.ApplyFlags()
        print("Level Set geometry initialized")
        
        
    def ExecuteInitialize(self):
        self.Execute()

    def InitializeSkinModelPart(self):
        if self.skin_model_part_name=='naca0012':
            angle=math.radians(-self.geometry_parameter)
            self.origin=[0.25+self.initial_point[0],0+self.initial_point[1]] #to be defined from skin model part          
            for node in self.skin_model_part.Nodes:
                node.X=self.initial_point[0]+node.X+1e-5
                node.Y=node.Y+1e-5
            RotateModelPart(self.origin,angle,self.skin_model_part)
        elif self.skin_model_part_name=='circle':
            for node in self.skin_model_part.Nodes:
                node.X=node.X+1e-5
                node.Y=node.Y+1e-5
        elif self.skin_model_part_name=='ellipse':
            angle=math.radians(-self.geometry_parameter)   
            self.origin=[0.0,0.0] #to be defined from skin model part  
            for node in self.skin_model_part.Nodes:
                node.X=node.X+1e-5
                node.Y=node.Y+1e-5
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
    def PerturbateDistanceNumericalGradient(self):
        if self.node_id>0.0:
            print("IM IN")
            print(self.node_id)
            for node in self.main_model_part.Nodes:
                if node.Id==self.node_id:
                    distance=node.GetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE)
                    node.SetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE,distance+self.epsilon)
                    break
            # distance=self.main_model_part.GetNodes(self.node_id).GetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE)
            print("distance value about to be perturbated:",distance)
            # self.main_model_part.GetNodes(node_id).SetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE,distance+self.epsilon)
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
                "echo_level"                       : 0  ,      
                "discretization_type"                  : "ISOSURFACE",
                "minimal_size"                         : 0.001, 
                "enforce_current"                      : true, 
                "anisotropy_remeshing"                 : true, 
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
        mmg.Execute()

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

        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.main_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.MetricParameters)
        metric_process.Execute()

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
        # print(self.main_model_part)
        
        mmg_process.Execute()

        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(self.main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISTANCE"]
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
                for node in element.GetNodes():
                    if node.X>max_x:
                        max_x=node.X
                        max_node=node
            elif IsPositive:
                element.Set(KratosMultiphysics.FLUID,True)
            else:
                element.Set(KratosMultiphysics.FLUID,False)
                element.Set(KratosMultiphysics.ACTIVE,False)
                # for node in element.GetNodes():
                #     if node.X>max_x:
                #         max_x=node.X
                #         max_node=node
        self.main_model_part.CreateSubModelPart('KuttaLS').AddNode(max_node,0)
        for node in self.main_model_part.GetSubModelPart('KuttaLS').Nodes:
            node.Set(KratosMultiphysics.STRUCTURE,True)
##################
        # max_node=self.main_model_part.GetNode(7069,0)
##################
        self.aux_model_part = self.main_model_part.CreateSubModelPart("aux_model_part")
        for elem in self.boundary_model_part.Elements:
            center_X=elem.GetGeometry().Center().X
            touching_kutta=False    
            for node in elem.GetNodes():
                if node.Id == max_node.Id:
                    touching_kutta=True
            if center_X<max_x and touching_kutta:
                elem.Set(KratosMultiphysics.ACTIVE,False)
                self.aux_model_part.Elements.append(elem)
        node_list=[]
        for elem in self.aux_model_part.Elements:
            for node in elem.GetNodes():
                if (not node in node_list):
                    if (not node.Id==max_node.Id):
                        node_list.append(node)
                elif (not node.Id==max_node.Id):
                    shared_node=node
        self.aux_model_part2 = self.main_model_part.CreateSubModelPart("aux_model_part2")
        for elem in self.boundary_model_part.Elements:
            # center_X=elem.GetGeometry().Center().X
            touching_shared=False    
            for node in elem.GetNodes():
                if node.Id == shared_node.Id:
                    touching_shared=True
            if touching_shared:
                elem.Set(KratosMultiphysics.ACTIVE,False)
        #         self.aux_model_part2.Elements.append(elem)
        # node_list=[]
        # for element in self.aux_model_part2.Elements:
        #     for node in element.GetNodes():
        #         if not node.Id in node_list:
        #             node_list.append(node.Id)
        # for element in self.boundary_model_part.Elements:
        #     counter=0
        #     for node in element.GetNodes():
        #         if node.Id in node_list:
        #             counter +=1
        #     if counter==2:
        #         element.Set(KratosMultiphysics.ACTIVE,False)

        # self.main_model_part.GetElement(142868,0).Set(KratosMultiphysics.ACTIVE,False)
        # self.main_model_part.GetElement(473307,0).Set(KratosMultiphysics.ACTIVE,False)

        
        
        # x0 = max_node.X
        # y0 = max_node.Y
        # direction = KratosMultiphysics.Vector(3)
        # angle=math.radians(self.geometry_parameter)
        # direction[0]=-math.cos(angle)
        # direction[1]=math.sin(angle)
        # n = KratosMultiphysics.Vector(3)
        # xn = KratosMultiphysics.Vector(3)
        # n[0] = -direction[1]
        # n[1] = direction[0]
        # n[2] = 0.0
        # for elem in self.boundary_model_part.Elements:         
        #     distances = KratosMultiphysics.Vector(len(elem.GetNodes()))
        #     '''with positive epsilon'''
        #     counter = 0
        #     pos_nodes=[]
        #     npos = 0
        #     nneg = 0
        #     for elnode in elem.GetNodes():
        #         xn[0] = elnode.X - x0
        #         xn[1] = elnode.Y - y0
        #         xn[2] = 0.0
        #         d =  xn[0]*n[0] + xn[1]*n[1]
        #         if(d < 0):
        #             nneg += 1
        #             pos_nodes.append(elnode.Id)
        #         else:
        #             npos += 1

        #         if(abs(d) < 1e-6):
        #             if d >= 0.0:
        #                 d = 1e-6
        #             else:
        #                 d = -1e-6
        #         counter += 1

        #     if(nneg>0 and npos>0):
        #         center_X=elem.GetGeometry().Center().X
        #         if (center_X>self.origin[0]):
        #             elem.Set(KratosMultiphysics.ACTIVE,False)
        #             elem.Set(KratosMultiphysics.THERMAL,True)