import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow
import KratosMultiphysics.MeshingApplication as MeshingApplication

import math

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess(Model, settings["Parameters"])
def RotateModelPart(origin, angle, model_part):
    ox,oy=origin        
    for node in model_part.Nodes:
        node.X = ox+math.cos(angle)*(node.X - ox)-math.sin(angle)*(node.Y - oy)
        node.Y = oy+math.sin(angle)*(node.X - ox)+math.cos(angle)*(node.Y - oy)

class DefineWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        print("Initialize Wake Process")
        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "model_part_name"           : "please specify the model part that contains the kutta nodes",
                "fluid_part_name"           : "MainModelPart",
                "direction"                 : [1.0,0.0,0.0],
                "stl_filename"              : "please specify name of stl file",
                "epsilon"    : 1e-9,
                "geometry_parameter": 0.0
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings) 
        self.geometry_parameter=settings["geometry_parameter"].GetDouble()
        self.direction = KratosMultiphysics.Vector(3)
        self.direction[0] = settings["direction"][0].GetDouble()
        self.direction[1] = settings["direction"][1].GetDouble()
        self.direction[2] = settings["direction"][2].GetDouble()
        dnorm = math.sqrt(self.direction[0]**2 + self.direction[1]**2 + self.direction[2]**2)
        self.direction[0] /= dnorm
        self.direction[1] /= dnorm
        self.direction[2] /= dnorm
        
        self.epsilon = settings["epsilon"].GetDouble()
        
        self.fluid_model_part = Model.GetModelPart(settings["fluid_part_name"].GetString()).GetRootModelPart()
        self.model=Model
        self.wake_model_part_name=settings["model_part_name"].GetString()
        self.wake_line_model_part=Model.CreateModelPart("wake")
        angle=math.radians(-self.geometry_parameter)
        for node in self.wake_line_model_part.Nodes:
            node.X = ox+math.cos(angle)*(node.X - ox)-math.sin(angle)*(node.Y - oy)
            node.Y = oy+math.sin(angle)*(node.X - ox)+math.cos(angle)*(node.Y - oy)
          # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(3)
        ZeroVector[0] = 0.0
        ZeroVector[1] = 0.0
        ZeroVector[2] = 0.0
        for node in self.fluid_model_part.Nodes:
            node.SetValue(MeshingApplication.METRIC_TENSOR_2D, ZeroVector)

        self.MetricParameters = KratosMultiphysics.Parameters("""
        {
            "minimal_size"                         : 2e-3, 
            "enforce_current"                      : true, 
            "anisotropy_remeshing"                 : true, 
            "anisotropy_parameters": 
            {   
                "reference_variable_name"          : "DISTANCE",
                "hmin_over_hmax_anisotropic_ratio"      : 0.1, 
                "boundary_layer_max_distance"           : 0.5, 
                "interpolation"                         : "Linear"
            }
        }
                       
        """)
        
        self.stl_filename = settings["stl_filename"].GetString()

    def Execute(self):
        self.FindWake()
        # self.RefineMesh()
        # self.FindWake()
    def DefineWakeFromLevelSet(self):
        KratosMultiphysics.ModelPartIO("wake").ReadModelPart(self.wake_line_model_part)   
        origin=[0.25,0]  
        angle=math.radians(-self.geometry_parameter) 
        RotateModelPart(origin,angle,self.wake_line_model_part)
        KratosMultiphysics.CalculateDiscontinuousDistanceToSkinProcess2D(self.fluid_model_part, self.wake_line_model_part).Execute()
        for elem in self.fluid_model_part.Elements:
            d=elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
            i=0            
            for node in elem.GetNodes():
                node.SetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE,d[i])
                i += 1
        # KratosMultiphysics.VariableUtils().CopyScalarVar(KratosMultiphysics.DISTANCE,CompressiblePotentialFlow.WAKE_DISTANCE, self.fluid_model_part.Nodes)



        # self.kutta_model_part=self.model.CreateModelPart(self.wake_model_part_name)
        # self.kutta_model_part.AddNode(kutta_node,0)
     
        
        
    def FindWake(self):
        print("Executing wake process")     
        
        if self.wake_model_part_name == "LevelSetWake":
            self.DefineWakeFromLevelSet()
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

            # Neigbour search tool instance
            AvgElemNum = 10
            AvgNodeNum = 10
            nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_model_part,AvgElemNum, AvgNodeNum)
            # Find neighbours
            nodal_neighbour_search.Execute()

            for elem in self.fluid_model_part.Elements:    
                distances=elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                npos = 0
                nneg = 0
                npos_ls = 0
                nneg_ls = 0
                for d in distances:
                    if(abs(d) < self.epsilon):
                        if d > 0.0:
                            d = self.epsilon
                        else:
                            d = -self.epsilon
                    if(d < 0):
                        nneg += 1
                    else:
                        npos += 1
                for node in elem.GetNodes():
                    d=node.GetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE)
                    if(d < 0):
                        nneg_ls += 1
                    else:
                        npos_ls += 1
                if(nneg>0 and npos>0):
                    elem.Set(KratosMultiphysics.MARKER,True)
                    counter = 0
                    for elnode in elem.GetNodes():
                        elnode.SetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE,0,distances[counter])
                        elnode.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,-1)
                        counter+=1
                    elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,distances)
                if(nneg>0 and npos>0) and (nneg_ls>0 and npos_ls>0):
                    elem.Set(KratosMultiphysics.INTERFACE) 
                    elem.Set(KratosMultiphysics.ACTIVE,True)
                    
            from gid_output_process import GiDOutputProcess
            gid_output = GiDOutputProcess(self.fluid_model_part,
                                        "wake_preview",
                                        KratosMultiphysics.Parameters("""
                                            {
                                                "result_file_configuration" : {
                                                    "gidpost_flags": {
                                                        "GiDPostMode": "GiD_PostBinary",
                                                        "WriteDeformedMeshFlag": "WriteUndeformed",
                                                        "WriteConditionsFlag": "WriteConditions",
                                                        "MultiFileFlag": "SingleFile"
                                                    },
                                                    "nodal_results"       : ["DISTANCE","LEVEL_SET_DISTANCE","WAKE_DISTANCE","TEMPERATURE"],
                                                    "elemental_conditional_flags_results": ["MARKER","INTERFACE"]
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
                                        
        else:
            self.kutta_model_part = self.model.GetModelPart(self.wake_model_part_name)         
            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

            # Neigbour search tool instance
            AvgElemNum = 10
            AvgNodeNum = 10
            nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_model_part,AvgElemNum, AvgNodeNum)
            # Find neighbours
            nodal_neighbour_search.Execute()

            #mark as STRUCTURE and deactivate the elements that touch the kutta node
            for node in self.kutta_model_part.Nodes:
                print("Wake Node:",node.Id,": (", node.X,",",node.Y,")")            
                node.Set(KratosMultiphysics.STRUCTURE) 
            #     master_node=node
            # count=1
            # for elem in self.fluid_model_part.Elements:
            #     for node in elem.GetNodes():
            #         if (node.Is(KratosMultiphysics.STRUCTURE)):
            #             elem.Set(KratosMultiphysics.SELECTED)
            #             # elem.Set(KratosMultiphysics.MARKER,True)
            # for elem in self.fluid_model_part.Elements:
            #     if elem.Is(KratosMultiphysics.SELECTED):
            #         print("kutta element:", elem.Id)
            #         for node in elem.GetNodes():
            #             if (node.IsNot(KratosMultiphysics.STRUCTURE) and node.IsNot(KratosMultiphysics.SLAVE)):
            #                 node.Set(KratosMultiphysics.SLAVE)
            #                 print('Master:',master_node.Id,', Slave:',node.Id)
            #                 self.fluid_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", count, master_node, CompressiblePotentialFlow.POSITIVE_POTENTIAL, node, CompressiblePotentialFlow.POSITIVE_POTENTIAL, 1.0, 0)
            #                 # self.fluid_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", count+1, master_node, CompressiblePotentialFlow.NEGATIVE_POTENTIAL, node, CompressiblePotentialFlow.NEGATIVE_POTENTIAL, 1.0, 0)
            #                 count += 1
            
            
            #compute the distances of the elements of the wake, and decide which ones are wak    
            if(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2): #2D case
                
                xn = KratosMultiphysics.Vector(3)
                
                self.n = KratosMultiphysics.Vector(3)
                self.n[0] = -self.direction[1]
                self.n[1] = self.direction[0]
                self.n[2] = 0.0
                print("normal =",self.n)
                
                for node in self.kutta_model_part.Nodes:
                    x0 = node.X
                    y0 = node.Y
                    for elem in self.fluid_model_part.Elements:
                    
                        
        
                        #check in the potentially active portion
                        potentially_active_portion = False
                        for elnode in elem.GetNodes():
                            xn[0] = elnode.X - x0
                            xn[1] = elnode.Y - y0
                            xn[2] = 0.0
                            dx = xn[0]*self.direction[0] + xn[1]*self.direction[1]
                            if(dx > 0): 
                                potentially_active_portion = True
                                break
                            if(elnode.Is(KratosMultiphysics.STRUCTURE)): ##all nodes that touch the kutta nodes are potentiallyactive
                                potentially_active_portion = True
                                break
                            
                            
                        if(potentially_active_portion):                   
                            distances = KratosMultiphysics.Vector(len(elem.GetNodes()))
                            
                            
                            counter = 0
                            for elnode in elem.GetNodes():
                                xn[0] = elnode.X - x0
                                xn[1] = elnode.Y - y0
                                xn[2] = 0.0
                                d =  xn[0]*self.n[0] + xn[1]*self.n[1]
                                if(abs(d) < self.epsilon):
                                    if d >= 0.0:
                                        d = self.epsilon
                                    else:
                                        d = -self.epsilon
                                distances[counter] = d
                                counter += 1

                            npos = 0
                            nneg = 0
                            for d in distances:
                                if(d < 0):
                                    nneg += 1
                                else:
                                    npos += 1
                                    
                            if(nneg>0 and npos>0):
                                elem.Set(KratosMultiphysics.MARKER,True)
                                counter = 0
                                for elnode in elem.GetNodes():
                                    elnode.SetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE,0,distances[counter])
                                    elnode.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,-1)
                                    counter+=1
                                elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,distances)
                                
                                #for elnode in elem.GetNodes():
                                    #if elnode.Is(KratosMultiphysics.STRUCTURE):
                                        #elem.Set(KratosMultiphysics.ACTIVE,False)
                                        #elem.Set(KratosMultiphysics.MARKER,False)
                                        
       
                            
            else: #3D case
                print("3D wake")
                import numpy
                from stl import mesh #this requires numpy-stl

                #initialize the distances to zero on all elements
                zero = KratosMultiphysics.Vector(4)
                zero[0] = 0
                zero[1] = 0
                zero[2] = 0
                zero[3] = 0
                for elem in self.fluid_model_part.Elements:
                    elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, zero)
                    # elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, zero)


                mesh = mesh.Mesh.from_multi_file(self.stl_filename)
                wake_mp = KratosMultiphysics.ModelPart("wake_stl")
                wake_mp.AddNodalSolutionStepVariable(CompressiblePotentialFlow.WAKE_DISTANCE)
                wake_mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
                prop = wake_mp.Properties[0]

                node_id = 1
                elem_id = 1

                for m in mesh:
                    print(m)

                    for vertex in m.points:
                        n1 = wake_mp.CreateNewNode(node_id, float(vertex[0]), float(vertex[1]), float(vertex[2]))
                        node_id+=1
                        n2 = wake_mp.CreateNewNode(node_id, float(vertex[3]), float(vertex[4]), float(vertex[5]))
                        node_id+=1
                        n3 = wake_mp.CreateNewNode(node_id, float(vertex[6]), float(vertex[7]), float(vertex[8]))
                        node_id+=1
                        
                        el = wake_mp.CreateNewElement("Element3D3N",elem_id,  [n1.Id, n2.Id, n3.Id], prop)
                        elem_id += 1
                        
                        
                #CHAPUZA! - to be removed
                #for node in wake_mp.Nodes:
                    #node.X = node.X - 1.0
                    #node.Y = node.Y -0.001
                    
                representation_mp = KratosMultiphysics.ModelPart("representation_mp")
                
                distance_calculator = KratosMultiphysics.CalculateSignedDistanceTo3DSkinProcess(wake_mp, self.fluid_model_part)
                distance_calculator.Execute()

                # KratosMultiphysics.VariableUtils().CopyScalarVar(CompressiblePotentialFlow.WAKE_DISTANCE,CompressiblePotentialFlow.WAKE_DISTANCE, self.fluid_model_part.Nodes)
                # KratosMultiphysics.VariableUtils().CopyVectorVar(KratosMultiphysics.ELEMENTAL_DISTANCES,KratosMultiphysics.ELEMENTAL_DISTANCES, self.fluid_model_part.Nodes)
                # for elem in self.fluid_model_part.Elements:
                #     d=elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                #     # print(d)
                #     elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,d)

                for elem in self.fluid_model_part.Elements:
                    if(elem.Is(KratosMultiphysics.TO_SPLIT)):
                        elem.Set(KratosMultiphysics.MARKER,True)
                
                #the following MORE OR LESS WORKS
                #for elem in self.fluid_model_part.Elements:
                    #kutta_elem = False
                    #for node in elem.GetNodes():
                        #if(node.Is(KratosMultiphysics.STRUCTURE)):
                            #kutta_elem = True

                    #if(kutta_elem == True and elem.IsNot(KratosMultiphysics.MARKER)):
                        #elem.Set(KratosMultiphysics.ACTIVE,False)
                        
                #for elem in self.fluid_model_part.Elements:
                    #kutta_elem = False
                    #for node in elem.GetNodes():
                        #if(node.Is(KratosMultiphysics.STRUCTURE)):
                            #kutta_elem = True

                    #if(kutta_elem == True and elem.IsNot(KratosMultiphysics.MARKER)):
                        #d = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                        #i = 0
                        #for node in elem.GetNodes():
                            #if(node.Is(KratosMultiphysics.STRUCTURE)):
                                #d[i] = -1e-4
                            #else:
                                #d[i] = 1.0
                            #i+=1
                        #elem.SetValue(KratosMultiphysics.WAKE_ELEMENTAL_DISTANCES, d)
                        #elem.Set(KratosMultiphysics.MARKER,True)
                        
                CompressiblePotentialFlow.KuttaConditionProcess(self.fluid_model_part).Execute()
                
                for elem in self.fluid_model_part.Elements:
                    d = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                    for i in range(len(d)):
                        if(abs(d[i]) < self.epsilon ):
                            if(d[i] < 0):
                                d[i] = -self.epsilon
                            else:
                                d[i] = self.epsilon
                    elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,d)
                    i=0
                    for node in elem.GetNodes():
                        i=i+1
                        node.SetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE,d[i])

                        

                    #if(len(d) == 4):
                        #npos = 0
                        #nneg = 0
                        #i = 0
                        #for node in elem.GetNodes():
                            ##d[i] = node.Z - 2502.5
                            #if(d[i] >= 0):
                                #npos += 1
                            #else: 
                                #nneg += 1
                                
                            #i += 1
                            
                        ##elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,d)
                        #if(npos > 0 and nneg>0):
                                #elem.Set(KratosMultiphysics.TO_SPLIT,True)
            
            
                ##chapuza
                #for elem in self.fluid_model_part.Elements:
                    #d = KratosMultiphysics.Vector(4)

                    #if(len(d) == 4):
                        #npos = 0
                        #nneg = 0
                        #i = 0
                        #for node in elem.GetNodes():
                            ##d[i] = node.Z - 2502.5
                            #if(d[i] >= 0):
                                #npos += 1
                            #else: 
                                #nneg += 1
                                
                            #i += 1
                            
                        ##elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,d)
                        #if(npos > 0 and nneg>0):
                                #elem.Set(KratosMultiphysics.TO_SPLIT,True)        
            
            
                #from here it is to output the wake
                distance_calculator.GenerateSkinModelPart(representation_mp);
            
                max_node_id = 1
                for node in self.fluid_model_part.Nodes:
                    if node.Id > max_node_id:
                        max_node_id = node.Id
                        
                max_el_id = 1
                for elem in self.fluid_model_part.Elements:
                    if elem.Id > max_el_id:
                        max_el_id = elem.Id    
                        
                i = max_node_id+1
                for node in representation_mp.Nodes:
                    node.Id = i
                    i += 1
                    
                i = max_el_id+1
                for elem in representation_mp.Elements:
                    elem.Id = i
                    i+=1
            
                from gid_output_process import GiDOutputProcess
                output_file = "representation_of_wake"
                gid_output =  GiDOutputProcess(self.fluid_model_part,
                                        output_file,
                                        KratosMultiphysics.Parameters("""
                                            {
                                                "result_file_configuration": {
                                                    "gidpost_flags": {
                                                        "GiDPostMode": "GiD_PostAscii",
                                                        "WriteDeformedMeshFlag": "WriteUndeformed",
                                                        "WriteConditionsFlag": "WriteConditions",
                                                        "MultiFileFlag": "SingleFile"
                                                    },
                                                    "file_label": "time",
                                                    "output_control_type": "step",
                                                    "output_frequency": 1.0,
                                                    "body_output": true,
                                                    "node_output": false,
                                                    "skin_output": false,
                                                    "plane_output": [],
                                                    "nodal_results": ["DISTANCE"],
                                                    "nodal_nonhistorical_results": [],
                                                    "nodal_flags_results": [],
                                                    "gauss_point_results": [],
                                                    "additional_list_files": []
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
                
    def ExecuteInitialize(self):
        self.Execute()

    def RefineMesh(self):
        
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess2D(self.fluid_model_part, CompressiblePotentialFlow.WAKE_DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.fluid_model_part)
        find_nodal_h.Execute()

        metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess2D(self.fluid_model_part,  KratosMultiphysics.DISTANCE_GRADIENT, self.MetricParameters)
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
        mmg_process = MeshingApplication.MmgProcess2D(self.fluid_model_part, mmg_parameters)

        # We remesh
        mmg_process.Execute()
