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
                "upper_surface_model_part_name" : "Body2D_UpperSurface",
                "lower_surface_model_part_name" : "Body2D_LowerSurface",
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
        for element in self.fluid_model_part.Elements:
            element.Set(KratosMultiphysics.STRUCTURE,False)
            element.Set(KratosMultiphysics.MARKER,False)
        for node in self.fluid_model_part.Nodes:
            node.Set(KratosMultiphysics.STRUCTURE,False)
        self.model=Model
        self.wake_model_part_name=settings["model_part_name"].GetString()
        self.upper_surface_model_part_name=settings["upper_surface_model_part_name"].GetString()
        self.lower_surface_model_part_name=settings["lower_surface_model_part_name"].GetString()

        self.wake_line_model_part=Model.CreateModelPart("wake")
        if (self.fluid_model_part.HasSubModelPart("trailing_edge_model_part")):
            self.fluid_model_part.RemoveSubModelPart("trailing_edge_model_part")
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")
        else:
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")
        angle=math.radians(-self.geometry_parameter)
        for node in self.wake_line_model_part.Nodes:
            node.X = ox+math.cos(angle)*(node.X - ox)-math.sin(angle)*(node.Y - oy)
            node.Y = oy+math.sin(angle)*(node.X - ox)+math.cos(angle)*(node.Y - oy)
        
        self.stl_filename = settings["stl_filename"].GetString()

    def Execute(self):
        self.FindWake()

    def DefineWakeFromLevelSet(self):
        KratosMultiphysics.ModelPartIO("wake").ReadModelPart(self.wake_line_model_part)  
        kutta_model_part_ls=self.fluid_model_part.GetSubModelPart('KuttaLS')
        for node in kutta_model_part_ls.Nodes:
            kutta_node=node 
        # origin=[0.0,0.0] 
        # angle=math.radians(-self.geometry_parameter) 
        for node in self.wake_line_model_part.Nodes:
            node.X=kutta_node.X-0.25+node.X+1e-4
            node.Y=kutta_node.Y+node.Y+1e-4
        # RotateModelPart(origin,angle,self.wake_line_model_part)

        # for node in self.wake_line_model_part.Nodes:
        #     print(node.X)
        #     print(node.Y)
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
    
    def MarkTrailingEdgeElements(self, elem):
        # This function marks the elements touching the trailing
        # edge and saves them in the trailing_edge_model_part for
        # further computations
        if (elem.Is(KratosMultiphysics.ACTIVE)):
            for elnode in elem.GetNodes():
                if(elnode.Is(KratosMultiphysics.STRUCTURE)):
                    elem.Set(KratosMultiphysics.STRUCTURE, True)
                    elem.Set(KratosMultiphysics.MARKER,True)
                    # elem.Set(KratosMultiphysics.BOUNDARY,False)
                    self.trailing_edge_model_part.Elements.append(elem)
                    break
        
    def FindWake(self):
        print("Executing wake process")     
        
        # self.DefineWakeFromLevelSet()
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.fluid_model_part,self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])

        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_model_part,AvgElemNum, AvgNodeNum)
        # Find neighbours
        nodal_neighbour_search.Execute()
        kutta_model_part_ls=self.fluid_model_part.GetSubModelPart('KuttaLS')
        for node in kutta_model_part_ls.Nodes:
            kutta_node=node
            print("Kutta Node Id:", node.Id)

        # for elem in self.fluid_model_part.Elements:    
        #     self.MarkTrailingEdgeElements(elem)     
        #     distances=elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
        #     npos = 0
        #     nneg = 0
        #     npos_ls = 0
        #     nneg_ls = 0
        #     for d in distances:
        #         if(abs(d) < self.epsilon):
        #             if d > 0.0:
        #                 d = self.epsilon
        #             else:
        #                 d = -self.epsilon
        #         if(d < 0):
        #             nneg += 1
        #         else:
        #             npos += 1
        #     for node in elem.GetNodes():
        #         d=node.GetSolutionStepValue(CompressiblePotentialFlow.LEVEL_SET_DISTANCE)
        #         if(d < 0):
        #             nneg_ls += 1
        #         else:
        #             npos_ls += 1
        #     if(nneg>0 and npos>0):
        #         elem.Set(KratosMultiphysics.MARKER,True)
        #         counter = 0
        #         for elnode in elem.GetNodes():
        #             elnode.SetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE,0,distances[counter])
        #             elnode.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,-1)
        #             counter+=1
        #         elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,distances)
        #     if(nneg>0 and npos>0) and (nneg_ls>0 and npos_ls>0):
        #         elem.Set(KratosMultiphysics.ACTIVE,True)

                # center_X=elem.GetGeometry().Center().X
                # if center_X>max_X_center:
                    # max_id=elem.Id
                    # max_X_center=center_X
                #only the max element of the trailing edge is considered part of the wake. 

        #find element with the maximum X value

        xn = KratosMultiphysics.Vector(3)
        
        self.n = KratosMultiphysics.Vector(3)
        self.n[0] = -self.direction[1]
        self.n[1] = self.direction[0]
        self.n[2] = 0.0
        print("normal =",self.n)
        
        x0 = kutta_node.X
        y0 = kutta_node.Y
        for elem in self.fluid_model_part.Elements:
            self.MarkTrailingEdgeElements(elem)      
            elem.Set(KratosMultiphysics.MARKER,False)                     
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
                        counter+=1
                    elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,distances)

        max_X_center = -1e10    
        for element in self.trailing_edge_model_part.Elements:
            # print(element.Id)
            element.Set(KratosMultiphysics.INTERFACE,True)
            if (element.Is(KratosMultiphysics.MARKER)):
                center_X=element.GetGeometry().Center().X
                nneg=0
                for node in element.GetNodes():
                    if node.GetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE)<0:
                        nneg += 1
                if center_X>max_X_center and nneg==1:
                    max_id=element.Id
                    max_X_center=center_X
        # print(max_id)
        #only the max element of the trailing edge is considered part of the wake. 
        for element in self.trailing_edge_model_part.Elements:
            if not (element.Id==max_id):
                element.Set(KratosMultiphysics.MARKER,False)
                nneg=0
                for node in element.GetNodes():
                    if node.GetSolutionStepValue(CompressiblePotentialFlow.WAKE_DISTANCE)<0:
                        nneg += 1
                if not nneg>0:
                    element.Set(KratosMultiphysics.STRUCTURE,False)      

                
                
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
                                               
    def ExecuteInitialize(self):
        self.Execute()

