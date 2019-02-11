import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CompressiblePotentialFlow

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

        if (self.fluid_model_part.HasSubModelPart("trailing_edge_model_part")):
            self.fluid_model_part.RemoveSubModelPart("trailing_edge_model_part")
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")
        else:
            self.trailing_edge_model_part = self.fluid_model_part.CreateSubModelPart("trailing_edge_model_part")

    def Execute(self):
        self.FindWake()

        # self.kutta_model_part=self.model.CreateModelPart(self.wake_model_part_name)
        # self.kutta_model_part.AddNode(kutta_node,0)
    def FindKuttaNode(self):
        max_x=-1e30
        upper_surface=self.model.GetModelPart(self.upper_surface_model_part_name)
        lower_surface=self.model.GetModelPart(self.lower_surface_model_part_name)
        for node in upper_surface.Nodes:
            if node.X > max_x:
                max_x=node.X 
                kutta_node=node
        for node in lower_surface.Nodes:
            if node.X > max_x:
                max_x=node.X 
                kutta_node=node
        return kutta_node
    
    def MarkTrailingEdgeElements(self, elem):
        # This function marks the elements touching the trailing
        # edge and saves them in the trailing_edge_model_part for
        # further computations
        for elnode in elem.GetNodes():
            if(elnode.Is(KratosMultiphysics.STRUCTURE)):
                elem.Set(KratosMultiphysics.STRUCTURE, True)
                self.trailing_edge_model_part.Elements.append(elem)
                break
        
        
    def FindWake(self): 
        if self.fluid_model_part.HasSubModelPart(self.wake_model_part_name):
            self.kutta_model_part = self.model.GetModelPart(self.wake_model_part_name)     
        else:
            if (self.model.HasModelPart("Kutta")):
                self.kutta_model_part = self.model.GetModelPart('Kutta')   
            else:
                kutta_node=self.FindKuttaNode()
                self.kutta_model_part = self.model.CreateModelPart('Kutta')                
                self.kutta_model_part.AddNode(kutta_node,0)

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
            node.Set(KratosMultiphysics.STRUCTURE,True) 
        
        
        #compute the distances of the elements of the wake, and decide which ones are wak    
        if(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3): 
            raise("Using the 2D Wake Process with a DOMAIN_SIZE=3!")
            
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
                        
                        #for elnode in elem.GetNodes():
                            #if elnode.Is(KratosMultiphysics.STRUCTURE):
                                #elem.Set(KratosMultiphysics.ACTIVE,False)
                                #elem.Set(KratosMultiphysics.MARKER,False)

        #find element with the maximum X value
        max_X_center = -1e10    
        for element in self.trailing_edge_model_part.Elements:
            # print(element.Id)
            element.Set(KratosMultiphysics.INTERFACE,True)
            if (element.Is(KratosMultiphysics.MARKER)):
                center_X=element.GetGeometry().Center().X
                if center_X>max_X_center:
                    max_id=element.Id
                    max_X_center=center_X
        # print(max_id)
        #only the max element of the trailing edge is considered part of the wake. 
        for element in self.trailing_edge_model_part.Elements:
            if not (element.Id==max_id):
                element.Set(KratosMultiphysics.MARKER,False)
                if len(element.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES))==0:
                    element.Set(KratosMultiphysics.STRUCTURE,False)          
    def ExecuteInitialize(self):
        self.Execute()

