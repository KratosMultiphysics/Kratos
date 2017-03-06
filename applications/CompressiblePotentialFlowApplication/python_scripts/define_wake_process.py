import KratosMultiphysics
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess(Model, settings["Parameters"])


class DefineWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                   : 0,
                "kutta_part_name"           : "please specify the model part that contains the kutta nodes",
                "fluid_part_name"           : "model_part",
                "direction"                  : [1.0,0.0,0.0],
                "epsilon"    : 1e-9
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings) 
 
        self.direction = KratosMultiphysics.Vector(3)
        self.direction[0] = settings["direction"][0].GetDouble()
        self.direction[1] = settings["direction"][1].GetDouble()
        self.direction[2] = settings["direction"][2].GetDouble()
        dnorm = math.sqrt(self.direction[0]**2 + self.direction[1]**2 + self.direction[2]**2)
        self.direction[0] /= dnorm
        self.direction[1] /= dnorm
        self.direction[2] /= dnorm
        print(self.direction)
        
        self.epsilon = settings["epsilon"].GetDouble()

        self.kutta_model_part = Model[settings["kutta_part_name"].GetString()]
        self.fluid_model_part = Model[settings["fluid_part_name"].GetString()]
        
    def Execute(self):
        #mark as STRUCTURE and deactivate the elements that touch the kutta node
        for node in self.kutta_model_part.Nodes:
            node.Set(KratosMultiphysics.STRUCTURE)
            #node.Fix(KratosMultiphysics.POSITIVE_FACE_PRESSURE)

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
                                d = self.epsilon
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
                                elnode.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,distances[counter])
                                counter+=1
                            #elem.SetValue(ELEMENTAL_DISTANCE,distances)
                            
                            #for elnode in elem.GetNodes():
                                #if elnode.Is(KratosMultiphysics.STRUCTURE):
                                    #elem.Set(KratosMultiphysics.ACTIVE,False)
                                    #elem.Set(KratosMultiphysics.MARKER,False)
                                    

                            
        else: #3D case
            raise Exception("wake detection not yet implemented in 3D")
        
    def ExecuteInitialize(self):
        self.Execute()
            

