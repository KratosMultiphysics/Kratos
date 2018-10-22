import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication
import math
from math import *

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess(Model, settings["Parameters"])


class DefineWakeProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"                       : 0,
                "model_part_name"               : "please specify the model part that contains the kutta nodes",
                "fluid_part_name"               : "MainModelPart",
                "direction"                     : [1.0,0.0,0.0],
                "stl_filename"                  : "please specify name of stl file",
                "epsilon"    : 1e-9,
                "AOAdeg" : 0
            }
            """)

        settings.ValidateAndAssignDefaults(default_settings)
        
        
        self.AOAdeg                 = settings["AOAdeg"].GetDouble()
        
        #convert angle to radians
        self.AOArad = self.AOAdeg*pi/180       
 
        '''
        self.direction = KratosMultiphysics.Vector(3)
        self.direction[0] = settings["direction"][0].GetDouble()*cos(self.AOArad)
        self.direction[2] = settings["direction"][0].GetDouble()*sin(self.AOArad)
        self.direction[1] = settings["direction"][2].GetDouble()  
        '''
        self.direction = KratosMultiphysics.Vector(3)
        self.direction[0] = settings["direction"][0].GetDouble()*cos(self.AOArad)
        self.direction[1] = settings["direction"][0].GetDouble()*sin(self.AOArad)
        self.direction[2] = settings["direction"][2].GetDouble()
        #'''
        dnorm = math.sqrt(self.direction[0]**2 + self.direction[1]**2 + self.direction[2]**2)
        self.direction[0] /= dnorm
        self.direction[1] /= dnorm
        self.direction[2] /= dnorm
        print('wake direction =', self.direction)
        
        self.epsilon = settings["epsilon"].GetDouble()

        self.kutta_model_part =         Model[settings["model_part_name"].GetString()]
        self.fluid_model_part =         Model[settings["fluid_part_name"].GetString()]
        #self.upper_surface_model_part = Model[settings["upper_surface_model_part_name"].GetString()]
        #self.lower_surface_model_part = Model[settings["lower_surface_model_part_name"].GetString()]
        
        # Neigbour search tool instance
        AvgElemNum = 10
        AvgNodeNum = 10
        nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fluid_model_part,AvgElemNum, AvgNodeNum)
        # Find neighbours
        nodal_neighbour_search.Execute()       
        
        self.stl_filename = settings["stl_filename"].GetString()
        
    def Execute(self):
        #mark as STRUCTURE and deactivate the elements that touch the kutta node
        for node in self.kutta_model_part.Nodes:
            node.Set(KratosMultiphysics.STRUCTURE)
            x1 = node.X
            y1 = node.Y
            z1 = node.Z
            
                   

        #find wake node in the outflow boundary
        pos = 0    
        #print(self.kutta_model_part)            
        

        #compute the distances of the elements of the wake, and decide which ones are wak    
        if(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2): #2D case
            #self.projection = KratosMultiphysics.Matrix(2,2)
            #self.projection[0,0] = 1.0
            #self.projection[0,1] = 0.0
            #self.projection[1,1] = 1.0
            #self.projection[1,0] = 0.0
            #self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.UPPER_PROJECTION,self.projection)
            
            xn = KratosMultiphysics.Vector(3)
            
            self.n = KratosMultiphysics.Vector(3)
            self.n[0] = -self.direction[1]
            self.n[1] = self.direction[0]
            self.n[2] = 0.0
            print("wake normal =",self.n)
            
            self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE_NORMAL,self.n)
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
                            #if(xn[0] < 0 and xn[1] > 0 and xn[1] < 0.0001):#for high angles of attack
                            if(xn[0] < 0 and d > 0 and d < 0.0003):#for high angles of attack (selecting nodes in the lower surface)
                                d = -self.epsilon
                                print(elnode)
                            #    print(elnode.X - x0)
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
                                #test to check whether the value of the distance affects the solution. IT DOES!
                                '''
                                if(distances[counter]>0):
                                    distances[counter] = 1.0
                                else:
                                    distances[counter] = -1.0
                                '''
                                elnode.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,distances[counter])
                                counter+=1
                                #dx = elnode.X - x1
                                #dy = elnode.Y - y1
                                #dz = elnode.Z - z1
                                #tmp = dx*self.direction[0] + dy*self.direction[1] + dz*self.direction[2]
                                #if(tmp > pos):
                                #    pos = tmp                               
                            
                            elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,distances)
                            
                            #for elnode in elem.GetNodes():
                                #if elnode.Is(KratosMultiphysics.STRUCTURE):
                                    #elem.Set(KratosMultiphysics.ACTIVE,False)
                                    #elem.Set(KratosMultiphysics.MARKER,False)
            
            '''
            #find the wake elements at the outflow boundary
            for elem in self.fluid_model_part.Elements:
                if(elem.Is(KratosMultiphysics.MARKER)):
                    counter = 0
                    for elnode in elem.GetNodes():
                        dx = elnode.X - x1
                        dy = elnode.Y - y1
                        dz = elnode.Z - z1
                        tmp = dx*self.direction[0] + dy*self.direction[1] + dz*self.direction[2]
                        
                        if(tmp > pos-1e-9):
                            counter +=1
                    if(counter > 0):
                        elem.Set(KratosMultiphysics.BOUNDARY,True)
                        
            '''
        
        else: #3D case
            import numpy
            print('\nImporting wake ...')
            from stl import mesh #this requires numpy-stl
            print('\nWake has been imported')
            #'''
            self.up_projection = KratosMultiphysics.Matrix(3,3)
            self.up_projection[0,0] = 1.0
            self.up_projection[0,1] = 0.0
            self.up_projection[0,2] = 0.0
            
            self.up_projection[1,0] = 0.0
            self.up_projection[1,1] = 1.0
            self.up_projection[1,2] = 0.0
            
            self.up_projection[2,0] = 0.0
            self.up_projection[2,1] = 0.0
            self.up_projection[2,2] = 1.0
            
            self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.UPPER_PROJECTION,self.up_projection)
            
            self.low_projection = KratosMultiphysics.Matrix(3,3)
            self.low_projection[0,0] = 1.0
            self.low_projection[0,1] = 0.0
            self.low_projection[0,2] = 0.0
            
            self.low_projection[1,0] = 0.0
            self.low_projection[1,1] = 1.0
            self.low_projection[1,2] = 0.0
            
            self.low_projection[2,0] = 0.0
            self.low_projection[2,1] = 0.0
            self.low_projection[2,2] = 1.0
            
            self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.LOWER_PROJECTION,self.low_projection)
            #'''
            
            self.wake_normal = KratosMultiphysics.Vector(3)
            self.wake_normal[0] = -self.direction[2]
            self.wake_normal[1] = 0.0
            self.wake_normal[2] = self.direction[0]
            print('wake normal =', self.wake_normal)        
            
            self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.CompressiblePotentialFlowApplication.WAKE_NORMAL,self.wake_normal)            

            #initialize the distances to zero on all elements
            zero = KratosMultiphysics.Vector(4)
            zero[0] = 0
            zero[1] = 0
            zero[2] = 0
            zero[3] = 0
            for elem in self.fluid_model_part.Elements:
                elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, zero)


            mesh = mesh.Mesh.from_multi_file(self.stl_filename)
            wake_mp = KratosMultiphysics.ModelPart("wake_stl")
            wake_mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
            wake_mp.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
            prop = wake_mp.Properties[0]
            
            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 2)#DOES NOT WORK BECAUSE IT IS INITIALIZE AFTERWARDS
                node.SetSolutionStepValue(KratosMultiphysics.WATER_PRESSURE, 0, 0)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0, 0)

            '''
            for node in self.fluid_model_part.Nodes:
                        if(node.GetSolutionStepValue(KratosMultiphysics.CompressiblePotentialFlowApplication.UPPER_SURFACE)==1):
                            node.SetSolutionStepValue(KratosMultiphysics.WATER_PRESSURE,0,5)
            '''
            
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
            
            print('\nComputing distance ...')
            distance_calculator = KratosMultiphysics.CalculateSignedDistanceTo3DSkinProcess(wake_mp, self.fluid_model_part)
            distance_calculator.Execute()
            print('\nDistance has been computed')
            
            for elem in self.fluid_model_part.Elements:
                if(elem.Is(KratosMultiphysics.TO_SPLIT)):
                    elem.Set(KratosMultiphysics.MARKER,True)
                    #for node in elem.GetNodes():
                    #    node.SetSolutionStepValue(KratosMultiphysics.WATER_PRESSURE,0,5)
            
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
                    #elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, d)
                    #elem.Set(KratosMultiphysics.MARKER,True)
                    
            
            print('\nExecuiting Kutta Condition Process')
            KratosMultiphysics.CompressiblePotentialFlowApplication.KuttaConditionProcess(self.fluid_model_part).Execute()            
            print('\nKutta Condition Process ended')
            
            
            for node in self.fluid_model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 0)#DOES WORK
                
                       
                
            #for node in self.upper_surface_model_part.Nodes:
            #    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,3)
            
            
            
            for elem in self.fluid_model_part.Elements:
                d = elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES)
                '''
                counter = 0
                for node in elem.GetNodes():                    
                    if(node.Is(KratosMultiphysics.STRUCTURE)):
                        d[counter] = -0.0001
                    node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,d[counter])
                    counter +=1
                '''
                '''
                for i in range(len(d)):
                    if(abs(d[i]) < self.epsilon ):
                        d[i] = self.epsilon
                        #if(d[i] < 0):
                        #    d[i] = -self.epsilon
                        #else:
                        #    d[i] = self.epsilon
                '''
                #elem.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES,d)
                if(elem.Is(KratosMultiphysics.MARKER)):
                    counter = 0
                    for node in elem.GetNodes():
                        node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0,d[counter])
                        if(d[counter]>0):
                            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,1)
                        else:
                            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,-1)
                        counter+=1
                
                     

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
            gid_output =  GiDOutputProcess(representation_mp,
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
                                                "nodal_results": [],
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
            
            print('Wake Process Done')


                    
                    #print(elem.Id, elem.GetValue(KratosMultiphysics.ELEMENTAL_DISTANCES))
                
    def ExecuteInitialize(self):
        self.Execute()
