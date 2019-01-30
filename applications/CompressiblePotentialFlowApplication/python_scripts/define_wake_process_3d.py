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
    
        self.stl_filename = settings["stl_filename"].GetString()

    def Execute(self):
        self.FindWake()       
        
    def FindWake(self):
        #compute the distances of the elements of the wake, and decide which ones are wak    
        if(self.fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2): #2D case
            raise("Using the 4D Wake Process for with a DOMAIN_SIZE=2!")         

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

        wake_mp = self.model.CreateModelPart("wake_stl")
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
            
        representation_mp = self.model.CreateModelPart("representation_mp")
        
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

