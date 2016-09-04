from KratosMultiphysics import *

#################33please do not change the following class
class KratosProcessFactory(object):
    def __init__(self, Model):
        self.Model = Model #model is a place
        
    def ConstructListOfProcesses( self, process_list ):
        constructed_processes = []
        for i in range(0,process_list.size()):
            item = process_list[i]
            if(item.Has("kratos_module")):
                kratos_module = __import__(item["kratos_module"].GetString())
            if(item.Has("python_module")):
                python_module = __import__(item["python_module"].GetString())
                p = python_module.Factory(item, self.Model)
                constructed_processes.append( p )
            elif(item.Has("implemented_in_module")):
                print("************************************************************************")
                print("******** WARNING USING THE OLD INTERFACE *******************************")
                print("****** kratos_module or python_module should be prescribed instead of **")
                print("****** implemented_in_module and implemented_in_file                  **")
                print("****** old interface is deprecated and will be removed                **")
                print("************************************************************************")
                module = __import__(item["implemented_in_module"].GetString())
                interface_file = __import__(item["implemented_in_file"].GetString())
                p = interface_file.Factory(item, self.Model)
                constructed_processes.append( p )

            #if( "implemented_in_python" in item): #check if implemented in python or in c++
                #if item["implemented_in_python"] == True: #here we treat the case of python implemented processes
                    #kratos_module = __import__(item["kratos_module"])
                    #python_module = __import__(item["python_module"])
                    #p = python_module.Factory(item, self.Model)
                    #constructed_processes.append( p )
                    
                #else: #here we create c++ processes
                    #kratos_module = __import__(item["kratos_module"])
                    #python_module = __import__(item["python_module"])
                    #p = python_module.Factory(item, self.Model)
                    #constructed_processes.append( p )
                    
            #else:
                #raise Exception( "process defined by settings ", item, " does not define if implemented in python or in c++")
            
    
        return constructed_processes
    
    
########## here we generate the common kratos processes --- IMPLEMENTED IN C++ ###################
def Factory(settings, Model):
    if(settings["process_name"].GetString() == "ApplyConstantScalarValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return ApplyConstantScalarValueProcess(model_part, settings["Parameters"])
        
    elif(settings["process_name"].GetString() == "ApplyConstantVectorValueProcess"):
        model_part = Model[settings["Parameters"]["model_part_name"].GetString()]
        return ApplyConstantVectorValueProcess(model_part, settings["Parameters"])
        
    raise Exception("process name not found ",)
        #params = settings["parameters"]
        #model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
        #mesh_id = int(params["mesh_id"])
        #variable = globals().get( params["variable_name"] ) 
        #value = params["value"]
        #factor = params["factor"]  
        
        #is_fixed_x = params["is_fixed_x"]
        #is_fixed_y = params["is_fixed_y"]
        #is_fixed_z = params["is_fixed_z"]
        
        #options = Flags()
        #options.Set(ApplyConstantVectorValueProcess.X_COMPONENT_FIXED,  is_fixed_x)
        #options.Set(ApplyConstantVectorValueProcess.Y_COMPONENT_FIXED,  is_fixed_y)
        #options.Set(ApplyConstantVectorValueProcess.Z_COMPONENT_FIXED,  is_fixed_z)
            
        #new_process = ApplyConstantVectorValueProcess(model_part, variable, factor, value,mesh_id, options)
        
        #return new_process
