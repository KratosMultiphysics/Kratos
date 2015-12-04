from KratosMultiphysics import *

#################33please do not change the following class
class KratosProcessFactory(object):
    def __init__(self, Model):
        self.Model = Model #model is a place
        
    def ConstructListOfProcesses( self, process_list ):
        constructed_processes = []
        for item in process_list:
            module = __import__(item["implemented_in_module"])
            interface_file = __import__(item["implemented_in_file"])
            p = interface_file.Factory(item, self.Model)
            constructed_processes.append( p )
            #if( "implemented_in_python" in item): #check if implemented in python or in c++
                #if item["implemented_in_python"] == True: #here we treat the case of python implemented processes
                    #module = __import__(item["implemented_in_module"])
                    #interface_file = __import__(item["implemented_in_file"])
                    #p = interface_file.Factory(item, self.Model)
                    #constructed_processes.append( p )
                    
                #else: #here we create c++ processes
                    #module = __import__(item["implemented_in_module"])
                    #interface_file = __import__(item["implemented_in_file"])
                    #p = interface_file.Factory(item, self.Model)
                    #constructed_processes.append( p )
                    
            #else:
                #raise Exception( "process defined by settings ", item, " does not define if implemented in python or in c++")
            
    
        return constructed_processes
    
    
########## here we generate the common kratos processes --- IMPLEMENTED IN C++ ###################
def Factory(settings, Model):
    
    if(settings["process_name"] == "ApplyConstantScalarValueProcess"):
        params = settings["parameters"]
        model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
        mesh_id = params["mesh_id"]
        variable_name = params["variable_name"] 
        value = params["value"] 
        is_fixed = params["is_fixed"]
    
        options = Flags()
        options.Set(ApplyConstantScalarValueProcess.VARIABLE_IS_FIXED,  is_fixed)
        
      
        return ApplyConstantScalarValueProcess(model_part, variable_name, value,mesh_id,options)
    
    elif(settings["process_name"] == "ApplyConstantVectorValueProcess"):
        params = settings["parameters"]
        model_part = Model.get(  params.get( "model_part_name", "not found!!" ) , "model part not found" )
        mesh_id = int(params["mesh_id"])
        variable_name = params["variable_name"]  
        value = params["value"]
        factor = params["factor"]  
        
        is_fixed_x = params["is_fixed_x"]
        is_fixed_y = params["is_fixed_y"]
        is_fixed_z = params["is_fixed_z"]
        
        options = Flags()
        options.Set(ApplyConstantVectorValueProcess.X_COMPONENT_FIXED,  is_fixed_x)
        options.Set(ApplyConstantVectorValueProcess.Y_COMPONENT_FIXED,  is_fixed_y)
        options.Set(ApplyConstantVectorValueProcess.Z_COMPONENT_FIXED,  is_fixed_z)
            
        new_process = ApplyConstantVectorValueProcess(model_part, variable_name, factor, value,mesh_id, options)
        
        return new_process