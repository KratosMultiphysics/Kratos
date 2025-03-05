# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# additional imports
# Kratos Core and Apps
import KratosMultiphysics as KM
from KratosMultiphysics import Logger

# Additional imports
import KratosMultiphysics.OptimizationApplication.responses.analysis_based_response_function_factory as analysis_based_response_function_factory
import KratosMultiphysics.OptimizationApplication.responses.analysis_free_response_function_factory as analysis_free_response_function_factory
import time as timer

# ==============================================================================
def CreateController(reponses_settings,model,model_parts_controller,analyses_controller):
    return ResponsesController(reponses_settings,model,model_parts_controller,analyses_controller)

# ==============================================================================
class ResponsesController:
    # --------------------------------------------------------------------------
    def __init__(self,reponses_settings,model,model_parts_controller,analyses_controller):
        
        self.reponses_settings = reponses_settings
        self.analyses_controller = analyses_controller
        self.model_parts_controller = model_parts_controller
        self.model = model

        required_settings = KM.Parameters("""
        {
            "name"                : "REPONSE_NAME",
            "type"                : "RESPONSE_TYPE",
            "settings"                : {
                "evaluated_objects": [],
                "controlled_objects": [],
                "control_types": []                
            }
        }""")

      
        for itr in range(self.reponses_settings.size()):
            for key in required_settings.keys():
                if not self.reponses_settings[itr].Has(key):
                    raise RuntimeError("ResponsesController: Required entry '{}' missing in 'response Nr.{}'!".format(key,itr+1))     

            for key in required_settings["settings"].keys():
                if not self.reponses_settings[itr]["settings"].Has(key):
                    raise RuntimeError("ResponsesController: Required entry '{}' missing in settings of 'response Nr.{}'!".format(key,itr+1)) 


        analysis_based_responses = ["partition_interface_stress"]
        structural_responses = ["strain_energy", "stress"]        
        analysis_based_responses.extend(structural_responses)

        analysis_free_responses = ["mass","linear","plane_symmetry","interface","partition_mass","am_max_overhang_angle", "dummy_response"]
        shape_opt_responses = ["plane_based_packaging","mesh_based_packaging","surface_normal_shape_change","face_angle","airfoil_chord_length",
                               "airfoil_perimeter","total_volume"]                                
        analysis_free_responses.extend(shape_opt_responses)

        self.supported_responses = []
        self.supported_responses.extend(analysis_based_responses) 
        self.supported_responses.extend(analysis_free_responses)


        self.responses = {}
        self.responses_type = {}
        self.responses_analyses = {}
        self.responses_evaluated_objects = {}
        self.responses_controlled_objects = {}
        self.responses_control_types = {}
        for itr in range(self.reponses_settings.size()):
            response_settings = self.reponses_settings[itr]
            response_name = response_settings["name"].GetString()            
            response_type = response_settings["type"].GetString()

            # check for name
            if  response_name in self.responses.keys():  
                raise RuntimeError("ResponsesController: Response name {} already exists.".format(response_name))

            # check for evaluated_objects
            evaluated_objects = response_settings["settings"]["evaluated_objects"].GetStringArray()                  
            if len(evaluated_objects)>0:
                self.model_parts_controller.CheckIfRootModelPartsExist(evaluated_objects,True)
            else:
                raise RuntimeError("ResponsesController: 'evaluated_objects' of response {} can not be empty.".format(response_name))

            self.responses_evaluated_objects[response_name] = evaluated_objects
            
            # check for controlled_objects   
            controlled_objects = response_settings["settings"]["controlled_objects"].GetStringArray()                      
            if len(controlled_objects)>0:
                self.model_parts_controller.CheckIfRootModelPartsExist(controlled_objects,True)
            else:
                raise RuntimeError("ResponsesController: 'controlled_objects' of response {} can not be empty.".format(response_name))                
            self.responses_controlled_objects[response_name] = controlled_objects

            # check if controlled_objects and evaluated_objects have the same root model parts
            first_root_evaluated_model_part = evaluated_objects[0].split(".")[0]

            for model_part in evaluated_objects:
                if model_part.split(".")[0] != first_root_evaluated_model_part:
                    raise RuntimeError("ResponsesController: 'evaluated_objects' of response {} must have the same root model part.".format(response_name))                     

            first_root_controlled_model_part = controlled_objects[0].split(".")[0]
            for model_part in controlled_objects:
                if model_part.split(".")[0] != first_root_controlled_model_part:
                    raise RuntimeError("ResponsesController: 'controlled_objects' of response {} must have the same root model part.".format(response_name))  

            if first_root_controlled_model_part != first_root_evaluated_model_part:
                raise RuntimeError("ResponsesController: 'controlled_objects' and 'evaluated_objects' of response {} must have the same root model parts.".format(response_name))                  

            # check for control_types
            control_types = response_settings["settings"]["control_types"].GetStringArray()
            if len(control_types) != len(controlled_objects):
                raise RuntimeError("ResponsesController: 'control_types' of response {} must be of the same size as controlled_objects .".format(response_name)) 
            self.responses_control_types[response_name] = control_types              

            response = None
            if response_type in analysis_based_responses:
                if not response_settings["settings"].Has("analysis_name"):
                    raise RuntimeError("ResponsesController: Response {} of type {} requires analysis with analysis_name entry in settings".format(response_name,response_type)) 
                response_analysis_name = response_settings["settings"]["analysis_name"].GetString()
                response_analysis = analyses_controller.GetAnalysis(response_analysis_name)
                response = analysis_based_response_function_factory.CreateResponseFunction(response_name,response_type,response_settings["settings"],response_analysis,model)
                self.responses_analyses[response_name] = response_analysis_name
            elif response_type in analysis_free_responses:
                response = analysis_free_response_function_factory.CreateResponseFunction(response_name,response_type,response_settings["settings"],model)
                self.responses_analyses[response_name] = None
            else:
                raise RuntimeError("ResponsesController: Response {} of type {} is not supported, supported responses are {}".format(response_name,response_type,self.supported_responses)) 
            
            self.responses[response_name] = response

            self.responses_type[response_name] = response_type
            
        # here we do dependency checks
        for itr1 in range(self.reponses_settings.size()):
            response_settings1 = self.reponses_settings[itr1]
            response_name1 = response_settings1["name"].GetString()            
            response_type1 = response_settings1["type"].GetString()
            response_analysis1 = None
            if response_settings1["settings"].Has("analysis_name"):
                response_analysis1 = response_settings1["settings"]["analysis_name"].GetString()
            evaluated_objects1 = response_settings1["settings"]["evaluated_objects"].GetStringArray()
            controlled_objects1 = response_settings1["settings"]["controlled_objects"].GetStringArray() 
            control_types1 = response_settings1["settings"]["control_types"].GetStringArray()

            for itr2 in range(self.reponses_settings.size()):
                if itr1 != itr2:
                    response_settings2 = self.reponses_settings[itr2]
                    response_name2 = response_settings2["name"].GetString()            
                    response_type2 = response_settings2["type"].GetString()    
                    response_analysis2 = None
                    if response_settings2["settings"].Has("analysis_name"):
                        response_analysis2 = response_settings2["settings"]["analysis_name"].GetString()                                    
                    evaluated_objects2 = response_settings2["settings"]["evaluated_objects"].GetStringArray()
                    controlled_objects2 = response_settings2["settings"]["controlled_objects"].GetStringArray() 
                    control_types2 = response_settings2["settings"]["control_types"].GetStringArray()
                    if response_name1 == response_name2:
                        raise RuntimeError("ResponsesController: Response name {} is duplicated.".format(response_name1))
                    if response_type1 == response_type2 and response_analysis1 == response_analysis2:
                        overlap_evaluated_objects =  list(set(evaluated_objects1) & set(evaluated_objects2))
                        if len(overlap_evaluated_objects)>0:
                            overlap_controlled_objects =  list(set(controlled_objects1) & set(controlled_objects2))
                            if len(overlap_controlled_objects)>0:
                                for control_obj in overlap_controlled_objects:
                                    index_1 = controlled_objects1.index(control_obj)
                                    index_2 = controlled_objects2.index(control_obj)
                                    type_1 = control_types1[index_1]
                                    type_2 = control_types2[index_2]
                                    if type_1 == type_2:
                                        raise RuntimeError("ResponsesController: found dependencies between Response {} and Response {}".format(response_name1,response_name2))


    # --------------------------------------------------------------------------
    def Initialize(self):
        for response in self.responses.values():
            response.Initialize()

    # --------------------------------------------------------------------------
    def CheckIfResponseExists(self,response_name,raise_error=True):
        if not response_name in self.responses.keys():
            if raise_error:
                raise RuntimeError("ResponsesController:CheckIfResponseExists: Response name {} does not exist.".format(response_name))
            else: 
                return False
        else:
            return True

    # --------------------------------------------------------------------------
    def CheckIfResponsesExist(self,responses_name,raise_error=True):
        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:CheckIfResponsesExist requires list of response names")
        
        if_exist = True
        for response_name in responses_name:
            if not response_name in self.responses.keys():
                if raise_error:
                    raise RuntimeError("ResponsesController:CheckIfResponsesExist: Response {} does not exist!".format(response_name))
                else:
                    if_exist = False
                    break

        return if_exist                    

    # --------------------------------------------------------------------------
    def CalculateResponseValue(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:CalculateResponseValue: Try to calculate response {} which does not exist.".format(response_name))
        else:
            response_value = self.responses[response_name].CalculateValue()
        return response_value

    # --------------------------------------------------------------------------
    def CalculateResponsesValue(self,responses_name):
        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:CalculateResponsesValue requires list of response names")

        responses_value = {}
        for response_name in responses_name:
            if not response_name in self.responses.keys():
                raise RuntimeError("ResponsesController:CalculateResponsesValue: Try to calculate response {} which does not exist.".format(response_name))            
            responses_value[response_name] = self.responses[response_name].CalculateValue()

        return responses_value        

    # --------------------------------------------------------------------------
    def CalculateAllValues(self):
        responses_values = []
        Logger.PrintInfo("ResponsesController:CalculateAllValues: Starting calculation of values of responses ")
        startTime = timer.time()          
        for name,response in self.responses.items():          
            response_value = response.CalculateValue()
            responses_values.append(response_value)
        Logger.PrintInfo("ResponsesController:CalculateAllValues: Time needed for all response values calculation ",round(timer.time() - startTime,2),"s")            
        return responses_values
    # --------------------------------------------------------------------------
    def CalculateResponseGradients(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:CalculateResponseGradient: Try to calculate all gradients of response {} which does not exist ".format(response_name))
        else:
            self.responses[response_name].CalculateGradients()

    # --------------------------------------------------------------------------
    def CalculateResponsesGradients(self,responses_name):

        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:CalculateResponsesGradients requires list of response names")

        for response_name in responses_name:
            if not response_name in self.responses.keys():
                raise RuntimeError("ResponsesController:CalculateResponsesGradients: Try to calculate gradients of response {} which does not exist.".format(response_name))            
            self.responses[response_name].CalculateGradients()

    # --------------------------------------------------------------------------
    def GetResponseGradients(self,response_name):
        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:GetResponseGradients: Try to calculate all gradients of response {} which does not exist ".format(response_name))
        else:
            return self.responses[response_name].GetGradients()

    # --------------------------------------------------------------------------
    def GetResponseAnalysis(self,response_name):

        if not response_name in self.responses.keys():
            raise RuntimeError("ResponsesController:GetResponseAnalysis: Try to get analysis of response {} which does not exist ".format(response_name))
        
        return self.responses_analyses[response_name]
    # --------------------------------------------------------------------------
    def GetResponsesAnalyses(self,responses_name):
        if type(responses_name) is not list:
            raise RuntimeError("ResponsesController:GetResponsesAnalyses requires list of response names")        

        analyses_list = []
        for response_name in responses_name:
            if self.responses_analyses[response_name] != None:
                analyses_list.append(self.responses_analyses[response_name])
        
        return list(set(analyses_list)) # here we remove duplicates     
    # --------------------------------------------------------------------------
    def GetResponsesForControl(self,control_type,controlled_objects):
   

        if type(controlled_objects) is not list:
            raise RuntimeError("ResponsesController:GetResponsesForControl requires a control type and list of controlled objects")  

        response_dict = {}
        for response_name in self.responses_controlled_objects.keys():
            response_control_types = self.responses_control_types[response_name]
            response_controlled_objects = self.responses_controlled_objects[response_name]
            overlap_objects =  list(set(controlled_objects) & set(response_controlled_objects))
            if len(overlap_objects)>0:
                for object in overlap_objects:
                    index = response_controlled_objects.index(object)
                    object_control_type = response_control_types[index]
                    if object_control_type == control_type:
                        if response_name in response_dict.keys():
                            response_dict[response_name].append(object)
                        else:
                            response_dict[response_name] = [object]


        return response_dict


    # --------------------------------------------------------------------------
    def CalculateResponseGradientsForTypesAndObjects(self,response_name,control_types,controlled_objects,raise_error=True):

        if type(controlled_objects) is not list:
            raise RuntimeError("ResponsesController:CalculateResponseGradientsForTypesAndObjects requires list of controlled objects") 

        if type(control_types) is not list:
            raise RuntimeError("ResponsesController:CalculateResponseGradientsForTypesAndObjects requires list of control types")

        if len(control_types) != len(controlled_objects):
            raise RuntimeError("ResponsesController:CalculateResponseGradientsForTypesAndObjects control types and control object lists should have the same size !")

        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:CalculateResponseGradientsForTypesAndObjects: response {} does not exist".format(response_name))
            
            for i_index in range(len(controlled_objects)):
                controlled_object = controlled_objects[i_index]
                control_type = control_types[i_index]
                found = False
                for r_index in range(len(self.responses_controlled_objects[response_name])):
                    response_controlled_object = self.responses_controlled_objects[response_name][r_index]
                    response_controlled_type = self.responses_control_types[response_name][r_index]
                    if response_controlled_object==controlled_object and control_type==response_controlled_type:
                        found = True
                        break
                if not found:
                    raise RuntimeError("ResponsesController:CalculateResponseGradientsForTypesAndObjects: response {} does not have controlled objects {}".format(response_name,controlled_objects))   

        self.responses[response_name].CalculateGradientsForTypesAndObjects(control_types,controlled_objects,not raise_error)
    # --------------------------------------------------------------------------
    def GetResponseType(self,response_name,raise_error=True): 
        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:CalculateResponseGradientsForTypeAndObjects: response {} does not exist".format(response_name))   

        return self.responses_type[response_name]
    # --------------------------------------------------------------------------
    def GetResponseVariableName(self, response_name, raise_error=True): 
        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:GetResponseVariableName: response {} does not exist".format(response_name))           
        return self.responses[response_name].GetVariableName()

    # --------------------------------------------------------------------------
    def GetResponseControlledObjects(self, response_name, raise_error=True): 
        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:GetResponseControlledObjects: response {} does not exist".format(response_name)) 

        return self.responses_controlled_objects[response_name]
    # --------------------------------------------------------------------------
    def GetResponseControlTypes(self, response_name, raise_error=True): 
        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:GetResponseControlTypes: response {} does not exist".format(response_name)) 

        return self.responses_control_types[response_name]
    # --------------------------------------------------------------------------
    def GetResponseGradientsVariablesName(self, response_name, raise_error=True): 
        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:GetResponseGradientsVariablesName: response {} does not exist".format(response_name))           
        return self.responses[response_name].GetGradientsVariablesName() 

    # --------------------------------------------------------------------------
    def GetResponseGradientVariableNameForType(self, response_name, control_type, raise_error=True): 
        if raise_error:
            if response_name not in self.responses.keys():
                raise RuntimeError("ResponsesController:GetGradientsVariablesName: response {} does not exist".format(response_name))           
        return self.responses[response_name].GetGradientVariableNameForType(control_type)      
                



            

               

