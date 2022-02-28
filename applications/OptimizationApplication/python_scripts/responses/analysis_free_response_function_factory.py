# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShapeOptimizationApplication.response_functions.response_function_factory as sho_response_factory

def CreateResponseFunction(response_name,response_type,response_settings,model):

    all_supported_responses = []

    shape_opt_responses = ["plane_based_packaging","mesh_based_packaging","surface_normal_shape_change","face_angle","airfoil_chord_length",
                            "airfoil_perimeter","total_volume"]
    all_supported_responses.extend(shape_opt_responses)

    if response_type == "total_volume": 
        settings_tobe_passed =  KM.Parameters()
        settings_tobe_passed.AddString("response_type",response_type)
        evaluate_model_parts_name_list = response_settings["evaluated_model_parts"].GetStringArray()
        if len(evaluate_model_parts_name_list)!=1:
            raise RuntimeError("analysis_free_response_function_factory: Response {} must have one evaluate_model_part ".format(response_name))
        settings_tobe_passed.AddString("model_part_name",evaluate_model_parts_name_list[0])
        return sho_response_factory.CreateResponseFunction(response_name,settings_tobe_passed,model)        
    else:
        raise RuntimeError("analysis_free_response_function_factory: Response name {} is not supported, the supported analysis free responses are {}.".format(response_name,all_supported_responses))