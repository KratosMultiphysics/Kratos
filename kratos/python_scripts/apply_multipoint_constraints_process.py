from KratosMultiphysics import *
import python_process
   
    
def Factory(settings, Model):
    params = settings["parameters"]
    main_model_part_name = params
    main_model_part = Model[settings["main_model_part_name"].GetString()]

    return ApplyMultipointConstraintsProcess(main_model_part, params)