from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import KratosMultiphysics

def Factory(custom_settings, Model):
    return KratosMultiphysics.AssignScalarFieldToConditionsProcess(Model[custom_settings["Parameters"]["model_part_name"].GetString()]
                                                                   ,custom_settings["Parameters"] )
                                        
                                                                   

    
