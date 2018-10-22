from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Import applications
import KratosMultiphysics.DamApplication as KratosDam

class MappingVariablesUtility:

    def __init__(self,domain_size):

        # Construct the utility
        self.domain_size = domain_size
        if domain_size==2:
            self.MappingVariablesUtility = KratosDam.MappingVariables2DUtilities()
        else:
            self.MappingVariablesUtility = KratosDam.MappingVariables3DUtilities()


    def AddPreviousModelPartMechanical(self,main_model_part,post_model_part_mechanical,add_displacement,add_stress):
        ### Mapping between old and new model parts ------------------------------------------------------------------
        self.MappingVariablesUtility.MappingMechanicalModelParts(post_model_part_mechanical,main_model_part,add_displacement,add_stress)

        # delete auxiliary model_part
        del post_model_part_mechanical

        return main_model_part

    def AddPreviousModelPartThermal(self,main_model_part,post_model_part_thermal,add_temperature,add_reference_temperature):

        ### Mapping between old and new model parts ------------------------------------------------------------------
        self.MappingVariablesUtility.MappingThermalModelParts(post_model_part_thermal,main_model_part,add_temperature,add_reference_temperature)

        # delete auxiliary model_part
        del post_model_part_thermal

        return main_model_part

    def AddPreviousModelPartThermoMechanical(self,main_model_part,post_model_part_mechanical,post_model_part_thermal,add_displacement,add_stress,add_temperature,add_reference_temperature):

        ### Mapping between old and new model parts ------------------------------------------------------------------
        self.MappingVariablesUtility.MappingMechanicalModelParts(post_model_part_mechanical,main_model_part,add_displacement,add_stress)
        self.MappingVariablesUtility.MappingThermalModelParts(post_model_part_thermal,main_model_part,add_temperature,add_reference_temperature)

        # delete auxiliary model_part
        del post_model_part_mechanical
        del post_model_part_thermal

        return main_model_part



