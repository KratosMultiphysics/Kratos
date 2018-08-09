from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Import system python
import os
# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.ExternalSolversApplication as KratosSolvers
import KratosMultiphysics.DamApplication as KratosDam
# Import shutil to manage file copying
import shutil

class MappingVariablesUtility:

    def __init__(self,domain_size,main_model_part,post_model_part_mechanical,post_model_part_thermal,type_of_results):

        # Construct the utility
        self.domain_size = domain_size
        if domain_size==2:
            self.MappingVariablesUtility = KratosDam.MappingVariables2DUtilities()
        else:
            self.MappingVariablesUtility = KratosDam.MappingVariables3DUtilities()

        main_model_part = self.GenereateNewModelPart(main_model_part,
                                                     post_model_part_mechanical,
                                                     post_model_part_thermal,
                                                     type_of_results)

    def GenereateNewModelPart(self,main_model_part,post_model_part_mechanical,post_model_part_thermal,type_of_results):

        ### Mapping between old and new model parts ------------------------------------------------------------------
        if type_of_results == "Mechanical":
            self.MappingVariablesUtility.MappingMechanicalModelParts(post_model_part_mechanical,main_model_part)
        if type_of_results == "Thermal":
            self.MappingVariablesUtility.MappingThermalModelParts(post_model_part_thermal,main_model_part)
        if type_of_results == "Thermo-Mechanical":
            self.MappingVariablesUtility.MappingMechanicalModelParts(post_model_part_mechanical,main_model_part)
            self.MappingVariablesUtility.MappingThermalModelParts(post_model_part_thermal,main_model_part)

        # delete auxiliary model_part
        del post_model_part_mechanical
        del post_model_part_thermal

        return main_model_part
