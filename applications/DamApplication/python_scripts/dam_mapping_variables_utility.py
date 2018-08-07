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

    def __init__(self,domain_size,main_model_part,post_model_part):

        # Construct the utility
        self.domain_size = domain_size
        if domain_size==2:
            self.MappingVariablesUtility = KratosDam.MappingVariables2DUtilities()
            self.tcl_proc = "Dam_Application::PropagateFractures2D"
        else:
            self.MappingVariablesUtility = KratosDam.MappingVariables3DUtilities()
            self.tcl_proc = "Dam_Application::PropagateFractures3D"

        main_model_part = self.GenereateNewModelPart(main_model_part,
                                                     post_model_part)

    def GenereateNewModelPart(self,main_model_part,post_model_part):

        ### Mapping between old and new model parts ------------------------------------------------------------------

        self.MappingVariablesUtility.MappingModelParts(post_model_part,main_model_part)

        # delete auxiliary model_part
        del post_model_part

        return main_model_part
