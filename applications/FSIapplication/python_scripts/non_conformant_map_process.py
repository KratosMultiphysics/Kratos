from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library // TODO: Remove this libraies
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
from KratosMultiphysics.FSIApplication import *

CheckForPreviousImport()

def Factory(settings, Model):
    if(type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return NonConformatMapProcess(Model, settings["Parameters"])

class NonConformatMapProcess(Process):

    def __init__(self,model_part,params):


        ##settings string in json format
        default_parameters = Parameters("""
        {
            "origin_model_part_name": "",
            "destination_model_part_name": "",
            "search_radius_factor": 1.0,
            "it_max": 15,
            "tol": 1.0e-3,
            "sign": "Positive",
            "s2f_origin_variables": ["VELOCITY"],
            "s2f_destination_variables": ["VELOCITY"],
            "f2s_origin_variables": ["PRESSURE"],
            "f2s_destination_variables": ["PRESSURE"]
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.params = params
        self.params.ValidateAndAssignDefaults(default_parameters)

        self.o_model_part = model_part[self.params["origin_model_part_name"].GetString()]
        self.d_model_part = model_part[self.params["destination_model_part_name"].GetString()]

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        search_radius_factor = self.params["search_radius_factor"].GetDouble()
        it_max = self.params["it_max"].GetInt()
        tol = self.params["tol"].GetDouble()

        # Transfer pressure to destination mesh
        self.mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.d_model_part, self.o_model_part, search_radius_factor, it_max, tol)

    def ExecuteInitializeSolutionStep(self):

        for i in range(self.params["s2f_origin_variables"].size()):
              origin = self.params["s2f_origin_variables"][i]
              destination = self.params["s2f_destination_variables"][i]
              if (self.params["sign"].GetString() == "Positive"):
                sign_pos = True
              else:
                sign_pos = False
              variable_origin = KratosGlobals.GetVariable( origin.GetString() )
              variable_destination = KratosGlobals.GetVariable( destination.GetString() )

              val = node.GetSolutionStepValue(variable_origin, 0)
              if isinstance(val,float):
                value_scalar = True
              else:
                value_scalar = False

              if value_scalar:
                self.mapper.StructureToFluid_ScalarMap(variable_origin, variable_destination, sign_pos)
              else: # It is a vector
                self.mapper.StructureToFluid_VectorMap(variable_origin, variable_destination, sign_pos)

        for i in range(self.params["f2s_destination_variables"].size()):
              origin = self.params["f2s_origin_variables"][i]
              destination = self.params["f2s_destination_variables"][i]
              variable_origin = KratosGlobals.GetVariable( origin.GetString() )
              variable_destination = KratosGlobals.GetVariable( destination.GetString() )

              val = node.GetSolutionStepValue(variable_origin, 0)
              if isinstance(val,float):
                value_scalar = True
              else:
                value_scalar = False

              if value_scalar:
                self.mapper.FluidToStructure_ScalarMap(variable_origin, variable_destination)
              else: # It is a vector
                self.mapper.FluidToStructure_VectorMap(variable_origin, variable_destination)

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
