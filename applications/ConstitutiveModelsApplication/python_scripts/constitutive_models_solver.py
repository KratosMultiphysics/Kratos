from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import os
#import kratos core and applications
import KratosMultiphysics

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(model_part, custom_settings):
    return MaterialsSolver(model_part, custom_settings)

#Base class to develop other solvers
class MaterialsSolver(object):

    def __init__(self, model_part, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "properties_id" : 1,
            "time_settings":{
                "time_step": 1.0,
                "start_time": 0.0,
                "end_time": 1.0
            },
            "integration_settings":{
                "integration_point": []
            }
        }
        """)


        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #validate and assign other values
        self.settings["time_settings"].ValidateAndAssignDefaults(default_settings["time_settings"])
        self.settings["integration_settings"].ValidateAndAssignDefaults(default_settings["integration_settings"])

        self.integration_settings = self.settings["integration_settings"]

        # Set model part
        self.model_part = model_part

        # Set material law
        self.properties   = self.model_part.Properties[self.settings["properties_id"].GetInt()]
        if( self.properties.Has(KratosMultiphysics.CONSTITUTIVE_LAW) ):
            self.material_law = self.properties.GetValue(KratosMultiphysics.CONSTITUTIVE_LAW)
        else:
            raise Exception("ConstitutiveLaw does not exist in properties:",properties_id," id")

        # Process information
        self.process_info = self.model_part.ProcessInfo

        # Set time parameters
        time_settings = self.settings["time_settings"]
        self.process_info.SetValue(KratosMultiphysics.DELTA_TIME, time_settings["time_step"].GetDouble())
        self.process_info.SetValue(KratosMultiphysics.TIME, time_settings["start_time"].GetDouble())

        # Echo level
        self.echo_level = 0

        print("  [Time Step:", self.process_info[KratosMultiphysics.DELTA_TIME]," End_time:", time_settings["end_time"].GetDouble(),"]")


    def SetLawParameters(self, parameters):
        self.parameters = parameters
        
    def GetEndTime(self):
        return (self.settings["time_settings"]["end_time"].GetDouble())

    def ExecuteInitialize(self):       
        print("::[Material_Solver]:: Solver Ready")


    #### Solve loop methods ####

    def InitializeSolutionStep(self):

        #set calculation options
        self._set_calculation_options()
        
        #check material parameters
        self._check_material_parameters()

        pass

    def Execute(self):

        #calculate material response
        self._calculate_material_response()

    def FinalizeSolutionStep(self):
        pass

    def ExecuteFinalize(self):
        
        #set strain parameters
        self._set_calculation_options()

        self.echo_level = 1
        self._calculate_material_response()
        

    #### Solver internal methods ####

    def _set_calculation_options(self):

        #set calculation options to parameters
        self.options = KratosMultiphysics.Flags()

        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        
        #self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRAIN_ENERGY, False)
        #self.options.Set(KratosMultiphysics.ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY, False)
        #self.options.Set(KratosMultiphysics.ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY, False)
        #self.options.Set(KratosMultiphysics.ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE, False)
        #self.options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)

        self.parameters.SetOptions( self.options )

    #
    def _check_material_parameters(self):
        #check parameters
        self.parameters.CheckAllParameters()
        self.parameters.CheckMechanicalVariables()
        self.parameters.CheckShapeFunctions()
                
    #
    def _calculate_material_response(self):

        geometry = self.parameters.GetElementGeometry()
        shape_N  = self.parameters.GetShapeFunctionsValues()

        self.material_law.CalculateMaterialResponsePK2( self.parameters )
        if( self.echo_level > 0 ):
            print("PK2 Material Response")
            print( "stress = ", self.parameters.GetStressVector() )
            print( "strain = ", self.parameters.GetStrainVector() )
            print( "C      = ", self.parameters.GetConstitutiveMatrix() )

        #self.material_law.FinalizeMaterialResponsePK2( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, geometry, shape_N, self.process_info )

        self.material_law.CalculateMaterialResponseKirchhoff( self.parameters )
        if( self.echo_level > 0 ):
            print("\n Kirchhoff Material Response")
            print( "stress = ", self.parameters.GetStressVector() )
            print( "strain = ", self.parameters.GetStrainVector() )
            print( "C      = ", self.parameters.GetConstitutiveMatrix() )

        self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, geometry, shape_N, self.process_info )

        self.material_law.CalculateMaterialResponseCauchy( self.parameters )
        if( self.echo_level > 0 ):
            print("\n Cauchy Material Response")
            print( "stress = ", self.parameters.GetStressVector() )
            print( "strain = ", self.parameters.GetStrainVector() )
            print( "C      = ", self.parameters.GetConstitutiveMatrix() )

        self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
        self.material_law.FinalizeSolutionStep( self.properties, geometry, shape_N, self.process_info )
