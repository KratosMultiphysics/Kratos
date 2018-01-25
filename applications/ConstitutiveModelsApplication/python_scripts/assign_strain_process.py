import KratosMultiphysics
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
import importlib

import sys
from math import *

class compiled_space_time_function:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def function(self,x,y,z,t):
        return eval(self.compiled_function)

def CreateProcess(custom_settings, model_part, properties_id):
    if(type(custom_settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignStrainProcess(custom_settings, model_part, properties_id)

class AssignStrainProcess(KratosMultiphysics.Process):
    def __init__(self, custom_settings, model_part, properties_id):

        KratosMultiphysics.Process.__init__(self)

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "description": "strain description",
	    "deformation_gradient" : [ [1,0,0], [0,1,0], [0,0,1] ]
        }
        """)

        #overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #build model part
        self.model_part = model_part

        #get material law
        self.properties   = self.model_part.Properties[properties_id]
        if( self.properties.Has(KratosMultiphysics.CONSTITUTIVE_LAW) ):
            self.material_law = self.properties.GetValue(KratosMultiphysics.CONSTITUTIVE_LAW)
        else:
            raise Exception("ConstitutiveLaw does not exist in properties:",properties_id," id")

        #build a dummy geometry
        self.dimension = self.material_law.WorkingSpaceDimension()
        self.nodes = []
        if( self.dimension == 3 ):
            self.number_of_nodes = 4
            self.nodes.append(self.model_part.CreateNewNode(0,0.0,0.0,0.0))
            self.nodes.append(self.model_part.CreateNewNode(1,1.0,0.0,0.0))
            self.nodes.append(self.model_part.CreateNewNode(2,0.0,1.0,0.0))
            self.nodes.append(self.model_part.CreateNewNode(3,0.0,0.0,1.0))
            self.geometry = KratosMultiphysics.Tetrahedra3D4(self.nodes[0],self.nodes[1],self.nodes[2],self.nodes[3])
        if( self.dimension == 2 ):
            self.number_of_nodes = 3
            self.nodes.append(self.model_part.CreateNewNode(0,0.0,0.0,0.0))
            self.nodes.append(self.model_part.CreateNewNode(1,1.0,0.0,0.0))
            self.nodes.append(self.model_part.CreateNewNode(2,0.0,1.0,0.0))
            self.geometry = KratosMultiphysics.Triangle2D3(self.nodes[0],self.nodes[1],self.nodes[2])


        self.parameters = KratosMultiphysics.ConstitutiveLawParameters()

        #set process_info to parameters
        self.process_info = self.model_part.ProcessInfo
        self.parameters.SetProcessInfo( self.process_info )

        #set material properties to parameters
        self.parameters.SetMaterialProperties( self.properties )

        #set geometry properties to parameters
        self.N     = KratosMultiphysics.Vector(self.number_of_nodes)
        self.DN_DX = KratosMultiphysics.Matrix(self.number_of_nodes, self.dimension)

        self.parameters.SetShapeFunctionsValues( self.N )
        self.parameters.SetShapeFunctionsDerivatives( self.DN_DX )
        self.parameters.SetElementGeometry( self.geometry )

        #set calculation variables to parameters
        self.stress_vector       = KratosMultiphysics.Vector(self.material_law.GetStrainSize())
        self.strain_vector       = KratosMultiphysics.Vector(self.material_law.GetStrainSize())
        self.constitutive_matrix = KratosMultiphysics.Matrix(self.material_law.GetStrainSize(),self.material_law.GetStrainSize())

        self.parameters.SetStrainVector( self.strain_vector )
        self.parameters.SetStressVector( self.stress_vector )
        self.parameters.SetConstitutiveMatrix( self.constitutive_matrix )
        
        print("::[Strain]:: Ready ")

    #
    def GetLawParameters(self):
        return self.parameters
    #
    def ExecuteInitialize(self):
        pass

    #
    def Execute(self):

        #set strain parameters
        self._set_strain_parameters()
            
    #
    def ExecuteFinalize(self):
        pass
        
    #
    def _set_strain_parameters(self):

        self.F = KratosMultiphysics.Matrix(3,3)
        self.detF = 1.0

        self._set_strain_matrix(self.F)
        self.detF = self._get_matrix_determinant(self.F)

        #print(" DeterminantF ", self.F)
        #print(" DeterminantF ", self.detF)

        self.parameters.SetDeformationGradientF(self.F)
        self.parameters.SetDeterminantF(self.detF)

    #
    def _set_strain_matrix(self, F):

        strain = self.settings["deformation_gradient"]

        F[0,0] = self._get_strain_value(strain[0][0])
        F[0,1] = self._get_strain_value(strain[0][1])
        F[0,2] = self._get_strain_value(strain[0][2])
        F[1,0] = self._get_strain_value(strain[1][0])
        F[1,1] = self._get_strain_value(strain[1][1])
        F[1,2] = self._get_strain_value(strain[1][2])
        F[2,0] = self._get_strain_value(strain[2][0])
        F[2,1] = self._get_strain_value(strain[2][1])
        F[2,2] = self._get_strain_value(strain[2][2])


    def _get_strain_value(self, value):

        self.value_is_numeric = False

        if value.IsNumber():
            value_is_numeric = True
            return value.GetDouble()
        else:
            function_expression = value.GetString()

            if (sys.version_info > (3, 0)):
                compiled_function = compiled_space_time_function(compile(function_expression, '', 'eval', optimize=2))
            else:
                compiled_function = compiled_space_time_function(compile(function_expression, '', 'eval'))

            time = self.process_info[KratosMultiphysics.TIME]

            # evolution parameters passed via process_info

            # movement equations position
            position = KratosMultiphysics.Vector(self.dimension)
            if( self.process_info.Has(KratosMultiphysics.INTEGRATION_COORDINATES) ):
                position = self.process_info[KratosMultiphysics.INTEGRATION_COORDINATES]

            # time
            time = 0.0
            if( self.process_info.Has(KratosMultiphysics.TIME) ):
                time = self.process_info[KratosMultiphysics.TIME]
            
            value = compiled_function.function(position[0],position[1],position[2],time)

            return value

    #
    def _get_matrix_determinant(self, F):
        value = ((F[0,0]*F[1,1]*F[2,2])+(F[0,1]*F[1,2]*F[2,0])+(F[0,2]*F[1,0]*F[2,1])-(F[0,2]*F[1,1]*F[2,0])-(F[0,0]*F[1,2]*F[2,1])-(F[0,1]*F[1,0]*F[2,2]))
        return value
    
