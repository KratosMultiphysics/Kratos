from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
#import kratos core and applications
import KratosMultiphysics

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

from math import *

class compiled_space_time_function:
    def __init__(self, compiled_function ):
        self.compiled_function = compiled_function

    def function(self,x,y,z,t):
        return eval(self.compiled_function)

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
            },
	    "strain_settings":{
		"description": "shear strain time dependent field",
		"deformation_gradient" : [ [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ]
	    },
            "stress_measure": "Kirchhoff",
            "print_output": false
        }
        """)


        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        #validate and assign other values
        self.settings["time_settings"].ValidateAndAssignDefaults(default_settings["time_settings"])
        self.settings["integration_settings"].ValidateAndAssignDefaults(default_settings["integration_settings"])
        self.settings["strain_settings"].ValidateAndAssignDefaults(default_settings["strain_settings"])

        self.integration_settings = self.settings["integration_settings"]
        self.strain_settings = self.settings["strain_settings"]

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


        # Set integration point
        point = self.integration_settings["integration_point"]
        position = KratosMultiphysics.Vector(point.size())
        for i in range(0,point.size()):
            position[i] = point[i].GetDouble()

        self.process_info.SetValue(KratosMultiphysics.INTEGRATION_COORDINATES, position)

        # Set basic parameters
        self.parameters = KratosMultiphysics.ConstitutiveLawParameters()
        self._set_basic_parameters()

        # Set stress measures
        self.PK2       = False
        self.Kirchhoff = False
        self.Cauchy    = False

        measure = self.settings["stress_measure"].GetString()
        self._set_stress_measure(measure)

        # Print output
        self.print_output = self.settings["print_output"].GetBool()

        # Echo level
        self.echo_level = 0

        print("  [Time Step:", self.process_info[KratosMultiphysics.DELTA_TIME]," End_time:", time_settings["end_time"].GetDouble(),"]")


    def GetEndTime(self):
        return (self.settings["time_settings"]["end_time"].GetDouble())

    def Initialize(self):

        if( self.print_output ):
            self._file_open()

        print("::[Material_Solver]:: Solver Ready -"+self.settings["strain_settings"]["description"].GetString()+"-")


    #### Solve loop methods ####

    def InitializeSolutionStep(self):

        #set calculation options
        self._set_calculation_options()

        #set strain parameters
        self._set_strain_parameters()

        #check material parameters
        self._check_material_parameters()

    def Solve(self):

        #calculate material response
        self._calculate_material_response()

    def FinalizeSolutionStep(self):

        #write data in file
        if( self.print_output ):
            self._write_data_in_file()

    def Finalize(self):

        if( self.print_output ):
            self._file_close()

        #set strain parameters
        self._set_calculation_options()

        #set strain parameters
        self._set_strain_parameters()

        self.echo_level = 1
        self._calculate_material_response()


    #### Solver internal methods ####

    #
    def _set_stress_measure(self, measure):

        if( measure == "All" ):
            self.PK2       = True
            self.Kirchhoff = True
            self.Cauchy    = True
        elif( measure == "PK2" ):
            self.PK2       = True
        elif( measure == "Kirchhoff" ):
            self.Kirchhoff = True
        elif( measure == "Cauchy" ):
            self.Cauchy    = True
        else:
            raise Exception("Not valid stress measure:",measure)

    #
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

        if( self.PK2 ):
            self.initialize_calculation_variables()
            self.material_law.CalculateMaterialResponsePK2( self.parameters )
            if( self.echo_level > 0 ):
                print("PK2 Material Response")
                print( "stress = ", self.parameters.GetStressVector() )
                print( "strain = ", self.parameters.GetStrainVector() )
                print( "C      = ", self.parameters.GetConstitutiveMatrix() )

            self.material_law.FinalizeMaterialResponsePK2( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, geometry, shape_N, self.process_info )

        if( self.Kirchhoff ):
            self.initialize_calculation_variables()
            self.material_law.CalculateMaterialResponseKirchhoff( self.parameters )
            if( self.echo_level > 0 ):
                print("\n Kirchhoff Material Response")
                print( "stress = ", self.parameters.GetStressVector() )
                print( "strain = ", self.parameters.GetStrainVector() )
                print( "C      = ", self.parameters.GetConstitutiveMatrix() )

            self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, geometry, shape_N, self.process_info )

        if( self.Cauchy ):
            self.initialize_calculation_variables()
            self.material_law.CalculateMaterialResponseCauchy( self.parameters )
            if( self.echo_level > 0 ):
                print("\n Cauchy Material Response")
                print( "stress = ", self.parameters.GetStressVector() )
                print( "strain = ", self.parameters.GetStrainVector() )
                print( "C      = ", self.parameters.GetConstitutiveMatrix() )

            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, geometry, shape_N, self.process_info )

    #
    def _set_basic_parameters(self):

        #build a dummy geometry
        self._build_dummy_geometry()

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
        strain_size = self.material_law.GetStrainSize()
        self.stress_vector       = KratosMultiphysics.Vector(strain_size)
        self.strain_vector       = KratosMultiphysics.Vector(strain_size)
        self.constitutive_matrix = KratosMultiphysics.Matrix(strain_size,strain_size)

        self.initialize_calculation_variables()

    #
    def initialize_calculation_variables(self):

        strain_size = self.material_law.GetStrainSize()

        #set to zero
        for i in range(0,strain_size):
            self.stress_vector[i] = 0.0;
            self.strain_vector[i] = 0.0;
            for j in range(0,strain_size):
                self.constitutive_matrix[i,j] = 0.0;

        self.parameters.SetStrainVector( self.strain_vector )
        self.parameters.SetStressVector( self.stress_vector )
        self.parameters.SetConstitutiveMatrix( self.constitutive_matrix )

    #
    def _build_dummy_geometry(self):

        #tables depent on this nodal variables (TODO: get them from assign_materials_process.py)
        nodal_variables = {"TEMPERATURE":293.15 , "PRESSURE":0.0}
        for variable in nodal_variables:
            self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))

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


        for node in self.nodes:
            for variable, value in nodal_variables.items():
                node.SetSolutionStepValue(KratosMultiphysics.KratosGlobals.GetVariable(variable), value)

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

        strain = self.strain_settings["deformation_gradient"]

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

        if value.IsNumber():
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
    @classmethod
    def _get_matrix_determinant(self, F):
        value = ((F[0,0]*F[1,1]*F[2,2])+(F[0,1]*F[1,2]*F[2,0])+(F[0,2]*F[1,0]*F[2,1])-(F[0,2]*F[1,1]*F[2,0])-(F[0,0]*F[1,2]*F[2,1])-(F[0,1]*F[1,0]*F[2,2]))
        return value

    #
    def _file_open(self):
        self.file = open("stress_strain_evolution.txt", 'w')
        self.file.write("Time Strain_X Stress_X Strain_Y Stress_Y\n")

    #
    def _write_data_in_file(self):

        time = 0.0
        if( self.process_info.Has(KratosMultiphysics.TIME) ):
            time = self.process_info[KratosMultiphysics.TIME]

        stress = self.parameters.GetStressVector()
        #strain = self.parameters.GetStrainVector()
        strain = self._compute_infinitessimal_strain()

        string = str(time)+" "+str(strain[0])+" "+str(stress[0])+" "+str(strain[1])+" "+str(stress[1])+"\n"

        self.file.write(string)

    #
    def _file_close(self):
        self.file.close()

    #
    def _compute_infinitessimal_strain(self):

        J = self.F

        J[0,0] = J[0,0]-1.0
        J[1,1] = J[1,1]-1.0
        J[1,1] = J[1,1]-1.0

        E = J
        E[0,1] = 0.5 * (E[0,1] + J[1,0])
        E[0,2] = 0.5 * (E[0,2] + J[2,0])
        E[1,2] = 0.5 * (E[1,2] + J[2,1])

        E[1,0] = E[0,1]
        E[1,0] = E[0,2]
        E[2,1] = E[1,2]

        strain = KratosMultiphysics.Vector(6)

        strain[0] = E[0,0]
        strain[1] = E[1,1]
        strain[2] = E[2,2]

        strain[3] = E[0,1]
        strain[4] = E[1,2]
        strain[5] = E[2,0]

        return strain
