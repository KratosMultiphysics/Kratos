#import KratosMultiphysics

#import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
#import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
from KratosMultiphysics import KratosUnittest

from KratosMultiphysics import ConstitutiveModelsApplication
from KratosMultiphysics import UmatApplication



class TestJ2Plasticity(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        try:
            import math as math
        except ImportError as e:
            self.skipTest("Missing python libraries ( importlib and math)")

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)

    def test_Yield(self):

        self._create_material_model_and_law()

        self.properties.SetValue(KratosMultiphysics.OVER_CONSOLIDATION_RATIO, 4.0)
        self.parameters.SetMaterialProperties( self.properties )

        NumberIncrements = 10

        self.strain_vector = self.parameters.GetStrainVector()
        delta_strain = 0.0*self.strain_vector
        delta_strain[0] = 0.01
        delta_strain[1] = -2*0.01
        delta_strain[2] = 0.01

        self._compute_strain_driven_problem(delta_strain, NumberIncrements)
        Pressure, DeviatoricQ = self._calculate_invariants()

        UndrainedShearStrength    = self.properties.GetValue(KratosMultiphysics.YIELD_STRESS)
        self.assertAlmostEqual(DeviatoricQ, UndrainedShearStrength,places=4)

    def _compute_strain_driven_problem(self, delta_strain, nIncr):

        self.parameters.SetDeformationGradientF(self.F)
        self.parameters.SetDeterminantF(self.detF)
        self.strain_vector = self.parameters.GetStrainVector()
        self.strain_vector = 0.0*self.strain_vector
        self.parameters.SetStrainVector(self.strain_vector)
        self.material_law.CalculateMaterialResponseKirchhoff(self.parameters)
        self.material_law.FinalizeMaterialResponseKirchhoff(self.parameters)


        self.strain_vector = self.parameters.GetStrainVector()

        for step in range(1, nIncr+1):

            #Actualize strain
            self.strain_vector = self.strain_vector + delta_strain

            #Compute
            self.parameters.SetStrainVector( self.strain_vector )

            self.material_law.CalculateMaterialResponseCauchy( self.parameters )
            self.material_law.FinalizeMaterialResponseCauchy( self.parameters )

            self.stress = self.parameters.GetStressVector()


    def _create_material_model_and_law(self):

        settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "MaterialDomain",
            "properties_id"   : 1,
            "material_name"   : "steel",
            "constitutive_law": {
                "law_name"   : "KratosMultiphysics.ConstitutiveModelsApplication.SmallStrain3DLaw",
                "model_name" : "KratosMultiphysics.UmatApplication.VonMisesSmallStrainUmatModel"
            },
            "variables": {
                "KratosMultiphysics.YOUNG_MODULUS": 10000.0,
                "KratosMultiphysics.POISSON_RATIO": 0.20,
                "KratosMultiphysics.YIELD_STRESS": 200.0
            },
            "element_type": "Tetrahedra3D4",
            "nodes" : [ [0.0,0.0,0.0], [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ],
            "strain": {
                "deformation_gradient" : [ [1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0] ],
                "jacobian": 1.0
            },
            "echo_level" : 0

        }
        """)
        self.model = KratosMultiphysics.Model()
        self.model_part = self.model.CreateModelPart(settings["model_part_name"].GetString())
        self.echo_level = settings["echo_level"].GetInt()

        #read nodes
        self.number_of_nodes = settings["nodes"].size()
        self.nodes = [] #self.model_part.GetNodes()
        for i in range(0, self.number_of_nodes):
            node = self.model_part.CreateNewNode(i,settings["nodes"][i][0].GetDouble(),settings["nodes"][i][1].GetDouble(),settings["nodes"][i][2].GetDouble())
            self.nodes.append(node)

        self.geometry  = KratosMultiphysics.Geometry()
        self.dimension = 3
        if( settings["element_type"].GetString() == "Tetrahedra3D4"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", settings["element_type"].GetString() )
            else:
                self.geometry = KratosMultiphysics.Tetrahedra3D4(self.nodes[0],self.nodes[1],self.nodes[2],self.nodes[3])


        if( settings["element_type"].GetString() == "Triangle2D3"):
            if( self.number_of_nodes != 4 ):
                print(" number of nodes:",self.number_of_nodes," do not matches geometry :", settings["element_type"].GetString() )
            else:
                self.geometry  = KratosMultiphysics.Triangle2D3(self.nodes[0],self.nodes[1],self.nodes[2])
                self.dimension = 2

        #material properties
        self.properties = self.model_part.Properties[settings["properties_id"].GetInt()]

        #read variables
        self.variables = settings["variables"]
        for key, value in self.variables.items():
            variable = self._GetItemFromModule(key)
            self.properties.SetValue(variable, value.GetDouble())


        #create constitutive law
        self.material_model = self._GetItemFromModule( settings["constitutive_law"]["model_name"].GetString())()
        self.material_law   = self._GetItemFromModule( settings["constitutive_law"]["law_name"].GetString())(self.material_model)

        #check dimension
        if(self.material_law.WorkingSpaceDimension() != self.dimension):
            raise Exception( "mismatch between the ConstitutiveLaw dimension and the dimension of the space")

        #set strain
        self.F = KratosMultiphysics.Matrix(3,3)
        self.strain_measure = settings["strain"]["deformation_gradient"]

        for i in range(0,3):
            for j in range(0,3):
                self.F[i,j] = self.strain_measure[i][j].GetDouble()

        self.detF = settings["strain"]["jacobian"].GetDouble()

        #element parameters
        self.N     = KratosMultiphysics.Vector(self.number_of_nodes)
        self.DN_DX = KratosMultiphysics.Matrix(self.number_of_nodes, self.dimension)

        #set calculation flags
        self.options = KratosMultiphysics.Flags()
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)
        self.options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)

        #set calculation variables
        self.stress_vector       = KratosMultiphysics.Vector(self.material_law.GetStrainSize())
        self.strain_vector       = KratosMultiphysics.Vector(self.material_law.GetStrainSize())
        self.constitutive_matrix = KratosMultiphysics.Matrix(self.material_law.GetStrainSize(),self.material_law.GetStrainSize())


        self.parameters = KratosMultiphysics.ConstitutiveLawParameters()

        self.parameters.SetOptions( self.options )
        self.parameters.SetDeformationGradientF( self.F )
        self.parameters.SetDeterminantF( self.detF )
        self.parameters.SetStrainVector( self.strain_vector )
        self.parameters.SetStressVector( self.stress_vector )
        self.parameters.SetConstitutiveMatrix( self.constitutive_matrix )
        self.parameters.SetShapeFunctionsValues( self.N )
        self.parameters.SetShapeFunctionsDerivatives( self.DN_DX )
        self.parameters.SetProcessInfo( self.model_part.ProcessInfo )
        self.parameters.SetMaterialProperties( self.properties )
        self.parameters.SetElementGeometry( self.geometry )


    def _GetItemFromModule(self,my_string):

        import importlib

        splitted = my_string.split(".")
        if(len(splitted) == 0):
            raise Exception("something wrong. Trying to split the string "+my_string)
        if(len(splitted) == 1):
            return eval(my_string)
        else:
            module_name = ""
            for i in range(len(splitted)-1):
                module_name += splitted[i]
                if i != len(splitted)-2:
                    module_name += "."

            module = importlib.import_module(module_name)
            return getattr(module,splitted[-1])

    def _calculate_invariants(self):

        import math

        #Compute invariants
        Pressure = -(self.stress[0] + self.stress[1] + self.stress[2])/3.0 #Geotech sign convention
        dev = self.stress
        for i in range(0,3):
             dev[i] = dev[i] + Pressure

        J2 = 0
        for i in range(0,3):
            J2 = J2 + dev[i] *dev[i]
        for i in range(3,6):
            J2 = J2 + 2.0*dev[i] *dev[i]
        J2 = math.sqrt(J2*0.5)
        J2 = math.sqrt(3.0) * J2
        DeviatoricQ = J2
        return Pressure, DeviatoricQ

if __name__ == '__main__':
    KratosUnittest.main()
