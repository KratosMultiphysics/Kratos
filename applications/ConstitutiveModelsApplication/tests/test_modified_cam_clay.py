from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestModifiedCamClayModel(KratosUnittest.TestCase):

    def __init__(self, *args, **kwargs):
        try:
            import importlib
            import math as math
        except ImportError as e:
            self.skipTest("Missing python libraries ( importlib and math)")

        super(KratosUnittest.TestCase, self).__init__(*args, **kwargs)


    def test_UndrainedStressPath_OC(self):


        self._create_material_model_and_law()

        self.properties.SetValue(KratosMultiphysics.OVER_CONSOLIDATION_RATIO, 4.0)
        self.parameters.SetMaterialProperties( self.properties )

        NumberIncrements = 100
        IncrementalF = self._set_identity_matrix()
        IncrementalF[0,1] = 1.0/float(NumberIncrements)

        self._compute_strain_driven_problem(IncrementalF, NumberIncrements)
        Pressure, DeviatoricQ = self._calculate_invariants()


        #Compute Analytical solution
        pc0    = self.properties.GetValue(KratosMultiphysics.PRE_CONSOLIDATION_STRESS)
        OCR    = self.properties.GetValue(KratosMultiphysics.OVER_CONSOLIDATION_RATIO)
        kappa  = self.properties.GetValue(KratosMultiphysics.SWELLING_SLOPE)
        landa  = self.properties.GetValue(KratosMultiphysics.NORMAL_COMPRESSION_SLOPE)
        M      = self.properties.GetValue(KratosMultiphysics.CRITICAL_STATE_LINE)
        alphaS = self.properties.GetValue(KratosMultiphysics.ALPHA_SHEAR)
        if ( abs(alphaS) > 1e-12):
            self.skipTest("The constitutive problem no longer has analytical solution")

        p0 = pc0 / OCR
        BigLambda = ( landa - kappa) / landa
        pressureFailure = p0 * (  (OCR / 2.0 ) ** BigLambda)
        UndrainedShearStrenght = 0.5*p0*M * ( (OCR/2.0)**BigLambda)

        self.assertAlmostEqual(Pressure, pressureFailure)
        self.assertAlmostEqual(0.5*DeviatoricQ, UndrainedShearStrenght)


    def test_UndrainedStressPath_NC(self):


        self._create_material_model_and_law()

        self.properties.SetValue(KratosMultiphysics.OVER_CONSOLIDATION_RATIO, 1.0)
        self.parameters.SetMaterialProperties( self.properties )

        NumberIncrements = 100
        IncrementalF = self._set_identity_matrix()
        IncrementalF[0,1] = 1.0/float(NumberIncrements)

        self._compute_strain_driven_problem(IncrementalF, NumberIncrements)
        Pressure, DeviatoricQ = self._calculate_invariants()


        #Compute Analytical solution
        pc0    = self.properties.GetValue(KratosMultiphysics.PRE_CONSOLIDATION_STRESS)
        OCR    = self.properties.GetValue(KratosMultiphysics.OVER_CONSOLIDATION_RATIO)
        kappa  = self.properties.GetValue(KratosMultiphysics.SWELLING_SLOPE)
        landa  = self.properties.GetValue(KratosMultiphysics.NORMAL_COMPRESSION_SLOPE)
        M      = self.properties.GetValue(KratosMultiphysics.CRITICAL_STATE_LINE)
        alphaS = self.properties.GetValue(KratosMultiphysics.ALPHA_SHEAR)
        if ( abs(alphaS) > 1e-12):
            self.skipTest("The constitutive problem no longer has analytical solution")

        p0 = pc0 / OCR
        BigLambda = ( landa - kappa) / landa
        pressureFailure = p0 * (  (OCR / 2.0 ) ** BigLambda)
        UndrainedShearStrenght = 0.5*p0*M * ( (OCR/2.0)**BigLambda)

        self.assertAlmostEqual(Pressure, pressureFailure)
        self.assertAlmostEqual(0.5*DeviatoricQ, UndrainedShearStrenght)

    def test_IsotropicLoading(self):
        import math

        self._create_material_model_and_law()

        # read material parameters
        pc0    = self.properties.GetValue(KratosMultiphysics.PRE_CONSOLIDATION_STRESS)
        OCR    = self.properties.GetValue(KratosMultiphysics.OVER_CONSOLIDATION_RATIO)
        kappa  = self.properties.GetValue(KratosMultiphysics.SWELLING_SLOPE)
        landa  = self.properties.GetValue(KratosMultiphysics.NORMAL_COMPRESSION_SLOPE)
        p0 = pc0 / OCR


        NumberIncrements = 100
        IncrementalF = self._set_identity_matrix()
        IncrementalF = self._multiply_matrix_by_number(IncrementalF, 0.999)

        for step in range(0, NumberIncrements):

            self._compute_strain_driven_problem(IncrementalF, 1)
            Pressure, DeviatoricQ = self._calculate_invariants()

            #Analytical solution
            epsi_v = math.log(self.detF)
            p = p0 * math.exp(- epsi_v / kappa)
            if ( p > pc0):
                lnp = epsi_v - kappa * math.log(p0) - ( landa - kappa) * math.log(pc0)
                lnp = -lnp / landa
                p = math.exp(lnp)

            self.assertAlmostEqual(Pressure, p)
            self.assertAlmostEqual(DeviatoricQ, 0.0)
            for i in range(3,6):
                self.assertAlmostEqual( self.stress[i], 0.0)

    def _compute_strain_driven_problem(self, IncrF, nIncr):

        self.parameters.SetDeformationGradientF(self.F)
        self.parameters.SetDeterminantF(self.detF)
        self.material_law.CalculateMaterialResponseKirchhoff(self.parameters)
        self.material_law.FinalizeMaterialResponseKirchhoff(self.parameters)
        self.material_law.FinalizeSolutionStep(self.properties, self.geometry, self.N, self.model_part.ProcessInfo)
        for step in range(0, nIncr):

            #Actualize strain
            self.F = IncrF * self.F
            self.detF = self._compute_determinant(self.F)

            #Compute
            self.parameters.SetDeformationGradientF( self.F )
            self.parameters.SetDeterminantF( self.detF )
            self.material_law.CalculateMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeMaterialResponseKirchhoff( self.parameters )
            self.material_law.FinalizeSolutionStep( self.properties, self.geometry, self.N, self.model_part.ProcessInfo )

        self.stress = self.parameters.GetStressVector()

    def _compute_determinant(self, A):

        det = 0

        det = det + A[0,0]*A[1,1]*A[2,2]
        det = det + A[1,0]*A[2,1]*A[0,2]
        det = det + A[2,0]*A[0,1]*A[1,2]

        det = det - A[0,2]*A[1,1]*A[2,0]
        det = det - A[1,2]*A[2,1]*A[0,0]
        det = det - A[2,2]*A[0,1]*A[1,0]

        return det


    def _set_identity_matrix(self):
        identity = KratosMultiphysics.Matrix(3,3)
        for i in range(0,3):
            for j in range(0,3):
                if ( i == j ):
                    identity[i,j] = 1.0
                else:
                    identity[i,j] = 0.0

        return identity

    def _multiply_matrix_by_number(self, A, alpha):
        for i in range(0,3):
            for j in range(0,3):
                A[i,j] = alpha*A[i,j]
        return A

    def _create_material_model_and_law(self):

        settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "MaterialDomain",
            "properties_id"   : 1,
            "material_name"   : "steel",
            "constitutive_law": {
                "law_name"   : "KratosMultiphysics.ConstitutiveModelsApplication.LargeStrain3DLaw",
                "model_name" : "KratosMultiphysics.ConstitutiveModelsApplication.CamClayModel"
            },
            "variables": {
                "KratosMultiphysics.PRE_CONSOLIDATION_STRESS": 80.0,
                "KratosMultiphysics.OVER_CONSOLIDATION_RATIO": 4.0,
                "KratosMultiphysics.SWELLING_SLOPE": 0.0078,
                "KratosMultiphysics.NORMAL_COMPRESSION_SLOPE": 0.085,
                "KratosMultiphysics.INITIAL_SHEAR_MODULUS": 1000.0,
                "KratosMultiphysics.ALPHA_SHEAR": 0.0,
                "KratosMultiphysics.CRITICAL_STATE_LINE": 0.90
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
            if( value.IsDouble() ):
                self.properties.SetValue(variable, value.GetDouble())
            elif( value.IsArray() ):
                vector_value = KratosMultiphysics.Vector(value.size())
                for i in range(0, value.size() ):
                    vector_value[i] = value[i].GetDouble()
                self.properties.SetValue(variable, vector_value)


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
