from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestConstitutiveLaw(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def _create_geometry(self, model_part, dim):
        # Create new nodes
        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)

        if (dim == 2):
            nnodes = 3
            node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

            # Allocate a geometry
            geom = KratosMultiphysics.Triangle2D3(node1,node2,node3)
        elif (dim == 3):
            nnodes = 4
            node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
            node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)

            # Allocate a geometry
            geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node3,node4)
        else:
            raise Exception("Error: bad dimension value: ", dim)
        return [geom, nnodes]

    def _create_properties_old(self, model_part, young_modulus, poisson_ratio):
        prop_id = 0
        properties = model_part.Properties[prop_id]
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, young_modulus)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, poisson_ratio)
        return properties

    def _set_cl_parameters(self, cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom):
        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = KratosMultiphysics.ConstitutiveLawParameters()
        cl_params.SetOptions(cl_options)
        cl_params.SetDeformationGradientF(F)
        cl_params.SetDeterminantF(detF)
        cl_params.SetStrainVector(strain_vector)
        cl_params.SetStressVector(stress_vector)
        cl_params.SetConstitutiveMatrix(constitutive_matrix)
        cl_params.SetShapeFunctionsValues(N)
        cl_params.SetShapeFunctionsDerivatives(DN_DX)
        cl_params.SetProcessInfo(model_part.ProcessInfo)
        cl_params.SetMaterialProperties(properties)
        cl_params.SetElementGeometry(geom)

        ## Do all sort of checks
        cl_params.CheckAllParameters() # Can not use this until the geometry is correctly exported to python
        cl_params.CheckMechanicalVariables()
        cl_params.CheckShapeFunctions()
        return cl_params

    def _cl_check(self, cl, properties, geom, model_part, dim):
        cl.Check(properties, geom, model_part.ProcessInfo)

        if(cl.WorkingSpaceDimension() != dim):
            raise Exception("Mismatch between the WorkingSpaceDimension of the Constitutive Law and the dimension of the space in which the test is performed")

    def _set_cl_options(self, dict_options):
        cl_options = KratosMultiphysics.Flags()
        if ("USE_ELEMENT_PROVIDED_STRAIN" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, dict_options["USE_ELEMENT_PROVIDED_STRAIN"])
        if ("COMPUTE_STRESS" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, dict_options["COMPUTE_STRESS"])
        if ("COMPUTE_CONSTITUTIVE_TENSOR" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, dict_options["COMPUTE_CONSTITUTIVE_TENSOR"])
        if ("COMPUTE_STRAIN_ENERGY" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRAIN_ENERGY, dict_options["COMPUTE_STRAIN_ENERGY"])
        if ("ISOCHORIC_TENSOR_ONLY" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY, dict_options["ISOCHORIC_TENSOR_ONLY"])
        if ("VOLUMETRIC_TENSOR_ONLY" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY, dict_options["VOLUMETRIC_TENSOR_ONLY"])
        if ("FINALIZE_MATERIAL_RESPONSE" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE, dict_options["FINALIZE_MATERIAL_RESPONSE"])

        # From here below it should be an otput not an input
        if ("FINITE_STRAINS" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.FINITE_STRAINS, dict_options["FINITE_STRAINS"]) 
        if ("INFINITESIMAL_STRAINS" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.INFINITESIMAL_STRAINS, dict_options["INFINITESIMAL_STRAINS"])
        if ("PLANE_STRAIN_LAW" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.PLANE_STRAIN_LAW, dict_options["PLANE_STRAIN_LAW"])
        if ("PLANE_STRESS_LAW" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.PLANE_STRESS_LAW, dict_options["PLANE_STRESS_LAW"])
        if ("AXISYMMETRIC_LAW" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.AXISYMMETRIC_LAW, dict_options["AXISYMMETRIC_LAW"])
        if ("U_P_LAW" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.U_P_LAW, dict_options["U_P_LAW"])
        if ("ISOTROPIC" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ISOTROPIC, dict_options["ISOTROPIC"])
        if ("ANISOTROPIC" in dict_options): 
            cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ANISOTROPIC, dict_options["ANISOTROPIC"])
        return cl_options

    def _print_cl_output(self, cl, cl_params, properties, geom, N, model_part):
        print("The Material Response PK2")
        cl.CalculateMaterialResponsePK2(cl_params)
        print("Stress = ", cl_params.GetStressVector())
        print("Strain = ", cl_params.GetStrainVector())
        print("C      = ", cl_params.GetConstitutiveMatrix())

        cl.FinalizeMaterialResponsePK2(cl_params)
        cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

        print("\n The Material Response Kirchhoff")
        cl.CalculateMaterialResponseKirchhoff(cl_params)
        print("Stress = ", cl_params.GetStressVector())
        print("Strain = ", cl_params.GetStrainVector())
        print("C      = ", cl_params.GetConstitutiveMatrix())

        cl.FinalizeMaterialResponseKirchhoff(cl_params)
        cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

        print("\n The Material Response Cauchy")
        cl.CalculateMaterialResponseCauchy(cl_params)
        print("Stress = ", cl_params.GetStressVector())
        print("Strain = ", cl_params.GetStrainVector())
        print("C      = ", cl_params.GetConstitutiveMatrix())

        cl.FinalizeMaterialResponseCauchy(cl_params)
        cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

    def _generic_constitutive_law_test(self, model_part, deformation_test):
        # Define geometry
        [geom, nnodes] = self._create_geometry(model_part, deformation_test.cl.dim)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, deformation_test.cl.dim)

        # Material properties
        properties = deformation_test.cl.create_properties(model_part)

        # Construct a constitutive law
        cl = deformation_test.cl.create_constitutive_Law()
        self._cl_check(cl, properties, geom, model_part, deformation_test.cl.dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                        'COMPUTE_STRESS': True,
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = deformation_test.get_init_deformation_gradientF()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)

        # Check the results
        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0

        for i in range(100):
            F = deformation_test.get_deformation_gradientF(i)
            detF = deformation_test.get_determinantF(i)
            cl_params.SetDeformationGradientF(F)
            cl_params.SetDeterminantF(detF)

            # Chauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
            deformation_test.set_reference_stress(reference_stress, i)

            stress = cl_params.GetStressVector()

            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2)

    def test_Uniaxial_KirchhoffSaintVenant_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = UniaxialKirchhoffSaintVenant3D(0.05)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_KirchhoffSaintVenant_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = SimpleShearKirchhoffSaintVenant3D(0.02)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_Plus_Strech_KirchhoffSaintVenant_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = ShearPlusStrechKirchhoffSaintVenant3D()

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Uniaxial_HyperElastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = UniaxialHyperElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_HyperElastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = SimpleShearHyperElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_Plus_Strech_HyperElastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = ShearPlusStrechHyperElastic3D()

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Uniaxial_Linear_Elastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = UniaxialLinearElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_Linear_Elastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = SimpleShearLinearElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_Plus_Strech_Linear_Elastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = ShearPlusStrechLinearElastic3D()

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Uniaxial_Linear_Elastic_Plane_Stress_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = UniaxialLinearElasticPlaneStress2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_Linear_Elastic_Plane_Stress_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = SimpleShearLinearElasticPlaneStress2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Uniaxial_Linear_Elastic_Plane_Stress_Uncoupled_Shear_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = UniaxialElasticPlaneStressUncoupledShear2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

    def test_Shear_Linear_Elastic_Plane_Stress_Uncoupled_Shear_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_Test = SimpleShearElasticPlaneStressUncoupledShear2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_Test)

class Deformation():
    def get_init_deformation_gradientF(self):
        self.F = KratosMultiphysics.Matrix(self.cl.dim,self.cl.dim)
        for i in range(self.cl.dim):
            for j in range(self.cl.dim):
                if(i==j):
                    self.F[i,j] = 1.0
                else:
                    self.F[i,j] = 0.0
        return self.F

class UniaxialDeformation(Deformation):
    def __init__(self, deltaDef):
        Deformation.__init__(self)
        self.deltaDef = deltaDef

    def get_deformation_gradientF(self, i):
        self.F[0,0] = 1.0 + self.deltaDef * i
        return self.F

    def get_determinantF(self, i):
        return 1.0 + self.deltaDef * i

class UniaxialKirchhoffSaintVenant3D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = KirchhoffSaintVenant3D()

    def set_reference_stress(self, reference_stress, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        detF = self.get_determinantF(i)

        reference_stress[0] =( (lame_lambda * 0.5 + lame_mu) * (detF ** 2.0 - 1.0)*(detF ** 2.0) ) / detF
        reference_stress[1] = 0.5*lame_lambda*(detF ** 2.0 - 1.0) / detF
        reference_stress[2] = reference_stress[1]

class UniaxialHyperElastic3D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = HyperElastic3D()

    def set_reference_stress(self, reference_stress, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        detF = self.get_determinantF(i)

        reference_stress[0] = (lame_lambda * math.log(detF) + lame_mu * (detF ** 2.0 - 1.0)) / detF
        reference_stress[1] = (lame_lambda * math.log(detF)) / detF
        reference_stress[2] = reference_stress[1]

class UniaxialLinearElastic3D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = LinearElastic3D()

    def set_reference_stress(self, reference_stress, i):
        c0 = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        F00 = self.get_deformation_gradientF(i)[0,0]

        reference_stress[0] = c0 * (1.0 - self.cl.poisson_ratio) * (F00**2.0-1.0)/2.0
        reference_stress[1] = c0 * self.cl.poisson_ratio * (F00**2.0-1.0)/2.0
        reference_stress[2] = reference_stress[1]

class UniaxialLinearElasticPlaneStress2D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = LinearElasticPlaneStress2D()

    def set_reference_stress(self, reference_stress, i):
        c0 = self.cl.young_modulus / (1.0 - self.cl.poisson_ratio**2)
        F00 = self.get_deformation_gradientF(i)[0,0]

        reference_stress[0] = c0 * (F00**2.0-1.0)/2.0
        reference_stress[1] = c0 * self.cl.poisson_ratio * (F00**2.0-1.0)/2.0

class UniaxialElasticPlaneStressUncoupledShear2D(UniaxialLinearElasticPlaneStress2D):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = ElasticPlaneStressUncoupledShear2D()

class SimpleShearDeformation(Deformation):
    def __init__(self, deltaDef):
        Deformation.__init__(self)
        self.deltaDef = deltaDef

    def get_deformation_gradientF(self, i):
        self.F[0,1] = self.deltaDef * i
        return self.F

    def get_determinantF(self, i):
        return 1.0

class SimpleShearKirchhoffSaintVenant3D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = KirchhoffSaintVenant3D()

    def set_reference_stress(self, reference_stress, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        F01 = self.get_deformation_gradientF(i)[0,1]

        reference_stress[0] = (0.5*lame_lambda + 2*lame_mu) * (F01)**2.0 + (0.5*lame_lambda + lame_mu) * (F01)**4.0
        reference_stress[1] = (0.5*lame_lambda + lame_mu) * (F01)**2.0
        reference_stress[2] = 0.5*lame_lambda  * (F01)**2.0
        reference_stress[3] = lame_mu * (F01) + (0.5*lame_lambda + lame_mu) * (F01)**3.0

class SimpleShearHyperElastic3D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = HyperElastic3D()

    def set_reference_stress(self, reference_stress, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))

        reference_stress[0] = lame_mu * (self.deltaDef * i)**2.0
        reference_stress[3] = lame_mu * (self.deltaDef * i)

class SimpleShearLinearElastic3D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = LinearElastic3D()

    def set_reference_stress(self, reference_stress, i):
        c0 = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        F01 = self.get_deformation_gradientF(i)[0,1]

        reference_stress[0] = c0 * self.cl.poisson_ratio * (F01**2.0)/2.0
        reference_stress[1] = c0 * (1.0 - self.cl.poisson_ratio) * (F01**2.0)/2.0
        reference_stress[2] = reference_stress[0]
        reference_stress[3] = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * 2.0) * F01

class SimpleShearLinearElasticPlaneStress2D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = LinearElasticPlaneStress2D()

    def set_reference_stress(self, reference_stress, i):
        c0 = self.cl.young_modulus / (1.0 - self.cl.poisson_ratio**2)
        F01 = self.get_deformation_gradientF(i)[0,1]

        reference_stress[0] = c0 * self.cl.poisson_ratio * (F01**2.0)/2.0
        reference_stress[1] = c0 * (F01**2.0)/2.0
        reference_stress[2] = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * 2.0) * F01

class SimpleShearElasticPlaneStressUncoupledShear2D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = ElasticPlaneStressUncoupledShear2D()

    def set_reference_stress(self, reference_stress, i):
        c0 = self.cl.young_modulus / (1.0 - self.cl.poisson_ratio**2)
        F01 = self.get_deformation_gradientF(i)[0,1]
        absGamma12 = abs(F01)

        reference_stress[0] = c0 * self.cl.poisson_ratio * (F01**2.0)/2.0
        reference_stress[1] = c0 * (F01**2.0)/2.0
        reference_stress[2] = (self.cl.shear_modulus + self.cl.shear_modulus_gamma12 * absGamma12 + self.cl.shear_modulus_gamma12_2 * absGamma12**2 + self.cl.shear_modulus_gamma12_3 * absGamma12**3 + self.cl.shear_modulus_gamma12_4 * absGamma12**4)* F01

class ShearPlusStrechDeformation(Deformation):
    def __init__(self):
        Deformation.__init__(self)
        self.x1beta = 1.0
        self.x2beta = 1.0
        self.x3beta = math.pi/200

    def get_deformation_gradientF(self, i):
        self.F[0,0] =  math.cos(self.x3beta * i)
        self.F[0,1] = -math.sin(self.x3beta * i)
        self.F[1,0] =  math.sin(self.x3beta * i)
        self.F[1,1] =  math.cos(self.x3beta * i)
        self.F[0,2] =  - self.x1beta * math.sin(self.x3beta * i) - self.x2beta * math.cos(self.x3beta * i)
        self.F[1,2] =    self.x1beta * math.cos(self.x3beta * i) - self.x2beta * math.sin(self.x3beta * i)
        return self.F

    def get_determinantF(self, i):
        return 1.0

class ShearPlusStrechKirchhoffSaintVenant3D(ShearPlusStrechDeformation):
    def __init__(self):
        ShearPlusStrechDeformation.__init__(self)
        self.cl = KirchhoffSaintVenant3D()

    def set_reference_stress(self, reference_stress, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        x1beta = self.x1beta
        x2beta = self.x2beta
        x3beta = self.x3beta

        reference_stress[0]= math.cos(x3beta * i)*(x2beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) + (x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + x2beta*lame_mu*math.cos(x3beta * i) + x1beta*lame_mu*math.sin(x3beta * i)) + math.sin(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 + x1beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)))
        reference_stress[1]= math.cos(x3beta * i)*(x1beta*lame_mu*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) + (x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + x1beta*lame_mu*math.cos(x3beta * i) - x2beta*lame_mu*math.sin(x3beta * i)) + math.sin(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 - x2beta*lame_mu*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)))
        reference_stress[2]=(lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)
        reference_stress[3]= math.sin(x3beta * i)*(x2beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) - math.cos(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 + x1beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i))) - (x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + x2beta*lame_mu*math.cos(x3beta * i) + x1beta*lame_mu*math.sin(x3beta * i))
        reference_stress[4]=(lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + x1beta*lame_mu*math.cos(x3beta * i) - x2beta*lame_mu*math.sin(x3beta * i)
        reference_stress[5]=- (lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) - x2beta*lame_mu*math.cos(x3beta * i) - x1beta*lame_mu*math.sin(x3beta * i)

class ShearPlusStrechHyperElastic3D(ShearPlusStrechDeformation):
    def __init__(self):
        ShearPlusStrechDeformation.__init__(self)
        self.cl = HyperElastic3D()

    def set_reference_stress(self, reference_stress, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        x1beta = self.x1beta
        x2beta = self.x2beta
        x3beta = self.x3beta

        reference_stress[0] = (x2beta * math.cos(i * x3beta) + x1beta * math.sin(i * x3beta))**2.0
        reference_stress[1] = (x1beta * math.cos(i * x3beta) - x2beta * math.sin(i * x3beta))**2.0
        reference_stress[3] = (x2beta * math.cos(i * x3beta) + x1beta * math.sin(i * x3beta)) * (- x1beta * math.cos(i * x3beta) + x2beta * math.sin(i * x3beta))
        reference_stress[4] = x1beta * math.cos(i * x3beta) - x2beta * math.sin(i * x3beta)
        reference_stress[5] = - x2beta * math.cos(i * x3beta) - x1beta * math.sin(i * x3beta)
        reference_stress *= lame_mu

class ShearPlusStrechLinearElastic3D(ShearPlusStrechDeformation):
    def __init__(self):
        ShearPlusStrechDeformation.__init__(self)
        self.cl = LinearElastic3D()

    def set_reference_stress(self, reference_stress, i):
        c0 = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        c1 = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        x1beta = self.x1beta
        x2beta = self.x2beta
        x3beta = self.x3beta

        reference_stress[0] = c0 * self.cl.poisson_ratio * (x1beta**2.0 + x2beta**2.0) / 2.0
        reference_stress[1] = reference_stress[0]
        reference_stress[2] = c0 * (1.0 - self.cl.poisson_ratio) * (x1beta**2.0 + x2beta**2.0) / 2.0
        reference_stress[4] = x2beta * c1
        reference_stress[5] = -c1 * x1beta

class LinearElastic():
    def __init__(self):
        self.young_modulus = 200e9
        self.poisson_ratio = 0.3

    def create_properties(self, model_part):
        prop_id = 0
        properties = model_part.Properties[prop_id]
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, self.young_modulus)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, self.poisson_ratio)
        return properties

class KirchhoffSaintVenant3D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 3

    def create_constitutive_Law(self):
        return StructuralMechanicsApplication.KirchhoffSaintVenant3DLaw()

class HyperElastic3D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 3

    def create_constitutive_Law(self):
        return StructuralMechanicsApplication.HyperElastic3DLaw()

class LinearElastic3D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 3

    def create_constitutive_Law(self):
        return StructuralMechanicsApplication.LinearElastic3DLaw()

class LinearElasticPlaneStress2D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 2

    def create_constitutive_Law(self):
        return StructuralMechanicsApplication.LinearElasticPlaneStress2DLaw()

class ElasticPlaneStressUncoupledShear2D(LinearElasticPlaneStress2D):
    def __init__(self):
        LinearElasticPlaneStress2D.__init__(self)
        self.shear_modulus = 0.2e6 #shear_modulus = 75e9
        self.shear_modulus_gamma12 = -1.6e6
        self.shear_modulus_gamma12_2 = 6.4e6
        self.shear_modulus_gamma12_3 = -9.8e6
        self.shear_modulus_gamma12_4 = 6.7e6

    def create_properties(self, model_part):
        properties = LinearElastic.create_properties(self, model_part)
        properties.SetValue(KratosMultiphysics.SHEAR_MODULUS, self.shear_modulus)
        properties.SetValue(KratosMultiphysics.SHEAR_MODULUS_GAMMA12, self.shear_modulus_gamma12)
        properties.SetValue(KratosMultiphysics.SHEAR_MODULUS_GAMMA12_2, self.shear_modulus_gamma12_2)
        properties.SetValue(KratosMultiphysics.SHEAR_MODULUS_GAMMA12_3, self.shear_modulus_gamma12_3)
        properties.SetValue(KratosMultiphysics.SHEAR_MODULUS_GAMMA12_4, self.shear_modulus_gamma12_4)
        return properties

    def create_constitutive_Law(self):
        return StructuralMechanicsApplication.ElasticPlaneStressUncoupledShear2DLaw()

if __name__ == '__main__':
    KratosUnittest.main()
