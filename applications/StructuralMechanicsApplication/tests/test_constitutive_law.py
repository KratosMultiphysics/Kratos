from __future__ import print_function, absolute_import, division
import KratosMultiphysics

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestConstitutiveLaw(KratosUnittest.TestCase):
    def setUp(self):
        pass

    @staticmethod
    def _create_geometry(model_part, dim):
        # Create new nodes
        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)

        if (dim == 2):
            nnodes = 3

            # Allocate a geometry
            geom = KratosMultiphysics.Triangle2D3(node1,node2,node3)
        elif (dim == 3):
            nnodes = 4
            node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)

            # Allocate a geometry
            geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node3,node4)
        else:
            raise Exception("Error: bad dimension value: ", dim)
        return [geom, nnodes]

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

        print("\nThe Material Response Kirchhoff")
        cl.CalculateMaterialResponseKirchhoff(cl_params)
        print("Stress = ", cl_params.GetStressVector())
        print("Strain = ", cl_params.GetStrainVector())
        print("C      = ", cl_params.GetConstitutiveMatrix())

        cl.FinalizeMaterialResponseKirchhoff(cl_params)
        cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

        print("\nThe Material Response Cauchy")
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
        cl.InitializeMaterial(properties, geom, N)

        # Check the results
        deformation_test.initialize_reference_stress(cl.GetStrainSize())

        for i in range(deformation_test.nr_timesteps):
            deformation_test.set_deformation(cl_params, i)

            # Chauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
            reference_stress = deformation_test.get_reference_stress(i)

            stress = cl_params.GetStressVector()

            tolerance = 1.0e-4
            for j in range(cl.GetStrainSize()):
                if (abs(stress[j]) > tolerance):
                    self.assertAlmostEqual((reference_stress[j] - stress[j])/stress[j], 0.0, msg=("Error checking solution " + str(stress[j]) + " different from " + str(reference_stress[j]) + " with tolerance of " + str(tolerance)), delta=tolerance)

    def test_Uniaxial_KirchhoffSaintVenant_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = UniaxialKirchhoffSaintVenant3D(0.05)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_KirchhoffSaintVenant_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = SimpleShearKirchhoffSaintVenant3D(0.02)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_Plus_Strech_KirchhoffSaintVenant_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = ShearPlusStrechKirchhoffSaintVenant3D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Uniaxial_HyperElastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = UniaxialHyperElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_HyperElastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = SimpleShearHyperElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_Plus_Strech_HyperElastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = ShearPlusStrechHyperElastic3D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Uniaxial_Linear_Elastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = UniaxialLinearElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_Linear_Elastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = SimpleShearLinearElastic3D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_Plus_Strech_Linear_Elastic_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = ShearPlusStrechLinearElastic3D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Uniaxial_Linear_Elastic_Plane_Stress_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = UniaxialLinearElasticPlaneStress2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_Linear_Elastic_Plane_Stress_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = SimpleShearLinearElasticPlaneStress2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Uniaxial_Linear_Elastic_Plane_Stress_Uncoupled_Shear_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = UniaxialElasticPlaneStressUncoupledShear2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Shear_Linear_Elastic_Plane_Stress_Uncoupled_Shear_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = SimpleShearElasticPlaneStressUncoupledShear2D(0.2)

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_J2_Plasticity_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = DeformationLinearJ2Plasticity3D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_J2_Plasticity_Plane_Strain_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = DeformationLinearJ2PlasticityPlaneStrain2D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Isotropic_Damage_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = DeformationLinearIsotropicDamage3D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Isotropic_Damage_Plane_Strain_2D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = DeformationLinearIsotropicDamagePlaneStrain2D()

        self._generic_constitutive_law_test(model_part, deformation_test)

    def test_Small_Strain_Isotropic_Plasticity_3D(self):
        # Define a model
        model_part = KratosMultiphysics.ModelPart("test")

        deformation_test = DeformationSmallStrainIsotropicPlasticity3D()

        self._generic_constitutive_law_test(model_part, deformation_test)

class Deformation():
    def __init__(self):
        self.nr_timesteps = 100

    def get_init_deformation_gradientF(self):
        self.F = KratosMultiphysics.Matrix(self.cl.dim,self.cl.dim)
        for i in range(self.cl.dim):
            for j in range(self.cl.dim):
                if(i==j):
                    self.F[i,j] = 1.0
                else:
                    self.F[i,j] = 0.0
        return self.F

    def initialize_reference_stress(self, strain_size):
        self.reference_stress = KratosMultiphysics.Vector(strain_size)
        for i in range(strain_size):
            self.reference_stress[i] = 0.0

    def set_deformation(self, cl_params, i):
        F = self.get_deformation_gradientF(i)
        detF = self.get_determinantF(i)
        cl_params.SetDeformationGradientF(F)
        cl_params.SetDeterminantF(detF)

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

    def get_reference_stress(self, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        detF = self.get_determinantF(i)

        self.reference_stress[0] =( (lame_lambda * 0.5 + lame_mu) * (detF ** 2.0 - 1.0)*(detF ** 2.0) ) / detF
        self.reference_stress[1] = 0.5*lame_lambda*(detF ** 2.0 - 1.0) / detF
        self.reference_stress[2] = self.reference_stress[1]
        return self.reference_stress

class UniaxialHyperElastic3D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = HyperElastic3D()

    def get_reference_stress(self, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        detF = self.get_determinantF(i)

        self.reference_stress[0] = (lame_lambda * math.log(detF) + lame_mu * (detF ** 2.0 - 1.0)) / detF
        self.reference_stress[1] = (lame_lambda * math.log(detF)) / detF
        self.reference_stress[2] = self.reference_stress[1]
        return self.reference_stress

class UniaxialLinearElastic3D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = LinearElastic3D()

    def get_reference_stress(self, i):
        c0 = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        F00 = self.get_deformation_gradientF(i)[0,0]

        self.reference_stress[0] = c0 * (1.0 - self.cl.poisson_ratio) * (F00**2.0-1.0)/2.0
        self.reference_stress[1] = c0 * self.cl.poisson_ratio * (F00**2.0-1.0)/2.0
        self.reference_stress[2] = self.reference_stress[1]
        return self.reference_stress

class UniaxialLinearElasticPlaneStress2D(UniaxialDeformation):
    def __init__(self, deltaDef):
        UniaxialDeformation.__init__(self, deltaDef)
        self.cl = LinearElasticPlaneStress2D()

    def get_reference_stress(self, i):
        c0 = self.cl.young_modulus / (1.0 - self.cl.poisson_ratio**2)
        F00 = self.get_deformation_gradientF(i)[0,0]

        self.reference_stress[0] = c0 * (F00**2.0-1.0)/2.0
        self.reference_stress[1] = c0 * self.cl.poisson_ratio * (F00**2.0-1.0)/2.0
        return self.reference_stress

class UniaxialElasticPlaneStressUncoupledShear2D(UniaxialLinearElasticPlaneStress2D):
    def __init__(self, deltaDef):
        UniaxialLinearElasticPlaneStress2D.__init__(self, deltaDef)
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

    def get_reference_stress(self, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        F01 = self.get_deformation_gradientF(i)[0,1]

        self.reference_stress[0] = (0.5*lame_lambda + 2*lame_mu) * (F01)**2.0 + (0.5*lame_lambda + lame_mu) * (F01)**4.0
        self.reference_stress[1] = (0.5*lame_lambda + lame_mu) * (F01)**2.0
        self.reference_stress[2] = 0.5*lame_lambda  * (F01)**2.0
        self.reference_stress[3] = lame_mu * (F01) + (0.5*lame_lambda + lame_mu) * (F01)**3.0
        return self.reference_stress

class SimpleShearHyperElastic3D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = HyperElastic3D()

    def get_reference_stress(self, i):
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))

        self.reference_stress[0] = lame_mu * (self.deltaDef * i)**2.0
        self.reference_stress[3] = lame_mu * (self.deltaDef * i)
        return self.reference_stress

class SimpleShearLinearElastic3D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = LinearElastic3D()

    def get_reference_stress(self, i):
        c0 = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        F01 = self.get_deformation_gradientF(i)[0,1]

        self.reference_stress[0] = c0 * self.cl.poisson_ratio * (F01**2.0)/2.0
        self.reference_stress[1] = c0 * (1.0 - self.cl.poisson_ratio) * (F01**2.0)/2.0
        self.reference_stress[2] = self.reference_stress[0]
        self.reference_stress[3] = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * 2.0) * F01
        return self.reference_stress

class SimpleShearLinearElasticPlaneStress2D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = LinearElasticPlaneStress2D()

    def get_reference_stress(self, i):
        c0 = self.cl.young_modulus / (1.0 - self.cl.poisson_ratio**2)
        F01 = self.get_deformation_gradientF(i)[0,1]

        self.reference_stress[0] = c0 * self.cl.poisson_ratio * (F01**2.0)/2.0
        self.reference_stress[1] = c0 * (F01**2.0)/2.0
        self.reference_stress[2] = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * 2.0) * F01
        return self.reference_stress

class SimpleShearElasticPlaneStressUncoupledShear2D(SimpleShearDeformation):
    def __init__(self, deltaDef):
        SimpleShearDeformation.__init__(self, deltaDef)
        self.cl = ElasticPlaneStressUncoupledShear2D()

    def get_reference_stress(self, i):
        c0 = self.cl.young_modulus / (1.0 - self.cl.poisson_ratio**2)
        F01 = self.get_deformation_gradientF(i)[0,1]
        absGamma12 = abs(F01)

        self.reference_stress[0] = c0 * self.cl.poisson_ratio * (F01**2.0)/2.0
        self.reference_stress[1] = c0 * (F01**2.0)/2.0
        self.reference_stress[2] = (self.cl.shear_modulus + self.cl.shear_modulus_gamma12 * absGamma12 + self.cl.shear_modulus_gamma12_2 * absGamma12**2 + self.cl.shear_modulus_gamma12_3 * absGamma12**3 + self.cl.shear_modulus_gamma12_4 * absGamma12**4)* F01
        return self.reference_stress

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

    def get_reference_stress(self, i):
        lame_lambda = (self.cl.young_modulus * self.cl.poisson_ratio) / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        x1beta = self.x1beta
        x2beta = self.x2beta
        x3beta = self.x3beta

        self.reference_stress[0]= math.cos(x3beta * i)*(x2beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) + (x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + x2beta*lame_mu*math.cos(x3beta * i) + x1beta*lame_mu*math.sin(x3beta * i)) + math.sin(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 + x1beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)))
        self.reference_stress[1]= math.cos(x3beta * i)*(x1beta*lame_mu*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) + (x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + x1beta*lame_mu*math.cos(x3beta * i) - x2beta*lame_mu*math.sin(x3beta * i)) + math.sin(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 - x2beta*lame_mu*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)))
        self.reference_stress[2]=(lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)
        self.reference_stress[3]= math.sin(x3beta * i)*(x2beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) - math.cos(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 + x1beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i))) - (x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + x2beta*lame_mu*math.cos(x3beta * i) + x1beta*lame_mu*math.sin(x3beta * i))
        self.reference_stress[4]=(lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + x1beta*lame_mu*math.cos(x3beta * i) - x2beta*lame_mu*math.sin(x3beta * i)
        self.reference_stress[5]=- (lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) - x2beta*lame_mu*math.cos(x3beta * i) - x1beta*lame_mu*math.sin(x3beta * i)
        return self.reference_stress

class ShearPlusStrechHyperElastic3D(ShearPlusStrechDeformation):
    def __init__(self):
        ShearPlusStrechDeformation.__init__(self)
        self.cl = HyperElastic3D()

    def get_reference_stress(self, i):
        lame_mu = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        x1beta = self.x1beta
        x2beta = self.x2beta
        x3beta = self.x3beta

        self.reference_stress[0] = (x2beta * math.cos(i * x3beta) + x1beta * math.sin(i * x3beta))**2.0
        self.reference_stress[1] = (x1beta * math.cos(i * x3beta) - x2beta * math.sin(i * x3beta))**2.0
        self.reference_stress[3] = (x2beta * math.cos(i * x3beta) + x1beta * math.sin(i * x3beta)) * (- x1beta * math.cos(i * x3beta) + x2beta * math.sin(i * x3beta))
        self.reference_stress[4] = x1beta * math.cos(i * x3beta) - x2beta * math.sin(i * x3beta)
        self.reference_stress[5] = - x2beta * math.cos(i * x3beta) - x1beta * math.sin(i * x3beta)
        self.reference_stress *= lame_mu
        return self.reference_stress

class ShearPlusStrechLinearElastic3D(ShearPlusStrechDeformation):
    def __init__(self):
        ShearPlusStrechDeformation.__init__(self)
        self.cl = LinearElastic3D()

    def get_reference_stress(self, i):
        c0 = self.cl.young_modulus / ((1.0 + self.cl.poisson_ratio) * (1.0 - 2.0 * self.cl.poisson_ratio))
        c1 = self.cl.young_modulus / (2.0 * (1.0 + self.cl.poisson_ratio))
        x1beta = self.x1beta
        x2beta = self.x2beta

        self.reference_stress[0] = c0 * self.cl.poisson_ratio * (x1beta**2.0 + x2beta**2.0) / 2.0
        self.reference_stress[1] = self.reference_stress[0]
        self.reference_stress[2] = c0 * (1.0 - self.cl.poisson_ratio) * (x1beta**2.0 + x2beta**2.0) / 2.0
        self.reference_stress[4] = x2beta * c1
        self.reference_stress[5] = -c1 * x1beta
        return self.reference_stress

class DeformationLinearJ2Plasticity(Deformation):
    def __init__(self):
        Deformation.__init__(self)
        self.nr_timesteps = 10

    def set_deformation(self, cl_params, i):
        self.strain = (i+1)/ self.nr_timesteps * self.initial_strain
        cl_params.SetStrainVector(self.strain)

class DeformationLinearJ2Plasticity3D(DeformationLinearJ2Plasticity):
    def __init__(self):
        DeformationLinearJ2Plasticity.__init__(self)
        self.cl = LinearJ2Plasticity3D()

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = 0.001
        self.initial_strain[1] = 0.001
        self.initial_strain[2] = 0.0
        self.initial_strain[3] = 0.001
        self.initial_strain[4] = 0.0
        self.initial_strain[5] = 0.001

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        r_stress[0][0] = 4.03846; r_stress[0][1] = 4.03846; r_stress[0][2] = 2.42308; r_stress[0][3] = 0.80769; r_stress[0][4] = 0.0; r_stress[0][5] = 0.80769
        r_stress[1][0] = 8.07692; r_stress[1][1] = 8.07692; r_stress[1][2] = 4.84615; r_stress[1][3] = 1.61538; r_stress[1][4] = 0.0; r_stress[1][5] = 1.61538
        r_stress[2][0] = 11.6595; r_stress[2][1] = 11.6595; r_stress[2][2] = 8.18099; r_stress[2][3] = 1.73926; r_stress[2][4] = 0.0; r_stress[2][5] = 1.73926
        r_stress[3][0] = 15.1595; r_stress[3][1] = 15.1595; r_stress[3][2] = 11.681 ; r_stress[3][3] = 1.73926; r_stress[3][4] = 0.0; r_stress[3][5] = 1.73926
        r_stress[4][0] = 18.6595; r_stress[4][1] = 18.6595; r_stress[4][2] = 15.181 ; r_stress[4][3] = 1.73926; r_stress[4][4] = 0.0; r_stress[4][5] = 1.73926
        r_stress[5][0] = 22.1595; r_stress[5][1] = 22.1595; r_stress[5][2] = 18.681 ; r_stress[5][3] = 1.73927; r_stress[5][4] = 0.0; r_stress[5][5] = 1.73927
        r_stress[6][0] = 25.6595; r_stress[6][1] = 25.6595; r_stress[6][2] = 22.181 ; r_stress[6][3] = 1.73927; r_stress[6][4] = 0.0; r_stress[6][5] = 1.73927
        r_stress[7][0] = 29.1595; r_stress[7][1] = 29.1595; r_stress[7][2] = 25.681 ; r_stress[7][3] = 1.73928; r_stress[7][4] = 0.0; r_stress[7][5] = 1.73928
        r_stress[8][0] = 32.6595; r_stress[8][1] = 32.6595; r_stress[8][2] = 29.181 ; r_stress[8][3] = 1.73928; r_stress[8][4] = 0.0; r_stress[8][5] = 1.73928
        r_stress[9][0] = 36.1595; r_stress[9][1] = 36.1595; r_stress[9][2] = 32.681; r_stress[9][3] = 1.73929; r_stress[9][4] = 0.0; r_stress[9][5] = 1.73929
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]

class DeformationLinearJ2PlasticityPlaneStrain2D(DeformationLinearJ2Plasticity):
    def __init__(self):
        DeformationLinearJ2Plasticity.__init__(self)
        self.cl = LinearJ2PlasticityPlaneStrain2D()

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = 0.001
        self.initial_strain[1] = 0.001
        self.initial_strain[2] = 0.0
        self.initial_strain[3] = 0.001

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        r_stress[0][0] = 4.03846; r_stress[0][1] = 4.03846; r_stress[0][2] = 2.42308; r_stress[0][3] = 0.807692;
        r_stress[1][0] = 8.07692; r_stress[1][1] = 8.07692; r_stress[1][2] = 4.84615; r_stress[1][3] = 1.61538;
        r_stress[2][0] = 11.8859; r_stress[2][1] = 11.8859; r_stress[2][2] = 7.72826; r_stress[2][3] = 2.07881;
        r_stress[3][0] = 15.3859; r_stress[3][1] = 15.3859; r_stress[3][2] = 11.2283; r_stress[3][3] = 2.07881;
        r_stress[4][0] = 18.8859; r_stress[4][1] = 18.8859; r_stress[4][2] = 14.7282; r_stress[4][3] = 2.07882;
        r_stress[5][0] = 22.3859; r_stress[5][1] = 22.3859; r_stress[5][2] = 18.2282; r_stress[5][3] = 2.07882;
        r_stress[6][0] = 25.8859; r_stress[6][1] = 25.8859; r_stress[6][2] = 21.7282; r_stress[6][3] = 2.07882;
        r_stress[7][0] = 29.3859; r_stress[7][1] = 29.3859; r_stress[7][2] = 25.2282; r_stress[7][3] = 2.07883;
        r_stress[8][0] = 32.8859; r_stress[8][1] = 32.8859; r_stress[8][2] = 28.7282; r_stress[8][3] = 2.07883;
        r_stress[9][0] = 36.3859; r_stress[9][1] = 36.3859; r_stress[9][2] = 32.2282; r_stress[9][3] = 2.07884;
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]

class DeformationLinearIsotropicDamage(Deformation):
    def __init__(self):
        Deformation.__init__(self)
        self.nr_timesteps = 10

    def set_deformation(self, cl_params, i):
        self.strain = (i+1)/ self.nr_timesteps * self.initial_strain
        cl_params.SetStrainVector(self.strain)

class DeformationLinearIsotropicDamage3D(DeformationLinearIsotropicDamage):
    def __init__(self):
        DeformationLinearIsotropicDamage.__init__(self)
        self.cl = LinearIsotropicDamage3D()

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = 0.001
        self.initial_strain[1] = 0.001
        self.initial_strain[2] = 0.0
        self.initial_strain[3] = 0.001
        self.initial_strain[4] = 0.0
        self.initial_strain[5] = 0.001

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        r_stress[0][0] = 0.57692; r_stress[0][1] = 0.57692; r_stress[0][2] = 0.34615; r_stress[0][3] = 0.11538; r_stress[0][4] = 0.0; r_stress[0][5] = 0.11538;
        r_stress[1][0] = 1.15384; r_stress[1][1] = 1.15384; r_stress[1][2] = 0.69231; r_stress[1][3] = 0.23077; r_stress[1][4] = 0.0; r_stress[1][5] = 0.23077;
        r_stress[2][0] = 1.73076; r_stress[2][1] = 1.73076; r_stress[2][2] = 1.03850; r_stress[2][3] = 0.34615; r_stress[2][4] = 0.0; r_stress[2][5] = 0.34615;
        r_stress[3][0] = 1.94550; r_stress[3][1] = 1.94550; r_stress[3][2] = 1.16730; r_stress[3][3] = 0.38910; r_stress[3][4] = 0.0; r_stress[3][5] = 0.38910;
        r_stress[4][0] = 2.11858; r_stress[4][1] = 2.11858; r_stress[4][2] = 1.27120; r_stress[4][3] = 0.42372; r_stress[4][4] = 0.0; r_stress[4][5] = 0.42372;
        r_stress[5][0] = 2.29166; r_stress[5][1] = 2.29166; r_stress[5][2] = 1.37500; r_stress[5][3] = 0.45833; r_stress[5][4] = 0.0; r_stress[5][5] = 0.45833;
        r_stress[6][0] = 2.46473; r_stress[6][1] = 2.46473; r_stress[6][2] = 1.47880; r_stress[6][3] = 0.49295; r_stress[6][4] = 0.0; r_stress[6][5] = 0.49295;
        r_stress[7][0] = 2.63781; r_stress[7][1] = 2.63781; r_stress[7][2] = 1.58270; r_stress[7][3] = 0.52756; r_stress[7][4] = 0.0; r_stress[7][5] = 0.52756;
        r_stress[8][0] = 2.68543; r_stress[8][1] = 2.68543; r_stress[8][2] = 1.61130; r_stress[8][3] = 0.53709; r_stress[8][4] = 0.0; r_stress[8][5] = 0.53709;
        r_stress[9][0] = 2.68543; r_stress[9][1] = 2.68543; r_stress[9][2] = 1.61130; r_stress[9][3] = 0.53709; r_stress[9][4] = 0.0; r_stress[9][5] = 0.53709;
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]

class DeformationLinearIsotropicDamagePlaneStrain2D(DeformationLinearIsotropicDamage):
    def __init__(self):
        DeformationLinearIsotropicDamage.__init__(self)
        self.cl = LinearIsotropicDamagePlaneStrain2D()

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = 0.001
        self.initial_strain[1] = 0.001
        self.initial_strain[2] = 0.001

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        r_stress[0][0] = 0.57692; r_stress[0][1] = 0.57692; r_stress[0][2] = 0.11538;
        r_stress[1][0] = 1.15384; r_stress[1][1] = 1.15384; r_stress[1][2] = 0.23077;
        r_stress[2][0] = 1.73076; r_stress[2][1] = 1.73076; r_stress[2][2] = 0.34615;
        r_stress[3][0] = 2.00123; r_stress[3][1] = 2.00123; r_stress[3][2] = 0.40025;
        r_stress[4][0] = 2.17431; r_stress[4][1] = 2.17431; r_stress[4][2] = 0.43486;
        r_stress[5][0] = 2.34738; r_stress[5][1] = 2.34738; r_stress[5][2] = 0.46948;
        r_stress[6][0] = 2.52046; r_stress[6][1] = 2.52046; r_stress[6][2] = 0.50409;
        r_stress[7][0] = 2.69354; r_stress[7][1] = 2.69354; r_stress[7][2] = 0.53871;
        r_stress[8][0] = 2.80484; r_stress[8][1] = 2.80484; r_stress[8][2] = 0.56097;
        r_stress[9][0] = 2.80484; r_stress[9][1] = 2.80484; r_stress[9][2] = 0.56097;
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]

# todo -****************************
class DeformationSmallStrainIsotropicPlasticity3D(DeformationLinearIsotropicDamage):
    def __init__(self):
        DeformationLinearIsotropicDamage.__init__(self)
        self.cl = LinearIsotropicDamage3D()

    def initialize_reference_stress(self, strain_size):
        self.initial_strain = KratosMultiphysics.Vector(strain_size)
        self.initial_strain[0] = 0.001
        self.initial_strain[1] = 0.001
        self.initial_strain[2] = 0.0
        self.initial_strain[3] = 0.001
        self.initial_strain[4] = 0.0
        self.initial_strain[5] = 0.001

        r_stress = []
        for i in range(self.nr_timesteps):
            r_stress.append(KratosMultiphysics.Vector(strain_size))
        r_stress[0][0] = 0.57692; r_stress[0][1] = 0.57692; r_stress[0][2] = 0.34615; r_stress[0][3] = 0.11538; r_stress[0][4] = 0.0; r_stress[0][5] = 0.11538;
        r_stress[1][0] = 1.15384; r_stress[1][1] = 1.15384; r_stress[1][2] = 0.69231; r_stress[1][3] = 0.23077; r_stress[1][4] = 0.0; r_stress[1][5] = 0.23077;
        r_stress[2][0] = 1.73076; r_stress[2][1] = 1.73076; r_stress[2][2] = 1.03850; r_stress[2][3] = 0.34615; r_stress[2][4] = 0.0; r_stress[2][5] = 0.34615;
        r_stress[3][0] = 1.94550; r_stress[3][1] = 1.94550; r_stress[3][2] = 1.16730; r_stress[3][3] = 0.38910; r_stress[3][4] = 0.0; r_stress[3][5] = 0.38910;
        r_stress[4][0] = 2.11858; r_stress[4][1] = 2.11858; r_stress[4][2] = 1.27120; r_stress[4][3] = 0.42372; r_stress[4][4] = 0.0; r_stress[4][5] = 0.42372;
        r_stress[5][0] = 2.29166; r_stress[5][1] = 2.29166; r_stress[5][2] = 1.37500; r_stress[5][3] = 0.45833; r_stress[5][4] = 0.0; r_stress[5][5] = 0.45833;
        r_stress[6][0] = 2.46473; r_stress[6][1] = 2.46473; r_stress[6][2] = 1.47880; r_stress[6][3] = 0.49295; r_stress[6][4] = 0.0; r_stress[6][5] = 0.49295;
        r_stress[7][0] = 2.63781; r_stress[7][1] = 2.63781; r_stress[7][2] = 1.58270; r_stress[7][3] = 0.52756; r_stress[7][4] = 0.0; r_stress[7][5] = 0.52756;
        r_stress[8][0] = 2.68543; r_stress[8][1] = 2.68543; r_stress[8][2] = 1.61130; r_stress[8][3] = 0.53709; r_stress[8][4] = 0.0; r_stress[8][5] = 0.53709;
        r_stress[9][0] = 2.68543; r_stress[9][1] = 2.68543; r_stress[9][2] = 1.61130; r_stress[9][3] = 0.53709; r_stress[9][4] = 0.0; r_stress[9][5] = 0.53709;
        self.reference_stress = r_stress

    def get_reference_stress(self, i):
        return self.reference_stress[i]
# todo -****************************


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

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.KirchhoffSaintVenant3DLaw()

class HyperElastic3D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 3

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.HyperElastic3DLaw()

class LinearElastic3D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 3

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.LinearElastic3DLaw()

class LinearElasticPlaneStress2D(LinearElastic):
    def __init__(self):
        LinearElastic.__init__(self)
        self.dim = 2

    @staticmethod
    def create_constitutive_Law():
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

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.ElasticPlaneStressUncoupledShear2DLaw()

class LinearJ2Plasticity(LinearElastic):
    def __init__(self):
        self.young_modulus = 21000
        self.poisson_ratio = 0.3
        self.yield_stress = 5.5
        self.reference_hardening_modulus = 1.0
        self.isotropic_hardening_modulus = 0.12924
        self.infinity_hardening_modulus = 0.0
        self.hardening_exponent = 1.0

    def create_properties(self, model_part):
        properties = LinearElastic.create_properties(self, model_part)
        properties.SetValue(KratosMultiphysics.YIELD_STRESS, self.yield_stress)
        properties.SetValue(KratosMultiphysics.REFERENCE_HARDENING_MODULUS, self.reference_hardening_modulus)
        properties.SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, self.isotropic_hardening_modulus)
        properties.SetValue(KratosMultiphysics.INFINITY_HARDENING_MODULUS, self.infinity_hardening_modulus)
        properties.SetValue(KratosMultiphysics.HARDENING_EXPONENT, self.hardening_exponent)
        return properties

class LinearJ2Plasticity3D(LinearJ2Plasticity):
    def __init__(self):
        LinearJ2Plasticity.__init__(self)
        self.dim = 3

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.LinearJ2Plasticity3DLaw()

class LinearJ2PlasticityPlaneStrain2D(LinearJ2Plasticity):
    def __init__(self):
        LinearJ2Plasticity.__init__(self)
        self.dim = 2

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.LinearJ2PlasticityPlaneStrain2DLaw()

class LinearIsotropicDamage(LinearElastic):
    def __init__(self):
        self.young_modulus = 3000
        self.poisson_ratio = 0.3
        self.yield_stress = 2.0
        self.infinity_yield_stress = 3.0
        self.isotropic_hardening_modulus = 0.3

    def create_properties(self, model_part):
        properties = LinearElastic.create_properties(self, model_part)
        properties.SetValue(KratosMultiphysics.YIELD_STRESS, self.yield_stress)
        properties.SetValue(StructuralMechanicsApplication.INFINITY_YIELD_STRESS, self.infinity_yield_stress)
        properties.SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, self.isotropic_hardening_modulus)
        return properties

class LinearIsotropicDamage3D(LinearIsotropicDamage):
    def __init__(self):
        LinearIsotropicDamage.__init__(self)
        self.dim = 3

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.LinearIsotropicDamage3DLaw()

class LinearIsotropicDamagePlaneStrain2D(LinearIsotropicDamage):
    def __init__(self):
        LinearIsotropicDamage.__init__(self)
        self.dim = 2

    @staticmethod
    def create_constitutive_Law():
        return StructuralMechanicsApplication.LinearIsotropicDamagePlaneStrain2DLaw()

if __name__ == '__main__':
    KratosUnittest.main()
