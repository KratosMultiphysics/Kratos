from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestConstitutiveLaw(KratosUnittest.TestCase):

    def setUp(self):
        pass
    
    def _create_geometry(self, model_part, dim, nnodes):
        
        # Create nodes
        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node3 = model_part.CreateNewNode(3, 1.0, 1.0, 0.0)
        node4 = model_part.CreateNewNode(4, 0.0, 1.0, 0.0)
        node5 = model_part.CreateNewNode(5, 0.0, 0.0, 1.0)
        node6 = model_part.CreateNewNode(6, 1.0, 0.0, 1.0)
        node7 = model_part.CreateNewNode(7, 1.0, 1.0, 1.0)
        node8 = model_part.CreateNewNode(8, 0.0, 1.0, 1.0)

        if (dim == 2):
            if (nnodes == 3):
                geom = KratosMultiphysics.Triangle2D3(node1,node2,node4)
            elif (nnodes == 4):
                geom = KratosMultiphysics.Quadrilateral2D4(node1,node2,node3, node4)
        else:
            if (nnodes == 4):
                geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node4,node5)
            elif (nnodes == 8):
                geom = KratosMultiphysics.Hexahedra3D8(node1, node2, node3, node4, node5, node6, node7, node8)
                
        return geom
    
    def _create_properties(self, model_part, young_modulus, poisson_ratio):
        
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

    def _create_F(self):
        F = KratosMultiphysics.Matrix(3,3)
        F[0,0] = 1.0; F[0,1] = 0.0; F[0,2] = 0.0;
        F[1,0] = 0.0; F[1,1] = 1.0; F[1,2] = 0.0;
        F[2,0] = 0.0; F[2,1] = 0.0; F[2,2] = 1.0;
        return F
    
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

    def test_J2_Plasticity_PlaneStrain_2D(self):
        nnodes = 4
        dim = 2

        # Define a model and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        prop_id = 0
        properties = model_part.Properties[prop_id]
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, 21000)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        properties.SetValue(KratosMultiphysics.YIELD_STRESS, 5.5)
        properties.SetValue(KratosMultiphysics.REFERENCE_HARDENING_MODULUS, 1.0)
        properties.SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, 0.12924)
        properties.SetValue(KratosMultiphysics.INFINITY_HARDENING_MODULUS, 0.0)
        properties.SetValue(KratosMultiphysics.HARDENING_EXPONENT, 1.0)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law
        cl = StructuralMechanicsApplication.LinearJ2PlasticityPlaneStrain2DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                        'COMPUTE_STRESS': True,
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        }
        cl_options = self._set_cl_options(dict_options)
        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        F = self._create_F()
        detF = 1.0
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector,
                                            stress_vector, constitutive_matrix,
                                            N, DN_DX, model_part, properties, geom)
        cl.InitializeMaterial(properties, geom, N)

        # Check the results
        nr_timesteps = 10
        rstress = []
        for i in range(nr_timesteps):
            rstress.append(KratosMultiphysics.Vector(cl.GetStrainSize()))
        rstress[0][0] = 4.03846; rstress[0][1] = 4.03846; rstress[0][2] = 2.42308; rstress[0][3] = 0.807692;
        rstress[1][0] = 8.07692; rstress[1][1] = 8.07692; rstress[1][2] = 4.84615; rstress[1][3] = 1.61538;
        rstress[2][0] = 11.8859; rstress[2][1] = 11.8859; rstress[2][2] = 7.72826; rstress[2][3] = 2.07881;
        rstress[3][0] = 15.3859; rstress[3][1] = 15.3859; rstress[3][2] = 11.2283; rstress[3][3] = 2.07881;
        rstress[4][0] = 18.8859; rstress[4][1] = 18.8859; rstress[4][2] = 14.7282; rstress[4][3] = 2.07882;
        rstress[5][0] = 22.3859; rstress[5][1] = 22.3859; rstress[5][2] = 18.2282; rstress[5][3] = 2.07882;
        rstress[6][0] = 25.8859; rstress[6][1] = 25.8859; rstress[6][2] = 21.7282; rstress[6][3] = 2.07882;
        rstress[7][0] = 29.3859; rstress[7][1] = 29.3859; rstress[7][2] = 25.2282; rstress[7][3] = 2.07883;
        rstress[8][0] = 32.8859; rstress[8][1] = 32.8859; rstress[8][2] = 28.7282; rstress[8][3] = 2.07883;
        rstress[9][0] = 36.3859; rstress[9][1] = 36.3859; rstress[9][2] = 32.2282; rstress[9][3] = 2.07884;

        initial_strain = KratosMultiphysics.Vector(cl.GetStrainSize())
        initial_strain[0] = 0.001
        initial_strain[1] = 0.001
        initial_strain[2] = 0.0
        initial_strain[3] = 0.001

        t = dt = 1. / nr_timesteps
        c = 0
        while(t <= 1. + dt /10.):
            strain = t * initial_strain
            cl_params.SetStrainVector(strain)

            # Chauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
            stress = cl_params.GetStressVector()
            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(rstress[c][j], stress[j], 4)
            t += dt
            c += 1

    def test_J2_Plasticity_3D(self):
        nnodes = 8
        dim = 3

        # Define a model and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        prop_id = 0
        properties = model_part.Properties[prop_id]
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, 21000)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, 0.3)
        properties.SetValue(KratosMultiphysics.YIELD_STRESS, 5.5)
        properties.SetValue(KratosMultiphysics.REFERENCE_HARDENING_MODULUS, 1.0)
        properties.SetValue(KratosMultiphysics.ISOTROPIC_HARDENING_MODULUS, 0.12924)
        properties.SetValue(KratosMultiphysics.INFINITY_HARDENING_MODULUS, 0.0)
        properties.SetValue(KratosMultiphysics.HARDENING_EXPONENT, 1.0)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law
        cl = StructuralMechanicsApplication.LinearJ2Plasticity3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                        'COMPUTE_STRESS': True,
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        }
        cl_options = self._set_cl_options(dict_options)
        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        F = self._create_F()
        detF = 1.0
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector,
                                            stress_vector, constitutive_matrix,
                                            N, DN_DX, model_part, properties, geom)
        cl.InitializeMaterial(properties, geom, N)

        # Check the results
        nr_timesteps = 10
        rstress = []
        for i in range(nr_timesteps):
            rstress.append(KratosMultiphysics.Vector(cl.GetStrainSize()))
        rstress[0][0] = 4.03846; rstress[0][1] = 4.03846; rstress[0][2] = 2.42308; rstress[0][3] = 0.80769; rstress[0][4] = 0.0; rstress[0][5] = 0.80769;
        rstress[1][0] = 8.07692; rstress[1][1] = 8.07692; rstress[1][2] = 4.84615; rstress[1][3] = 1.61538; rstress[1][4] = 0.0; rstress[1][5] = 1.61538;
        rstress[2][0] = 11.6595; rstress[2][1] = 11.6595; rstress[2][2] = 8.18099; rstress[2][3] = 1.73926; rstress[2][4] = 0.0; rstress[2][5] = 1.73926;
        rstress[3][0] = 15.1595; rstress[3][1] = 15.1595; rstress[3][2] = 11.681 ; rstress[3][3] = 1.73926; rstress[3][4] = 0.0; rstress[3][5] = 1.73926;
        rstress[4][0] = 18.6595; rstress[4][1] = 18.6595; rstress[4][2] = 15.181 ; rstress[4][3] = 1.73926; rstress[4][4] = 0.0; rstress[4][5] = 1.73926;
        rstress[5][0] = 22.1595; rstress[5][1] = 22.1595; rstress[5][2] = 18.681 ; rstress[5][3] = 1.73927; rstress[5][4] = 0.0; rstress[5][5] = 1.73927;
        rstress[6][0] = 25.6595; rstress[6][1] = 25.6595; rstress[6][2] = 22.181 ; rstress[6][3] = 1.73927; rstress[6][4] = 0.0; rstress[6][5] = 1.73927;
        rstress[7][0] = 29.1595; rstress[7][1] = 29.1595; rstress[7][2] = 25.681 ; rstress[7][3] = 1.73928; rstress[7][4] = 0.0; rstress[7][5] = 1.73928;
        rstress[8][0] = 32.6595; rstress[8][1] = 32.6595; rstress[8][2] = 29.181 ; rstress[8][3] = 1.73928; rstress[8][4] = 0.0; rstress[8][5] = 1.73928;
        rstress[9][0] = 36.1595; rstress[9][1] = 36.1595; rstress[9][2] = 32.681; rstress[9][3] = 1.73929; rstress[9][4] = 0.0; rstress[9][5] = 1.73929;

        initial_strain = KratosMultiphysics.Vector(cl.GetStrainSize())
        initial_strain[0] = 0.001
        initial_strain[1] = 0.001
        initial_strain[2] = 0.0
        initial_strain[3] = 0.001
        initial_strain[4] = 0.0
        initial_strain[5] = 0.001

        t = dt = 1. / nr_timesteps
        c = 0
        while(t <= 1. + dt /10.):
            strain = t * initial_strain
            cl_params.SetStrainVector(strain)

            # Chauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
            stress = cl_params.GetStressVector()
            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(rstress[c][j], stress[j], 4)
            t += dt
            c += 1

    def test_Uniaxial_KirchhoffSaintVenant_3D(self):
        nnodes = 4
        dim = 3

        # Define a model  and geometry
        model_part = KratosMultiphysics.ModelPart("test")

        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties = self._create_properties(model_part, young_modulus, poisson_ratio)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law
        cl = StructuralMechanicsApplication.KirchhoffSaintVenant3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                        'COMPUTE_STRESS': True,
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = self._create_F()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)

        # Check the results
        lame_lambda = (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
        lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio))

        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0

        for i in range(100):
            F[0,0] = 1.0 + 0.05 * i
            detF = 1.0 + 0.05 * i
            cl_params.SetDeformationGradientF(F)
            cl_params.SetDeterminantF(detF)

            # Chauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
            reference_stress[0] =( (lame_lambda *  0.5  + lame_mu) * (detF ** 2.0 - 1.0)*(detF ** 2.0) ) / detF
            reference_stress[1] = 0.5*lame_lambda*(detF ** 2.0 - 1.0) / detF
            reference_stress[2] = reference_stress[1]

            stress = cl_params.GetStressVector()

            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2)

    def test_Shear_KirchhoffSaintVenant_3D(self):
        nnodes = 4
        dim = 3

        # Define a model  and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties = self._create_properties(model_part, young_modulus, poisson_ratio)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law
        cl = StructuralMechanicsApplication.KirchhoffSaintVenant3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                        'COMPUTE_STRESS': True,
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = self._create_F()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)

        # Check the results
        lame_lambda = (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
        lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio))

        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0

        for i in range(100):
            F[0,1] = 0.02 * i
            cl_params.SetDeformationGradientF(F)

            # Cauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

            reference_stress[0] = (0.5*lame_lambda + 2*lame_mu) * (F[0,1])**2.0 + (0.5*lame_lambda + lame_mu) * (F[0,1])**4.0
            reference_stress[1] = (0.5*lame_lambda + lame_mu) * (F[0,1])**2.0
            reference_stress[2] = 0.5*lame_lambda  * (F[0,1])**2.0
            reference_stress[3] = lame_mu * (F[0,1]) + (0.5*lame_lambda + lame_mu) * (F[0,1])**3.0

            stress = cl_params.GetStressVector()

            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2)

    def test_Shear_Plus_Strech_KirchhoffSaintVenant_3D(self):
        nnodes = 4
        dim = 3

        # Define a model  and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties = self._create_properties(model_part, young_modulus, poisson_ratio)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law
        cl = StructuralMechanicsApplication.KirchhoffSaintVenant3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False,
                        'COMPUTE_STRESS': True,
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = self._create_F()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)

        # Check the results
        lame_lambda = (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
        lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio))

        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0

        x1beta = 1.0
        x2beta = 1.0
        x3beta = math.pi/200
        for i in range(100):
            F[0,0] =  math.cos(x3beta * i)
            F[0,1] = -math.sin(x3beta * i)
            F[1,0] =  math.sin(x3beta * i)
            F[1,1] =  math.cos(x3beta * i)
            F[0,2] =  - x1beta * math.sin(x3beta * i) - x2beta * math.cos(x3beta * i)
            F[1,2] =    x1beta * math.cos(x3beta * i) - x2beta * math.sin(x3beta * i)

            cl_params.SetDeformationGradientF(F)

            # Cauchy
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

            reference_stress[0]= math.cos(x3beta * i)*(x2beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) + (x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + x2beta*lame_mu*math.cos(x3beta * i) + x1beta*lame_mu*math.sin(x3beta * i)) + math.sin(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 + x1beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)))
            reference_stress[1]= math.cos(x3beta * i)*(x1beta*lame_mu*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) + (x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + x1beta*lame_mu*math.cos(x3beta * i) - x2beta*lame_mu*math.sin(x3beta * i)) + math.sin(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 - x2beta*lame_mu*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)))
            reference_stress[2]=(lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)
            reference_stress[3]= math.sin(x3beta * i)*(x2beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + (lame_lambda*math.cos(x3beta * i)*(x1beta**2 + x2beta**2))/2) - math.cos(x3beta * i)*((lame_lambda*math.sin(x3beta * i)*(x1beta**2 + x2beta**2))/2 + x1beta*lame_mu*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i))) - (x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i))*((lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) + x2beta*lame_mu*math.cos(x3beta * i) + x1beta*lame_mu*math.sin(x3beta * i))
            reference_stress[4]=(lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x1beta*math.cos(x3beta * i) - x2beta*math.sin(x3beta * i)) + x1beta*lame_mu*math.cos(x3beta * i) - x2beta*lame_mu*math.sin(x3beta * i)
            reference_stress[5]=- (lame_lambda/2 + lame_mu)*(x1beta**2 + x2beta**2)*(x2beta*math.cos(x3beta * i) + x1beta*math.sin(x3beta * i)) - x2beta*lame_mu*math.cos(x3beta * i) - x1beta*lame_mu*math.sin(x3beta * i)

            stress = cl_params.GetStressVector()

            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2)

    def test_Uniaxial_HyperElastic_3D(self):
        nnodes = 4
        dim = 3

        # Define a model  and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties = self._create_properties(model_part, young_modulus, poisson_ratio)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law 
        cl = StructuralMechanicsApplication.HyperElastic3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False, 
                        'COMPUTE_STRESS': True, 
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = self._create_F()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)
        
        # Check the results
        lame_lambda = (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
        lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio))
        
        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0
        
        for i in range(100):
            F[0,0] = 1.0 + 0.2 * i
            detF = 1.0 + 0.2 * i
            cl_params.SetDeformationGradientF(F)
            cl_params.SetDeterminantF(detF)
            
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
        
            reference_stress[0] = (lame_lambda * math.log(detF) + lame_mu * (detF ** 2.0 - 1.0)) / detF
            reference_stress[1] = (lame_lambda * math.log(detF)) / detF
            reference_stress[2] = reference_stress[1]
        
            stress = cl_params.GetStressVector()
            
            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2) 
                
    def test_Shear_HyperElastic_3D(self):
        nnodes = 4
        dim = 3

        # Define a model  and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties = self._create_properties(model_part, young_modulus, poisson_ratio)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law 
        cl = StructuralMechanicsApplication.HyperElastic3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False, 
                        'COMPUTE_STRESS': True, 
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = self._create_F()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)
        
        # Check the results
        lame_lambda = (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
        lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio))
        
        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0
        
        for i in range(100):
            F[0,1] = 0.2 * i
            cl_params.SetDeformationGradientF(F)
            
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
        
            reference_stress[0] = lame_mu * (0.2 * i)**2.0
            reference_stress[3] = lame_mu * (0.2 * i)
        
            stress = cl_params.GetStressVector()
            
            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2)
                
    def test_Shear_Plus_Strech_HyperElastic_3D(self):
        nnodes = 4
        dim = 3

        # Define a model  and geometry
        model_part = KratosMultiphysics.ModelPart("test")
        geom = self._create_geometry(model_part, dim, nnodes)

        # Material properties
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties = self._create_properties(model_part, young_modulus, poisson_ratio)

        N = KratosMultiphysics.Vector(nnodes)
        DN_DX = KratosMultiphysics.Matrix(nnodes, dim)

        # Construct a constitutive law 
        cl = StructuralMechanicsApplication.HyperElastic3DLaw()
        self._cl_check(cl, properties, geom, model_part, dim)

        # Set the parameters to be employed
        dict_options = {'USE_ELEMENT_PROVIDED_STRAIN': False, 
                        'COMPUTE_STRESS': True, 
                        'COMPUTE_CONSTITUTIVE_TENSOR': True,
                        'FINITE_STRAINS': True,
                        'ISOTROPIC': True,
                        }
        cl_options = self._set_cl_options(dict_options)

        # Define deformation gradient
        F = self._create_F()
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        # Setting the parameters - note that a constitutive law may not need them all!
        cl_params = self._set_cl_parameters(cl_options, F, detF, strain_vector, stress_vector, constitutive_matrix, N, DN_DX, model_part, properties, geom)
        
        # Check the results
        lame_lambda = (young_modulus * poisson_ratio) / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio))
        lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio))
        
        reference_stress = KratosMultiphysics.Vector(cl.GetStrainSize())
        for i in range(cl.GetStrainSize()):
            reference_stress[i] = 0.0
        
        x1beta = 1.0 
        x2beta = 1.0 
        x3beta = math.pi/200
        for i in range(100):
            F[0,0] =  math.cos(x3beta * i)
            F[0,1] = -math.sin(x3beta * i)
            F[1,0] =  math.sin(x3beta * i)
            F[1,1] =  math.cos(x3beta * i)
            F[0,2] =  - x1beta * math.sin(x3beta * i) - x2beta * math.cos(x3beta * i)
            F[1,2] =    x1beta * math.cos(x3beta * i) - x2beta * math.sin(x3beta * i)
                        
            cl_params.SetDeformationGradientF(F)
            
            cl.CalculateMaterialResponseCauchy(cl_params)
            cl.FinalizeMaterialResponseCauchy(cl_params)
            cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
        
            reference_stress[0] = (x2beta * math.cos(i * x3beta) + x1beta * math.sin(i * x3beta))**2.0
            reference_stress[1] = (x1beta * math.cos(i * x3beta) - x2beta * math.sin(i * x3beta))**2.0
            reference_stress[3] = (x2beta * math.cos(i * x3beta) + x1beta * math.sin(i * x3beta)) * (- x1beta * math.cos(i * x3beta) + x2beta * math.sin(i * x3beta))
            reference_stress[4] = x1beta * math.cos(i * x3beta) - x2beta * math.sin(i * x3beta)
            reference_stress[5] = - x2beta * math.cos(i * x3beta) - x1beta * math.sin(i * x3beta)
        
            reference_stress *= lame_mu
        
            stress = cl_params.GetStressVector()
            
            for j in range(cl.GetStrainSize()):
                self.assertAlmostEqual(reference_stress[j], stress[j], 2) 
        
if __name__ == '__main__':
    KratosUnittest.main()

