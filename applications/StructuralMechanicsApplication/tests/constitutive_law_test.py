from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestConstitutiveLaw(KratosUnittest.TestCase):

    def setUp(self):
        pass
    
    def _create_geometry(self, model_part, dim, nnodes):
        
        # Create new nodes
        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
            
        if (dim == 2):
            if (nnodes == 3):
                node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
                
                # Allocate a geometry
                geom = KratosMultiphysics.Triangle2D3(node1,node2,node3)
        else:
            if (nnodes == 4):
                node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
                node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)
                
                # Allocate a geometry
                geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node3,node4)
                
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

