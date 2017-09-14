from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest

import math

class TestConstitutiveLaw(KratosUnittest.TestCase):

    def setUp(self):
        pass
    
    def test_Uniaxial_HyperElastic_3D(self):
        nnodes = 4
        dim = 3

        # Define a model part and create new nodes
        model_part = KratosMultiphysics.ModelPart("test")
        node1 = model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        node2 = model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        node3 = model_part.CreateNewNode(3, 0.0, 1.0, 0.0)
        node4 = model_part.CreateNewNode(4, 0.0, 0.0, 1.0)

        # Material properties
        prop_id = 0
        properties = model_part.Properties[prop_id]
        young_modulus = 200e9
        poisson_ratio = 0.3
        properties.SetValue(KratosMultiphysics.YOUNG_MODULUS, young_modulus)
        properties.SetValue(KratosMultiphysics.POISSON_RATIO, poisson_ratio)

        # Allocate a geometry
        geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node3, node4)

        N = KratosMultiphysics.Vector(4)
        DN_DX = KratosMultiphysics.Matrix(4,3)

        # Construct a constitutive law 
        cl = StructuralMechanicsApplication.HyperElastic3DLaw()
        cl.Check(properties, geom, model_part.ProcessInfo)

        if(cl.WorkingSpaceDimension() != dim):
            raise Exception("Mismatch between the WorkingSpaceDimension of the Constitutive Law and the dimension of the space in which the test is performed")

        ## Set the parameters to be employed
        #note that here i am adding them all to check that this does not fail
        cl_options = KratosMultiphysics.Flags()
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, False)
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRAIN_ENERGY, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE, False)

        # From here below it should be an otput not an input
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.FINITE_STRAINS, True) 
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.INFINITESIMAL_STRAINS, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.PLANE_STRAIN_LAW, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.PLANE_STRESS_LAW, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.AXISYMMETRIC_LAW, False)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.U_P_LAW, False)
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ISOTROPIC, True)
        #cl_options.Set(KratosMultiphysics.ConstitutiveLaw.ANISOTROPIC, False)

        F = KratosMultiphysics.Matrix(3,3)
        F[0,0] = 1.0; F[0,1] = 0.0; F[0,2] = 0.0;
        F[1,0] = 0.0; F[1,1] = 1.0; F[1,2] = 0.0;
        F[2,0] = 0.0; F[2,1] = 0.0; F[2,2] = 1.0;
        detF = 1.0

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())

        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

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

        #print("The Material Response PK2")
        #cl.CalculateMaterialResponsePK2(cl_params)
        #print("Stress = ", cl_params.GetStressVector())
        #print("Strain = ", cl_params.GetStrainVector())
        #print("C      = ", cl_params.GetConstitutiveMatrix())

        #cl.FinalizeMaterialResponsePK2(cl_params)
        #cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

        #print("\n The Material Response Kirchhoff")
        #cl.CalculateMaterialResponseKirchhoff(cl_params)
        #print("Stress = ", cl_params.GetStressVector())
        #print("Strain = ", cl_params.GetStrainVector())
        #print("C      = ", cl_params.GetConstitutiveMatrix())

        #cl.FinalizeMaterialResponseKirchhoff(cl_params)
        #cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)

        #print("\n The Material Response Cauchy")
        #cl.CalculateMaterialResponseCauchy(cl_params)
        #print("Stress = ", cl_params.GetStressVector())
        #print("Strain = ", cl_params.GetStrainVector())
        #print("C      = ", cl_params.GetConstitutiveMatrix())
        
        #cl.FinalizeMaterialResponseCauchy(cl_params)
        #cl.FinalizeSolutionStep(properties, geom, N, model_part.ProcessInfo)
        
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
        
if __name__ == '__main__':
    KratosUnittest.main()

