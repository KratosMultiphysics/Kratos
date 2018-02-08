from __future__ import print_function, absolute_import, division
import KratosMultiphysics 

import KratosMultiphysics.FluidDynamicsApplication 
import KratosMultiphysics.KratosUnittest as KratosUnittest


class TestTaiteEquation(KratosUnittest.TestCase):
    def setUp(self):
        pass
        
    def test_bulk_modulus_calculation(self):

        ######################### general parameters
        nnodes = 4
        dim = 3

        #define a model part and create new nodes
        model_part = KratosMultiphysics.ModelPart("test")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        
        node1 = model_part.CreateNewNode(1,0.0,0.0,0.0)
        node2 = model_part.CreateNewNode(2,1.0,0.0,0.0)
        node3 = model_part.CreateNewNode(3,0.0,1.0,0.0)
        node4 = model_part.CreateNewNode(4,0.0,0.0,1.0)
        
        #allocate a geometry
        geom = KratosMultiphysics.Tetrahedra3D4(node1,node2,node3,node4)
        print(geom)

        N = KratosMultiphysics.Vector(4)
        DN_DX = KratosMultiphysics.Matrix(4,3)
        for i in range(4):
            N[i] = 1.0/4.0

        #material properties
        prop_id = 0
        properties = model_part.Properties[prop_id]
        properties.SetValue(KratosMultiphysics.DYNAMIC_VISCOSITY, 1.0)
        properties.SetValue(KratosMultiphysics.FluidDynamicsApplication.REGULARIZATION_COEFFICIENT, 1.0)
        
        #b1m = 0.001295        
        #b2m = 8.588E-07   
        #b3m = 77200000
        #b4m = 0.003487      
        #b1s =  0.001242       
        #b2s = 9.153E-07  
        #b3s =  63200000   
        #b4s = 0.006881   
        #b5 =     403.0
        #b6 = 0.00000015  
        #b7  = 0.00002616 
        #b8  = 0.0714                     
        #b9  = 3.355E-09   
        
        #Coefficients by Ravi
        b1m = 1.01e-3   
        b2m = 6.4677e-7   
        b3m = 1.91e8
        b4m = 4.7537e-3   
        b1s =  1.01e-3     
        b2s = 2.0123e-7
        b3s =  2.9835e8   
        b4s = 2.0032e-3  
        b5 =    408.12
        b6 = 3.8307e-7 
        b7  = 0 
        b8  = 0                  
        b9  = 0 
        
        properties.SetValue(KratosMultiphysics.TAIT_PARAMETERS_MOLTEN_STATE, [0.0,b1m, b2m, b3m, b4m, b5, b6, b7, b8, b9 ])
        properties.SetValue(KratosMultiphysics.TAIT_PARAMETERS_SOLID_STATE,  [0.0,b1s, b2s, b3s, b4s, b5, b6, b7, b8, b9 ])
        

        


        ######################################## here we choose the constitutive law #########################
        #construct a constitutive law 
        cl = KratosMultiphysics.FluidDynamicsApplication.Bingham3DLaw()
        cl.Check( properties, geom, model_part.ProcessInfo )

        if(cl.WorkingSpaceDimension() != dim):
            raise Exception( "mismatch between the WorkingSpaceDimension of the Constitutive Law and the dimension of the space in which the test is performed")

        ##set the parameters to be employed
        #note that here i am adding them all to check that this does not fail
        cl_options = KratosMultiphysics.Flags()
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, False)
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_STRESS, True)
        cl_options.Set(KratosMultiphysics.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        #cl_options.Set(ConstitutiveLaw.COMPUTE_STRAIN_ENERGY, False)
        #cl_options.Set(ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY, False)
        #cl_options.Set(ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY, False)
        #cl_options.Set(ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE, False)

        ##from here below it should be an otput not an input
        #cl_options.Set(ConstitutiveLaw.FINITE_STRAINS, False) 
        #cl_options.Set(ConstitutiveLaw.INFINITESIMAL_STRAINS, False)
        #cl_options.Set(cl_params.SetDeformationGradientF( F )
        #cl_options.Set(ConstitutiveLaw.PLANE_STRESS_LAW, False)
        #cl_options.Set(ConstitutiveLaw.AXISYMMETRIC_LAW, False)
        #cl_options.Set(ConstitutiveLaw.U_P_LAW, False)
        #cl_options.Set(ConstitutiveLaw.ISOTROPIC, False)
        #cl_options.Set(ConstitutiveLaw.ANISOTROPIC, False)

        stress_vector = KratosMultiphysics.Vector(cl.GetStrainSize() )
        strain_vector = KratosMultiphysics.Vector(cl.GetStrainSize())
        constitutive_matrix = KratosMultiphysics.Matrix(cl.GetStrainSize(),cl.GetStrainSize())

        #setting the parameters - note that a constitutive law may not need them all!
        cl_params = KratosMultiphysics.ConstitutiveLawParameters()
        cl_params.SetOptions( cl_options )
        cl_params.SetStrainVector( strain_vector )
        cl_params.SetStressVector( stress_vector )
        cl_params.SetConstitutiveMatrix( constitutive_matrix )
        cl_params.SetShapeFunctionsValues( N )
        cl_params.SetShapeFunctionsDerivatives( DN_DX )
        cl_params.SetProcessInfo( model_part.ProcessInfo )
        cl_params.SetMaterialProperties( properties )
        cl_params.SetElementGeometry(geom)
        
        F = KratosMultiphysics.Matrix(3,3)
        detF = 1.0
        cl_params.SetDeformationGradientF( F )
        cl_params.SetDeterminantF( detF )

        ##do all sort of checks
        cl_params.CheckAllParameters() #can not use this until the geometry is correctly exported to python
        cl_params.CheckMechanicalVariables()
        cl_params.CheckShapeFunctions()

        cl.CalculateMaterialResponseCauchy( cl_params )
        print( "stress = ", cl_params.GetStressVector() )
        print( "strain = ", cl_params.GetStrainVector() )
        print( "C      = ", cl_params.GetConstitutiveMatrix() )

        #cl.FinalizeMaterialResponseCauchy( cl_params )
        cl.FinalizeSolutionStep( properties, geom, N, model_part.ProcessInfo )
        
        p0 = 0 #200000
        for node in model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 400.0)
            node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, p0)
            
        outfile = open("out.csv",'w')
        
        for step in range(1000):
            p = p0 + step*10000
            for node in model_part.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, p)
                
            #k = cl.GetValue(    KratosMultiphysics.BULK_MODULUS, k)
            k = cl.CalculateValue(cl_params, KratosMultiphysics.BULK_MODULUS)
            outfile.write( str(p) + " " + str(k) + "\n")
        
        outfile.close()

if __name__ == '__main__':
    KratosUnittest.main()
