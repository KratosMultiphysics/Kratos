from KratosMultiphysics import *
import KratosMultiphysics.ConstitutiveModelsApplication as KratosMaterialModels

######################### general parameters
nnodes = 3
dim = 3

#define a model part and create new nodes
model = Model()
model_part = model.ModelPart("test")
node1 = model_part.CreateNewNode(1,0.0,0.0,0.0)
node2 = model_part.CreateNewNode(2,1.0,0.0,0.0)
node3 = model_part.CreateNewNode(3,0.0,1.0,0.0)

#material properties
prop_id = 0
properties = model_part.Properties[prop_id]
properties.SetValue(YOUNG_MODULUS, 200e9)
properties.SetValue(POISSON_RATIO, 0.3)

C10 = 200e9/(4*(1+0.3))
properties.SetValue(KratosMaterialModels.C10, C10)

#allocate a geometry
#a = PointerVector()
#a.append(node1)
#geom = Geometry(a)
#geom = Geometry()
geom = Triangle2D3(node1,node2,node3)
print(geom)

N = Vector(3)
DN_DX = Matrix(3,2)

######################################## here we choose the constitutive law #########################
#construct a constitutive law
elasticity_model = KratosMaterialModels.SaintVenantKirchhoffModel()
cl  = KratosMaterialModels.LargeStrain3DLaw(elasticity_model)

#plasticity_model = KratosMaterialModels.VonMisesNeoHookeanPlasticityModel()
#cl  = KratosMaterialModels.LargeStrain3DLaw(plasticity_model)

cl.Check( properties, geom, model_part.ProcessInfo )

if(cl.WorkingSpaceDimension() != dim):
    raise Exception( "mismatch between the WorkingSpaceDimension of the Constitutive Law and the dimension of the space in which the test is performed")

##set the parameters to be employed
#note that here i am adding them all to check that this does not fail
cl_options = Flags()
cl_options.Set(ConstitutiveLaw.USE_ELEMENT_PROVIDED_STRAIN, True)
cl_options.Set(ConstitutiveLaw.COMPUTE_STRESS, True)
cl_options.Set(ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
#cl_options.Set(ConstitutiveLaw.COMPUTE_STRAIN_ENERGY, False)
#cl_options.Set(ConstitutiveLaw.ISOCHORIC_TENSOR_ONLY, False)
#cl_options.Set(ConstitutiveLaw.VOLUMETRIC_TENSOR_ONLY, False)
#cl_options.Set(ConstitutiveLaw.FINALIZE_MATERIAL_RESPONSE, False)

##from here below it should be an output not an input
#cl_options.Set(ConstitutiveLaw.FINITE_STRAINS, False)
#cl_options.Set(ConstitutiveLaw.INFINITESIMAL_STRAINS, False)
#cl_options.Set(ConstitutiveLaw.PLANE_STRAIN_LAW, False)
#cl_options.Set(ConstitutiveLaw.PLANE_STRESS_LAW, False)
#cl_options.Set(ConstitutiveLaw.AXISYMMETRIC_LAW, False)
#cl_options.Set(ConstitutiveLaw.U_P_LAW, False)
#cl_options.Set(ConstitutiveLaw.ISOTROPIC, False)
#cl_options.Set(ConstitutiveLaw.ANISOTROPIC, False)

from math import sqrt

F = Matrix(3,3)
F[0,0] = 1.0; F[0,1] = 0.0; F[0,2] = 2.0;
F[1,0] = 0.0; F[1,1] = 0.9; F[1,2] = 0.0;
F[2,0] = 0.0; F[2,1] = 0.0; F[2,2] = 0.1;
detF = 0.09

stress_vector = Vector(cl.GetStrainSize())
strain_vector = Vector(cl.GetStrainSize())


constitutive_matrix = Matrix(cl.GetStrainSize(),cl.GetStrainSize())

#setting the parameters - note that a constitutive law may not need them all!
cl_params = ConstitutiveLawParameters()
cl_params.SetOptions( cl_options )
cl_params.SetDeformationGradientF( F )
cl_params.SetDeterminantF( detF )
cl_params.SetStrainVector( strain_vector )
cl_params.SetStressVector( stress_vector )
cl_params.SetConstitutiveMatrix( constitutive_matrix )
cl_params.SetShapeFunctionsValues( N )
cl_params.SetShapeFunctionsDerivatives( DN_DX )
cl_params.SetProcessInfo( model_part.ProcessInfo )
cl_params.SetMaterialProperties( properties )
cl_params.SetElementGeometry(geom)

##do all sort of checks
cl_params.CheckAllParameters() #can not use this until the geometry is correctly exported to python
cl_params.CheckMechanicalVariables()
cl_params.CheckShapeFunctions()

print("The Material Response PK2")
cl.CalculateMaterialResponsePK2( cl_params )
print( "stress = ", cl_params.GetStressVector() )
print( "strain = ", cl_params.GetStrainVector() )
print( "C      = ", cl_params.GetConstitutiveMatrix() )

#cl.FinalizeMaterialResponsePK2( cl_params )
cl.FinalizeSolutionStep( properties, geom, N, model_part.ProcessInfo )


print("\n The Material Response Kirchhoff")
cl.CalculateMaterialResponseKirchhoff( cl_params )
print( "stress = ", cl_params.GetStressVector() )
print( "strain = ", cl_params.GetStrainVector() )
print( "C      = ", cl_params.GetConstitutiveMatrix() )

cl.FinalizeMaterialResponseKirchhoff( cl_params )
cl.FinalizeSolutionStep( properties, geom, N, model_part.ProcessInfo )

print("\n The Material Response Cauchy")
cl.CalculateMaterialResponseCauchy( cl_params )
print( "stress = ", cl_params.GetStressVector() )
print( "strain = ", cl_params.GetStrainVector() )
print( "C      = ", cl_params.GetConstitutiveMatrix() )

cl.FinalizeMaterialResponseCauchy( cl_params )
cl.FinalizeSolutionStep( properties, geom, N, model_part.ProcessInfo )
