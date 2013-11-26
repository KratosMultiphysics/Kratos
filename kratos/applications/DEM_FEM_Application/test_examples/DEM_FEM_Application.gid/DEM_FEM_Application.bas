Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

*# Property blocks
*loop materials
Begin Properties  1
End Properties
*end materials

*tcl(DEMFEM::SetSphereAndCircleNodes)
Begin Nodes
*tcl(DEMFEM::WriteSphereAndCircleNodes *FileId)*\
End Nodes


*set elems(sphere)
*Add elems(circle)
Begin Elements *GenData(DEM_Element_Type)
*loop elems 
*ElemsNum 1 *elemsconec(1)
*end elems
End Elements


Begin NodalData RADIUS
*set var iterator=1
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsradius
*endif
*end elems 
End NodalData


Begin NodalData PARTICLE_CONTINUUM
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(PARTICLE_CONTINUUM)
*endif
*end elems 
End NodalData

Begin NodalData PARTICLE_DENSITY
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(DENSITY)
*endif
*end elems 
End NodalData

Begin NodalData YOUNG_MODULUS
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(YOUNG_MODULUS)
*endif
*end elems 
End NodalData

Begin NodalData POISSON_RATIO
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(POISSON_RATIO)
*endif
*end elems 
End NodalData



Begin NodalData PARTICLE_FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *operation(tan(elemsmatprop(Friction(Deg),real)*3.141592653589793238462643383279502884197/180.0))
*endif
*end elems 
End NodalData


*if(strcmp(GenData(Rotation),"ON")==0)
*if(strcmp(GenData(Rota_Damp_Type),"RollingFric")==0)
Begin NodalData ROLLING_FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(LocalDampRatio)
*endif
*end elems 
End NodalData
*endif
*endif

*if(strcmp(GenData(Rotation),"ON")==0)
*if(strcmp(GenData(Rota_Damp_Type),"LocalDamp")==0)
Begin NodalData PARTICLE_ROTATION_DAMP_RATIO
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(LocalDampRatio)
*endif
*end elems 
End NodalData
*endif
*endif

Begin NodalData LN_OF_RESTITUTION_COEFF
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*if(elemsmatprop(Restitution_coef,real)==0.0)
*elemsconec(1) 0 1
*else
*elemsconec(1) 0 *operation(log(elemsmatprop(Restitution_coef,real)))
*endif
*endif
*end elems 
End NodalData


Begin NodalData PARTICLE_MATERIAL
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Colour)
*endif
*end elems 
End NodalData


Begin NodalData PARTICLE_COHESION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Cohesion)
*endif
*end elems 
End NodalData

Begin NodalData PARTICLE_TENSION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Tension)
*endif
*end elems 
End NodalData

Begin NodalData LOCAL_DAMP_RATIO
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(LocalDampRatio)
*endif
*end elems 
End NodalData
