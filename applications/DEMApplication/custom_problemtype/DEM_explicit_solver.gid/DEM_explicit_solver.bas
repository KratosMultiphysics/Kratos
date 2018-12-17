*# Define a condition index, which will be used to enforce that condition numbering begins from 1
*set var condbase=-1
*intformat "%i"        
*realformat "%10.5e"

Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

*# Property blocks
*loop materials
Begin Properties  1
End Properties
*end materials

*tcl(DEM::SetSphereAndCircleNodes)
Begin Nodes
*tcl(DEM::WriteSphereAndCircleNodes *FileId)*\
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
*tcl(DemElems_Addtolist *iterator *elemsconec(1))*\
*set var iterator=iterator+1
*endif
*end elems 
End NodalData

*Set cond group_CONDITIONS *groups
Begin NodalData VELOCITY_X
*loop groups *OnlyInCond
*if(cond(VELOCITY_X,int))
*set group *GroupName *elems
*loop elems *onlyingroup
*elemsconec(1) *cond(X_fixed) *cond(X_Value)
*end elems
*endif
*end groups
End NodalData

Begin NodalData VELOCITY_Y
*loop groups *OnlyInCond
*if(cond(VELOCITY_Y,int))
*set group *GroupName *elems
*loop elems *onlyingroup
*elemsconec(1) *cond(Y_fixed) *cond(Y_Value)
*end elems
*endif
*end groups
End NodalData

Begin NodalData VELOCITY_Z
*loop groups *OnlyInCond
*if(cond(VELOCITY_Z,int))
*set group *GroupName *elems
*loop elems *onlyingroup
*elemsconec(1) *cond(Z_fixed) *cond(Z_Value)
*end elems
*endif
*end groups
End NodalData

Begin NodalData ANGULAR_VELOCITY_X
*loop groups *OnlyInCond
*if(cond(ANGULAR_VELOCITY_X,int))
*set group *GroupName *elems
*loop elems *onlyingroup
*elemsconec(1) *cond(X_ang_vel_fixed) *cond(X_ang_vel_Value)
*end elems
*endif
*end groups
End NodalData

Begin NodalData ANGULAR_VELOCITY_Y
*loop groups *OnlyInCond
*if(cond(ANGULAR_VELOCITY_Y,int))
*set group *GroupName *elems
*loop elems *onlyingroup
*elemsconec(1) *cond(Y_ang_vel_fixed) *cond(Y_ang_vel_Value)
*end elems
*endif
*end groups
End NodalData

Begin NodalData ANGULAR_VELOCITY_Z
*loop groups *OnlyInCond
*if(cond(ANGULAR_VELOCITY_Z,int))
*set group *GroupName *elems
*loop elems *onlyingroup
*elemsconec(1) *cond(Z_ang_vel_fixed) *cond(Z_ang_vel_Value)
*end elems
*endif
*end groups
End NodalData

Begin NodalData GROUP_ID
*loop groups *OnlyInCond
*set group *GroupName *elems
*loop elems *onlyingroup
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *cond(GROUP_ID)
*endif
*end elems
*end groups
End NodalData

Begin NodalData COHESIVE_GROUP
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Continuum_group)
*endif
*end elems 
End NodalData

Begin NodalData PARTICLE_DENSITY
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Particle_density)
*endif
*end elems 
End NodalData

Begin NodalData YOUNG_MODULUS
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Young_modulus)
*endif
*end elems 
End NodalData

Begin NodalData POISSON_RATIO
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Poisson_ratio)
*endif
*end elems 
End NodalData

Begin NodalData FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*#elemsconec(1) 0 *operation(tan(elemsmatprop(Friction(Deg),real)*3.141592653589793238462643383279502884197/180.0))
*elemsconec(1) 0 *elemsmatprop(Dynamic_Friction)

*endif
*end elems 
End NodalData

*if(strcmp(GenData(Rota_Damp_Type),"RollingFric")==0)
Begin NodalData ROLLING_FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Rolling_Friction)
*endif
*end elems 
End NodalData
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

*Set cond volume_SET_SKIN_MANUALLY *elems
*Add cond surface_SET_SKIN_MANUALLY *elems
Begin NodalData PREDEFINED_SKIN
*loop elems *OnlyInCond
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *cond(PREDEFINED_SKIN)
*endif
*end elems 
End NodalData
