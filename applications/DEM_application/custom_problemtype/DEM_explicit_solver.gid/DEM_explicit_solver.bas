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


Begin Nodes
*#// id	  X	Y	Z
*loop nodes
*nodesnum	*NodesCoord(1)	*NodesCoord(2)	*NodesCoord(3)
*end nodes
End Nodes
*# Element blocks
*# Condition Blocks

*set elems(sphere)
Begin Elements *GenData(DEM_Element_Type)
*loop elems 
*ElemsNum 1 *elemsconec(1)
*end elems
End Elements

*Set cond volume_VELOCITY *elems
*Add cond surface_VELOCITY *elems
*Add cond line_VELOCITY *elems
*Add cond test_PARTICLES *elems
*#Add cond point_VELOCITY *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop elems *OnlyInCond
*if(cond(VELOCITY_X,int)==1)
*set var Xset=1
*break
*endif
*end elems
*if(Xset == 1)
Begin NodalData VELOCITY_X
*loop elems *OnlyInCond
*if(cond(VELOCITY_X,int)==1)
*elemsconec(1) *cond(X_fixed) *cond(X_Value)
*endif
*end elems
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop elems *OnlyInCond
*if(cond(VELOCITY_Y,int)==1)
*set var Yset=1
*break
*endif
*end elems
*if(Yset == 1)
Begin NodalData VELOCITY_Y
*loop elems *OnlyInCond
*if(cond(VELOCITY_Y,int)==1)
*elemsconec(1) *cond(Y_fixed) *cond(Y_Value)
*endif
*end elems

End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop elems *OnlyInCond
*if(cond(VELOCITY_Z,int)==1)
*set var Zset=1
*break
*endif
*end elems
*if(Zset == 1)
Begin NodalData VELOCITY_Z
*loop elems *OnlyInCond
*if(cond(VELOCITY_Z,int)==1)
*elemsconec(1) *cond(Z_fixed) *cond(Z_Value)
*endif
*end elems

End NodalData
*#
*endif
*# Check if some node has its X value set
*set var X_ang_vel_set=0
*loop elems *OnlyInCond
*if(cond(ANGULAR_VELOCITY_X,int)==1)
*set var X_ang_vel_set=1
*break
*endif
*end elems
*if(X_ang_vel_set == 1)
Begin NodalData ANGULAR_VELOCITY_X
*loop elems *OnlyInCond
*if(cond(ANGULAR_VELOCITY_X,int)==1)
*elemsconec(1) *cond(X_ang_vel_fixed) *cond(X_ang_vel_Value)
*endif
*end elems
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Y_ang_vel_set=0
*loop elems *OnlyInCond
*if(cond(ANGULAR_VELOCITY_Y,int)==1)
*set var Y_ang_vel_set=1
*break
*endif
*end elems
*if(Y_ang_vel_set == 1)
Begin NodalData ANGULAR_VELOCITY_Y
*loop elems *OnlyInCond
*if(cond(ANGULAR_VELOCITY_Y,int)==1)
*elemsconec(1) *cond(Y_ang_vel_fixed) *cond(Y_ang_vel_Value)
*endif
*end elems

End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Z_ang_vel_set=0
*loop elems *OnlyInCond
*if(cond(ANGULAR_VELOCITY_Z,int)==1)
*set var Z_ang_vel_set=1
*break
*endif
*end elems
*if(Z_ang_vel_set == 1)
Begin NodalData ANGULAR_VELOCITY_Z
*loop elems *OnlyInCond
*if(cond(ANGULAR_VELOCITY_Z,int)==1)
*elemsconec(1) *cond(Z_ang_vel_fixed) *cond(Z_ang_vel_Value)
*endif
*end elems
End NodalData

*endif

*endif

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

*Set cond volume_GROUP_ID *elems
*Add cond test_PARTICLES *elems

Begin NodalData GROUP_ID
*loop elems *OnlyInCond
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *cond(GROUP_ID)
*endif
*end elems 
End NodalData

Begin NodalData PARTICLE_CONTINUUM
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

Begin NodalData PARTICLE_FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Friction(Deg))
*endif
*end elems 
End NodalData

*if(strcmp(GenData(Rota_Damp_Id),"RollingFric")==0)
Begin NodalData ROLLING_FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Rolling_Friction)
*endif
*end elems 
End NodalData
*endif

*if(strcmp(GenData(Rota_Damp_Id),"LocalDamp")==0)
Begin NodalData PARTICLE_ROTATION_DAMP_RATIO
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(RotaDampRatio)
*endif
*end elems 
End NodalData
*endif


Begin NodalData RESTITUTION_COEFF
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Restitution_coef)
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

*Set cond SET_SKIN *elems
Begin ElementalData PREDEFINED_SKIN
*loop elems *OnlyInCond
*elemsnum 1.0
*end elems
End ElementalData

*#Begin NodalData PARTICLE_STATIC_FRICTION_COEF
*#*loop elems *all
*#*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*#*elemsconec(1) 0 *elemsmatprop(Static_friction_coef)
*#*endif
*#*end elems 
*#End NodalData

*#Begin NodalData PARTICLE_DYNAMIC_FRICTION_COEF
*#*loop elems *all
*#*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*#*elemsconec(1) 0 *elemsmatprop(Dynamic_friction_coef)
*#*endif
*#*end elems 
*#End NodalData
