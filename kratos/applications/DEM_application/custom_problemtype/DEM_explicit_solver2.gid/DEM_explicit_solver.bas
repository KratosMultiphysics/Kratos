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

Begin Elements SphericParticle
*loop elems 
*format "%i%i%i%i%i%i%i%i"
*ElemsNum 1 *elemsconec(1)
*end elems
End Elements

*Set cond volume_VELOCITY *elems
*Add cond surface_VELOCITY *elems
*Add cond line_VELOCITY *elems
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
*elemsconec(1) 1 *cond(X_Value)
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
*elemsconec(1) 1 *cond(Y_Value)
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
*elemsconec(1) 1 *cond(Z_Value)
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

Begin NodalData PARTICLE_FRICTION
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Friction(Deg))
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

Begin NodalData PARTICLE_LOCAL_DAMP_RATIO
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(LocalDampRatio)
*endif
*end elems 
End NodalData


Begin NodalData PARTICLE_ZETA
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(Zeta_value)
*endif
*end elems 
End NodalData


Begin NodalData PARTICLE_COEF_RESTITUTION
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