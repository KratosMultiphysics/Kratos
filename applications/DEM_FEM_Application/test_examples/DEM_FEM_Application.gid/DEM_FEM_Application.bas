Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

*loop materials
*#*if(strcmp(MatProp(Type),"Elements")==0)
Begin Properties  *MatNum
*format "%f"
DENSITY *MatProp(DENSITY,real)
*format "%f"
YOUNG_MODULUS *MatProp(YOUNG_MODULUS,real)
*format "%f"
POISSON_RATIO *MatProp(POISSON_RATIO,real)
*format "%f%f%f"
GRAVITY [3] (*MatProp(Gravity_X,real),*MatProp(Gravity_Y,real),*MatProp(Gravity_Z,real))
*format "%f%f%f"
BODY_FORCE [3] (*operation(MatProp(Gravity_X,real)*MatProp(DENSITY,real)),*operation(MatProp(Gravity_Y,real)*MatProp(DENSITY,real)),*operation(MatProp(Gravity_Z,real)*MatProp(DENSITY,real)))
*format "%f"
DP_ALPHA1 *operation(6.0*sin(MatProp(Friction(Deg),real)/180.0*3.1415926)/(sqrt(3.0)*(3.0+sin(MatProp(Friction(Deg),real)/180.0*3.1415926))))
*format "%f"
DP_K *operation(6.0*MatProp(Cohesion,real)*cos(MatProp(Friction(Deg),real)/180.0*3.1415926)/(sqrt(3.0)*(3.0+sin(MatProp(Friction(Deg),real)/180.0*3.1415926))))
*format "%f"
THICKNESS *MatProp(THICKNESS,real)
*format "%f"
PARTICLE_COHESION *MatProp(Cohesion,real)
*format "%f"
PARTICLE_FRICTION *MatProp(Friction(Deg),real)
*format "%f"
PARTICLE_TENSION *MatProp(Tension,real)

End Properties

*#//*endif
*end materials
*# Property blocks

Begin Nodes
*#// id	  X	Y	Z
*loop nodes
*nodesnum	*NodesCoord(1)	*NodesCoord(2)	*NodesCoord(3)
*end nodes
End Nodes


Begin NodalData RADIUS
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
  *ElemsConec(1) 0  *elemsradius
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

Begin NodalData PARTICLE_CONTINUUM
*loop elems *all
*if(strcmp(ElemsTypeName,"Sphere")==0 || strcmp(ElemsTypeName,"Circle")==0)
*elemsconec(1) 0 *elemsmatprop(PARTICLE_CONTINUUM)
*endif
*end elems 
End NodalData

*#//////////////////////////////////////////////////////////////

*Set cond surface_TotalLagrangian2D3N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangian2D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif

*#//////////////////////////////////////////////////////////////

*Set cond surface_TotalLagrangian2D4N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangian2D4N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif

*#//////////////////////////////////////////////////////////////

*Set cond volume_TotalLagrangian3D4N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangian3D4N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif

*#//////////////////////////////////////////////////////////////

*Set cond volume_TotalLagrangian3D8N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangian3D8N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif

*#//////////////////////////////////////////////////////////////

*Set cond surface_DEM_FEM_Particle2D *elems
*if(CondNumEntities > 0)
Begin Elements DEM_FEM_Particle2D
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif

*#//////////////////////////////////////////////////////////////

*Set cond volume_DEM_FEM_Particle3D *elems
*if(CondNumEntities > 0)
Begin Elements DEM_FEM_Particle3D
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif

*#//////////////////////////////////////////////////////////////

*# Element blocks

*set cond element_PointForce3D elem
*if(CondNumEntities > 0)
Begin Conditions PointForce3D
*loop elems *OnlyInCond
*format "%i%i%i"
*ElemsNum *ElemsMat	*ElemsConec(1)
*end elems
End Conditions

*endif
*set cond element_PointForce2D elem
*if(CondNumEntities > 0)
Begin Conditions PointForce2D
*loop elems *OnlyInCond
*format "%i%i%i"
*ElemsNum *ElemsMat	*ElemsConec(1)
*end elems
End Conditions

*endif
*Set cond line_Face2D *elems
*if(CondNumEntities > 0)
Begin Conditions Face2D
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Conditions

*endif
*Set cond surface_Face3D3N *elems
*if(CondNumEntities > 0)
Begin Conditions Face3D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Conditions

*endif
*Set cond surface_Face3D4N *elems
*if(CondNumEntities > 0)
Begin Conditions Face3D4N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Conditions

*endif
*# Condition Blocks

*Set cond volume_DISPLACEMENT *nodes
*Add cond surface_DISPLACEMENT *nodes
*Add cond line_DISPLACEMENT *nodes
*Add cond point_DISPLACEMENT *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(DISPLACEMENT_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData DISPLACEMENT_X
*loop nodes *OnlyInCond
*if(cond(DISPLACEMENT_X,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_X) *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop nodes *OnlyInCond
*if(cond(DISPLACEMENT_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData DISPLACEMENT_Y
*loop nodes *OnlyInCond
*if(cond(DISPLACEMENT_Y,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Y) *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop nodes *OnlyInCond
*if(cond(DISPLACEMENT_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData DISPLACEMENT_Z
*loop nodes *OnlyInCond
*if(cond(DISPLACEMENT_Z,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond volume_VELOCITY *nodes
*Add cond surface_VELOCITY *nodes
*Add cond line_VELOCITY *nodes
*Add cond point_VELOCITY *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(VELOCITY_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData VELOCITY_X
*loop nodes *OnlyInCond
*if(cond(VELOCITY_X,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_X) *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop nodes *OnlyInCond
*if(cond(VELOCITY_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData VELOCITY_Y
*loop nodes *OnlyInCond
*if(cond(VELOCITY_Y,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Y) *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop nodes *OnlyInCond
*if(cond(VELOCITY_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData VELOCITY_Z
*loop nodes *OnlyInCond
*if(cond(VELOCITY_Z,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond volume_FORCE *nodes
*Add cond surface_FORCE *nodes
*Add cond line_FORCE *nodes
*Add cond point_FORCE *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(FORCE_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData FORCE_X
*loop nodes *OnlyInCond
*if(cond(FORCE_X,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_X) *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop nodes *OnlyInCond
*if(cond(FORCE_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData FORCE_Y
*loop nodes *OnlyInCond
*if(cond(FORCE_Y,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Y) *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop nodes *OnlyInCond
*if(cond(FORCE_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData FORCE_Z
*loop nodes *OnlyInCond
*if(cond(FORCE_Z,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond volume_PRESSURE *nodes
*Add cond surface_PRESSURE *nodes
*Add cond line_PRESSURE *nodes
*Add cond point_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(PRESSURE)
*end nodes
End NodalData

*endif
*Set cond volume_NEGATIVE_FACE_PRESSURE *nodes
*Add cond surface_NEGATIVE_FACE_PRESSURE *nodes
*Add cond line_NEGATIVE_FACE_PRESSURE *nodes
*Add cond point_NEGATIVE_FACE_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData NEGATIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(NEGATIVE_FACE_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond volume_POSITIVE_FACE_PRESSURE *nodes
*Add cond surface_POSITIVE_FACE_PRESSURE *nodes
*Add cond line_POSITIVE_FACE_PRESSURE *nodes
*Add cond point_POSITIVE_FACE_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData POSITIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(POSITIVE_FACE_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond volume_IS_INTERFACE *nodes
*Add cond surface_IS_INTERFACE *nodes
*Add cond line_IS_INTERFACE *nodes
*Add cond point_IS_INTERFACE *nodes
*if(CondNumEntities > 0)
Begin NodalData IS_INTERFACE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(IS_INTERFACE)
*end nodes
End NodalData

*endif
*# Nodal Variable blocks

*# Elemental Variable blocks

*# Conditional Variable blocks
