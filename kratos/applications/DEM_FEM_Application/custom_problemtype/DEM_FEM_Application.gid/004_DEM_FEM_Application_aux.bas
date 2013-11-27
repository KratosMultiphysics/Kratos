Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData



Begin Properties  0
End Properties



*loop materials
*#*if(strcmp(MatProp(Type),"Elements")==0)
*format "%4i"
Begin Properties  *MatNum
*format "%f"
DENSITY *MatProp(DENSITY,real)
*format "%f"
YOUNG_MODULUS *MatProp(YOUNG_MODULUS,real)
*format "%f"
POISSON_RATIO *MatProp(POISSON_RATIO,real)
*if(strcmp(GenData(GravityOption),"ON") == 0)
*format "%f%f%f"
BODY_FORCE [3] (*operation(GenData(Gravity_x,real)*MatProp(DENSITY,real)),*operation(GenData(Gravity_y,real)*MatProp(DENSITY,real)),*operation(GenData(Gravity_z,real)*MatProp(DENSITY,real)))
*endif
*format "%f"
DP_ALPHA1 *operation(6.0*sin(MatProp(Friction(Deg),real)/180.0*3.1415926)/(sqrt(3.0)*(3.0+sin(MatProp(Friction(Deg),real)/180.0*3.1415926))))
*format "%f"
DP_K *operation(6.0*MatProp(Cohesion,real)*cos(MatProp(Friction(Deg),real)/180.0*3.1415926)/(sqrt(3.0)*(3.0+sin(MatProp(Friction(Deg),real)/180.0*3.1415926))))
*format "%f"
PARTICLE_COHESION *MatProp(Cohesion,real)
*format "%f"
PARTICLE_FRICTION *MatProp(Friction(Deg),real)
*format "%f"
PARTICLE_TENSION *MatProp(Tension,real)
*format "%f"
THICKNESS 1.0
End Properties

*#//*endif
*end materials
*# Property blocks



Begin Nodes
*tcl(DEMFEM::WriteNonSphereAndCircleNodes *FileId)*\
End Nodes
*tcl(DEMFEM::ReleaseSphereAndCircleNodes)


*set elems(All)
*if( GenData(Domain_Dimension,int) == 3 )
Begin Conditions RigidFace3D3N
*loop elems *all
*if(strcmp(ElemsTypeName,"Triangle") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Conditions
*endif




*set elems(All)
*if( GenData(Domain_Dimension,int) == 3 )
Begin Conditions RigidFace3D4N
*loop elems *all
*if(strcmp(ElemsTypeName,"Quadrilateral") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Conditions
*endif



*set elems(All)
*if( GenData(Domain_Dimension,int) == 2 )
Begin Conditions RigidEdge3D2N
*loop elems *all
*if(strcmp(ElemsTypeName,"Linear") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Conditions
*endif

*#/////////////////////////////////////////////////////////FEM Element

*#//////////////////////////////////////////////////////////////

*set elems(All)
*if(strcmp(GenData(FEM_Option),"ON") == 0)
*if( GenData(Domain_Dimension,int) == 3 )
Begin Elements TotalLagrangian3D4N
*loop elems *all
*if(strcmp(ElemsTypeName,"Tetrahedra") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Elements
*endif
*endif

*#//////////////////////////////////////////////////////////////

*set elems(All)
*if(strcmp(GenData(FEM_Option),"ON") == 0)
*if( GenData(Domain_Dimension,int) == 3 )
Begin Elements TotalLagrangian3D8N
*loop elems *all
*if(strcmp(ElemsTypeName,"Hexahedra") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Elements
*endif
*endif

*#//////////////////////////////////////////////////////////////

*set elems(All)
*if(strcmp(GenData(FEM_Option),"ON") == 0)
*if( GenData(Domain_Dimension,int) == 2 )
Begin Elements TotalLagrangian2D3N
*loop elems *all
*if(strcmp(ElemsTypeName,"Triangle") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Elements
*endif
*endif


*#//////////////////////////////////////////////////////////////

*set elems(All)
*if(strcmp(GenData(FEM_Option),"ON") == 0)
*if( GenData(Domain_Dimension,int) == 2 )
Begin Elements TotalLagrangian2D4N
*loop elems *all
*if(strcmp(ElemsTypeName,"Quadrilateral") == 0)
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
    *operation(ElemsConec(*i))*\
*end

*endif
*end elems
End Elements
*endif
*endif

*#/////////////////////////////////////////////////////////////////////////
*set cond element_PointForce3D elem
*if(CondNumEntities > 0)
Begin Conditions PointForce3D
*loop elems *OnlyInCond
*format "%i%i%i"
*ElemsNum *ElemsMat *ElemsConec(1)
*end elems
End Conditions
*endif

*set cond element_PointForce2D elem
*if(CondNumEntities > 0)
Begin Conditions PointForce2D
*loop elems *OnlyInCond
*format "%i%i%i"
*ElemsNum *ElemsMat *ElemsConec(1)
*end elems
End Conditions
*endif

*#//////////////////////////////////////////////////////////////


*Set cond line_Face2D *elems
*if(CondNumEntities > 0)
Begin Conditions Face2D
*#// id prop_id  n1 n2  n3  ...
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

*#//////////////////////////////////////////////////////////////


*Set cond surface_Face3D3N *elems
*if(CondNumEntities > 0)
Begin Conditions Face3D3N
*#// id prop_id  n1 n2  n3  ...
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

*#//////////////////////////////////////////////////////////////


*Set cond surface_Face3D4N *elems
*if(CondNumEntities > 0)
Begin Conditions Face3D4N
*#// id prop_id  n1 n2  n3  ...
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
