*# Element and condition indices. We renumber them so each type is numbered from one.
*set var ielem=0
*set var icond=0
*# Define a condition index, which will be used to enforce that condition numbering begins from 1 
Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

*loop materials
*if(strcmp(MatProp(Type),"Elastic")==0)
*format "%i"
Begin Properties *MatNum
*format "%10.5e"
 DENSITY *MatProp(DENSITY,real)
*format "%10.5e"
 YOUNG_MODULUS *MatProp(YOUNG_MODULUS,real)
*format "%10.5e"
 POISSON_RATIO *MatProp(POISSON_RATIO,real)
*format "%10.5e"
 THICKNESS *MatProp(THICKNESS,real)
End Properties

*endif
*end materials
*loop materials
*if(strcmp(MatProp(Type),"Plastic")==0)
*format "%i"
Begin Properties *MatNum
*format "%10.5e"
 DENSITY *MatProp(DENSITY,real)
*format "%10.5e"
 YOUNG_MODULUS *MatProp(YOUNG_MODULUS,real)
*format "%10.5e"
 POISSON_RATIO *MatProp(POISSON_RATIO,real)
*format "%10.5e"
 YIELD_STRESS *MatProp(YIELD_STRESS,real)
*format "%10.5e"
 KINEMATIC_HARDENING_MODULUS *MatProp(KINEMATIC_HARDENING_MODULUS,real)
*format "%10.5e"
 HARDENING_EXPONENT *MatProp(HARDENING_EXPONENT,real)
*format "%10.5e"
 REFERENCE_HARDENING_MODULUS *MatProp(REFERENCE_HARDENING_MODULUS,real)
*format "%10.5e"
 INFINITY_HARDENING_MODULUS *MatProp(INFINITY_HARDENING,real)
*format "%10.5e"
 THICKNESS *MatProp(THICKNESS,real)
End Properties

*endif
*end materials
*# Property blocks

Begin Nodes
*#// id	  X	Y	Z
*loop nodes
*format "%i%10.5e%10.5e%10.5e"
*nodesnum	*NodesCoord(1)	*NodesCoord(2)	*NodesCoord(3)
*end nodes
End Nodes


*Set cond surface_SpatialLagrangianElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement2D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var ielem=operation(ielem+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif
*Set cond surface_SpatialLagrangianUPElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianUPElement2D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var ielem=operation(ielem+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif
*Set cond surface_AxisymSpatialLagrangianElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSpatialLagrangianElement2D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var ielem=operation(ielem+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif
*Set cond surface_AxisymSpatialLagrangianUPElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSpatialLagrangianUPElement2D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var ielem=operation(ielem+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i%i%i%i"
*ElemsNum *ElemsMat*\
*for(i=1;i<=j;i=i+1)*\
	*ElemsConec(*i)*\
*end

*end elems
End Elements

*endif
*# Element blocks

*# Line Condition Blocks

*Set cond line_LineLoadCondition2D2N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions LineLoadCondition2D2N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond line_LineLoadAxisymCondition2D2N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions LineLoadAxisymCondition2D2N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond line_WallTip2DCondition *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions WallTip2DCondition
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif

*# Point Condition Blocks

*Set cond point_PointLoad2DCondition *nodes
*if(CondNumEntities > 0)
Begin Conditions PointLoad2DCondition
*loop nodes *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i%i"
*icond 0 *NodesNum
*end nodes
End Conditions

*endif
*Set cond point_PointLoadAxisym2DCondition *nodes
*if(CondNumEntities > 0)
Begin Conditions PointLoadAxisym2DCondition
*loop nodes *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i%i"
*icond 0 *NodesNum
*end nodes
End Conditions

*endif

*# Variable Blocks

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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond surface_FACE_LOAD *nodes
*Add cond line_FACE_LOAD *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData FACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_X) *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData FACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Y) *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData FACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond point_FORCE *nodes
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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
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
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond surface_POSITIVE_FACE_PRESSURE *nodes
*Add cond line_POSITIVE_FACE_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData POSITIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(POSITIVE_FACE_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond surface_WALL_REFERENCE_POINT *nodes
*Add cond line_WALL_REFERENCE_POINT *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(WALL_REFERENCE_POINT_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData WALL_REFERENCE_POINT_X
*loop nodes *OnlyInCond
*if(cond(WALL_REFERENCE_POINT_X,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_X) *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop nodes *OnlyInCond
*if(cond(WALL_REFERENCE_POINT_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData WALL_REFERENCE_POINT_Y
*loop nodes *OnlyInCond
*if(cond(WALL_REFERENCE_POINT_Y,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Y) *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop nodes *OnlyInCond
*if(cond(WALL_REFERENCE_POINT_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData WALL_REFERENCE_POINT_Z
*loop nodes *OnlyInCond
*if(cond(WALL_REFERENCE_POINT_Z,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond surface_WALL_VELOCITY *nodes
*Add cond line_WALL_VELOCITY *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(WALL_VELOCITY_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData WALL_VELOCITY_X
*loop nodes *OnlyInCond
*if(cond(WALL_VELOCITY_X,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_X) *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Y value set
*set var Yset=0
*loop nodes *OnlyInCond
*if(cond(WALL_VELOCITY_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData WALL_VELOCITY_Y
*loop nodes *OnlyInCond
*if(cond(WALL_VELOCITY_Y,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Y) *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*#
*# Check if some node has its Z value set
*set var Zset=0
*loop nodes *OnlyInCond
*if(cond(WALL_VELOCITY_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData WALL_VELOCITY_Z
*loop nodes *OnlyInCond
*if(cond(WALL_VELOCITY_Z,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond surface_RIGID_WALL *nodes
*Add cond line_RIGID_WALL *nodes
*if(CondNumEntities > 0)
Begin NodalData RIGID_WALL
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(RIGID_WALL)
*end nodes
End NodalData

*endif
*Set cond surface_WALL_TIP_RADIUS *nodes
*Add cond line_WALL_TIP_RADIUS *nodes
*if(CondNumEntities > 0)
Begin NodalData WALL_TIP_RADIUS
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(WALL_TIP_RADIUS)
*end nodes
End NodalData

*endif

*#Tcl( WriteMeshGroups *FileId ) *\

*set var igroup=0
*set var icond=0
        
*loop groups
*if(strcmp(groupparentname,"Domains")==0)
*set var igroup=operation(igroup+1)
Begin Mesh *igroup
*#groupnum "*GroupFullName" ("*groupname" parent:*groupparentnum) *groupcolorrgb
*set group *GroupName *nodes

 Begin MeshNodes
*if(GroupNumEntities)
*#nodes: *GroupNumEntities
*loop nodes *onlyingroup
 *nodesnum
*end nodes
*end if
 End MeshNodes

 Begin MeshElements
*if(strcmp(groupname,"RigidWall"))
*set group *GroupName *elems 
*if(GroupNumEntities) 
*#elements: *GroupNumEntities 
*loop elems *onlyingroup 
 *elemsnum
*end elems
*end if
*end if
 End MeshElements

 Begin MeshConditions
*# Element Condition Blocks
*if(strcmp(groupname,"RigidWall"))
*# Line Condition Blocks
*set cond line_LineLoadCondition2D2N *OverFaceElements *CanRepeat
*add cond line_LineLoadAxisymCondition2D2N *OverFaceElements *CanRepeat
*add cond line_WallTip2DCondition *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
*loop elems *onlyincond *onlyingroup
*format "%i"
 *Tcl( getCondId *elemsnum *condelemface )
*end elems
*endif
*endif
*# Point Condition Blocks
*set group *GroupName *nodes
*set cond point_PointLoad2DCondition *nodes
*add cond point_PointLoadAxisym2DCondition *nodes
*if(CondNumEntities > 0)	    
*loop nodes *onlyincond *onlyingroup
*set var icond=operation(icond+1)
*format "%i"
 *icond
*end nodes
*endif
 End MeshConditions

End Mesh

*endif
*end groups              
 

*# Nodal Variable blocks

*# Elemental Variable blocks

*# Conditional Variable blocks

*# Note: About elements/conditions: it is important that point elements/conditions are added AFTER regular points/conditions to keep numeration of elemental/conditional data consistent.
*# This is why point elements/conditions get their own blocks.
*#
*Tcl(resetCondId) *\
*# Clear list of condition Ids
