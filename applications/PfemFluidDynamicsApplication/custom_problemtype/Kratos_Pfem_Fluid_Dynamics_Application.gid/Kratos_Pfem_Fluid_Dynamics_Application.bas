*# Element and condition indices. We renumber them so each type is numbered from one.
*set var ielem=0
*set var icond=0

*# ModelPart block

Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

*# Property blocks

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
*if(strcmp(MatProp(Type),"Damage")==0)
*format "%i"
Begin Properties *MatNum
*format "%10.5e"
 DENSITY *MatProp(DENSITY,real)
*format "%10.5e"
 YOUNG_MODULUS *MatProp(YOUNG_MODULUS,real)
*format "%10.5e"
 POISSON_RATIO *MatProp(POISSON_RATIO,real)
*if(strcmp(MatProp(HARDENING_MODEL),"SIMO_JU")==0)
*format "%10.5e"
 STRENGTH_RATIO *MatProp(STRENGTH_RATIO,real)
*format "%10.5e"
 FRACTURE_ENERGY *MatProp(FRACTURE_ENERGY,real)
*format "%10.5e"
 DAMAGE_THRESHOLD *MatProp(DAMAGE_THRESHOLD,real)
*endif
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
*if(strcmp(MatProp(HARDENING_MODEL),"SIMO")==0)
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
*elseif(strcmp(MatProp(HARDENING_MODEL),"BASIC_TRESCA")==0)
*format "%10.5e"
 YIELD_STRESS *MatProp(YIELD_STRESS,real)
*format "%10.5e"
 DENSITY_WATER *MatProp(WATER_DENSITY,real)
*format "%10.5e"
 WATER_BULK_MODULUS *MatProp(WATER_BULK_MODULUS,real)
*format "%10.5e"
 PERMEABILITY *MatProp(PERMEABILITY,real)
*format "%10.5e"
 STABILIZATION_FACTOR *MatProp(STABILIZATION_FACTOR,real)
*format "%10.5e"
 CONTACT_ADHESION *MatProp(CONTACT_ADHESION,real)
*format "%10.5e"
 CONTACT_FRICTION_ANGLE *MatProp(CONTACT_FRICTION_ANGLE,real)
*format "%10.5e"
 K0 *MatProp(K0,real)
*elseif(strcmp(MatProp(HARDENING_MODEL),"MOHR_C")==0)
*format "%10.5e"
 INTERNAL_FRICTION_ANGLE *MatProp(INTERNAL_FRICTION_ANGLE,real)
*format "%10.5e"
 INTERNAL_DILATANCY_ANGLE *MatProp(INTERNAL_DILATANCY_ANGLE,real)
*format "%10.5e"
 COHESION *MatProp(COHESION,real)
*format "%10.5e"
 DENSITY_WATER *MatProp(WATER_DENSITY,real)
*format "%10.5e"
 WATER_BULK_MODULUS *MatProp(WATER_BULK_MODULUS,real)
*format "%10.5e"
 PERMEABILITY *MatProp(PERMEABILITY,real)
*format "%10.5e"
 STABILIZATION_FACTOR *MatProp(STABILIZATION_FACTOR,real)
*format "%10.5e"
 CONTACT_ADHESION *MatProp(CONTACT_ADHESION,real)
*format "%10.5e"
 CONTACT_FRICTION_ANGLE *MatProp(CONTACT_FRICTION_ANGLE,real)
*format "%10.5e"
 K0 *MatProp(K0,real)
*endif
*format "%10.5e"
 THICKNESS *MatProp(THICKNESS,real)
End Properties

*endif
*end materials
*loop materials
*if(strcmp(MatProp(Type),"CriticalStatePlasticity")==0)
*format "%i"
Begin Properties *MatNum
*format "%10.5e"
 DENSITY *MatProp(DENSITY,real)
*format "%10.5e"
 YOUNG_MODULUS *MatProp(YOUNG_MODULUS,real)
*format "%10.5e"
 SWELLING_SLOPE *MatProp(SWELLING_SLOPE,real)
*format "%10.5e"
 INITIAL_SHEAR_MODULUS *MatProp(INITIAL_SHEAR_MODULUS,real)
*format "%10.5e"
 ALPHA_SHEAR *MatProp(ALPHA_SHEAR,real)
*format "%10.5e"
 PRE_CONSOLIDATION_STRESS *MatProp(PRE_CONSOLIDATION_STRESS,real)
*format "%10.5e"
 OVER_CONSOLIDATION_RATIO *MatProp(OVER_CONSOLIDATION_RATIO,real)
*format "%10.5e"
 NORMAL_COMPRESSION_SLOPE *MatProp(NORMAL_COMPRESSION_SLOPE,real)
*format "%10.5e"
 CRITICAL_STATE_LINE *MatProp(CRITICAL_STATE_LINE,real)
*format "%10.5e"
 INTERNAL_FRICTION_ANGLE *MatProp(INTERNAL_FRICTION_ANGLE,real)
*format "%10.5e"
 DENSITY_WATER *MatProp(WATER_DENSITY,real)
*format "%10.5e"
 WATER_BULK_MODULUS *MatProp(WATER_BULK_MODULUS,real)
*format "%10.5e"
 PERMEABILITY *MatProp(PERMEABILITY,real)
*format "%10.5e"
 STABILIZATION_FACTOR *MatProp(STABILIZATION_FACTOR,real)
*format "%10.5e"
 CONTACT_ADHESION *MatProp(CONTACT_ADHESION,real)
*format "%10.5e"
 CONTACT_FRICTION_ANGLE *MatProp(CONTACT_FRICTION_ANGLE,real)
*format "%10.5e"
 K0 *MatProp(K0,real)
*format "%10.5e"
 THICKNESS *MatProp(THICKNESS,real)
End Properties

*endif
*end materials


*# Mesh0 block

*# Nodes block
Begin Nodes
*#// id	  X	Y	Z
*loop nodes
*format "%i%10.5e%10.5e%10.5e"
*NodesNum	*NodesCoord(1)	*NodesCoord(2)	*NodesCoord(3)
*end nodes
End Nodes


*# Element blocks

*Set cond surface_UpdatedLagrangianElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement2D3N
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
*Set cond surface_UpdatedLagrangianUwPElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianUwPElement2D3N
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
*Set cond surface_UpdatedLagrangianUwPStabElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianUwPStabElement2D3N
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
*Set cond surface_AxisymUpdatedLagrangianUwPStabElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymUpdatedLagrangianUwPStabElement2D3N
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
*Set cond surface_UpdatedLagrangianUPElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianUPElement2D3N
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
*Set cond surface_UpdatedLagrangianUPStabElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianUPStabElement2D3N
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
*Set cond surface_AxisymUpdatedLagrangianUPStabElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymUpdatedLagrangianUPStabElement2D3N
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
*Set cond surface_AxisymUpdatedLagrangianElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymUpdatedLagrangianElement2D3N
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
*Set cond surface_AxisymUpdatedLagrangianUPElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymUpdatedLagrangianUPElement2D3N
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
*# Condition Blocks

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

*# Point Condition Blocks

*Set cond point_PointLoadCondition2D1N *nodes
*if(CondNumEntities > 0)
Begin Conditions PointLoadCondition2D1N
*loop nodes *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i%i"
*Tcl( setCondId *NodesNum 0 ) 0 *NodesNum
*end nodes
End Conditions

*endif
*Set cond point_AxisymPointLoadCondition2D1N *nodes
*if(CondNumEntities > 0)
Begin Conditions AxisymPointLoadCondition2D1N
*loop nodes *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i%i"
*Tcl( setCondId *NodesNum 0 ) 0 *NodesNum
*end nodes
End Conditions

*endif
*# Group Condition Blocks

*# Set the start number for each condition:
*set var RigidWallsstart  = icond

*Set cond group_RigidWalls *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
*if(strcmp(cond(Contact_Condition),"3D")==0)
Begin Conditions WallCondition3D3N
*else
Begin Conditions WallCondition2D2N
*endif
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *elems
*loop elems *onlyingroup
*set var icond=operation(icond+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i"
*icond *ElemsMat *\
*for(i=1;i<=j;i=i+1)*\
 *ElemsConec(*i)*\
*end

*end elems
End Conditions

*end groups


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
*Set cond volume_VOLUME_ACCELERATION *nodes
*Add cond surface_VOLUME_ACCELERATION *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(VOLUME_ACCELERATION_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData VOLUME_ACCELERATION_X
*loop nodes *OnlyInCond
*if(cond(VOLUME_ACCELERATION_X,int)==1)
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
*if(cond(VOLUME_ACCELERATION_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData VOLUME_ACCELERATION_Y
*loop nodes *OnlyInCond
*if(cond(VOLUME_ACCELERATION_Y,int)==1)
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
*if(cond(VOLUME_ACCELERATION_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData VOLUME_ACCELERATION_Z
*loop nodes *OnlyInCond
*if(cond(VOLUME_ACCELERATION_Z,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond line_LINE_LOAD *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData LINE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_X,int)==1)
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
*if(cond(LINE_LOAD_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData LINE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_Y,int)==1)
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
*if(cond(LINE_LOAD_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData LINE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_Z,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond surface_SURFACE_LOAD *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData SURFACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_X,int)==1)
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
*if(cond(SURFACE_LOAD_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData SURFACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_Y,int)==1)
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
*if(cond(SURFACE_LOAD_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData SURFACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_Z,int)==1)
*format "%i%i%10.5e"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond point_POINT_LOAD *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(FORCE_LOAD_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData FORCE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(FORCE_LOAD_X,int)==1)
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
*if(cond(FORCE_LOAD_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData FORCE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(FORCE_LOAD_Y,int)==1)
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
*if(cond(FORCE_LOAD_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData FORCE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(FORCE_LOAD_Z,int)==1)
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
*Set cond surface_WATER_PRESSURE *nodes
*Add cond line_WATER_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData WATER_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(WATER_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond surface_WALL_TIP *nodes
*Add cond line_WALL_TIP *nodes
*if(CondNumEntities > 0)
Begin NodalData WALL_TIP_RADIUS
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(WALL_TIP_RADIUS)
*end nodes
End NodalData

Begin NodalData WALL_REFERENCE_POINT_X
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(X_Value)
*end nodes
End NodalData

Begin NodalData WALL_REFERENCE_POINT_Y
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(Y_Value)
*end nodes
End NodalData

Begin NodalData WALL_REFERENCE_POINT_Z
*loop nodes *OnlyInCond
*format "%i%i%10.5e"
*NodesNum *cond(Fixed) *cond(Z_Value)
*end nodes
End NodalData

*endif
*Set cond group_RigidWalls *groups
*if(CondNumEntities > 0)
Begin NodalData RIGID_WALL
*loop groups *OnlyInCond
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
*format "%i%i%i"
*NodesNum 0 *cond(Group_ID,int)
*end nodes
*endif
*end groups
End NodalData

*endif


*# Mesh Blocks


*# Mesh Blocks for Domain Remeshing and Contact

*Set cond group_DeformableBodies *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Mesh *cond(Group_ID)

 Begin MeshNodes
*set group *GroupName *nodes
*if(GroupNumEntities > 0)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End MeshNodes

 Begin MeshElements
*set group *GroupName *elems 
*if(GroupNumEntities > 0) 
*loop elems *onlyingroup 
*format "%i"
 *ElemsNum
*end elems
*endif
 End MeshElements
      
 Begin MeshConditions
*set group *GroupName *elems 
*if(GroupNumEntities > 0) 
*# Line Condition Blocks
*set cond line_LineLoadCondition2D2N *OverFaceElements *CanRepeat
*add cond line_LineLoadAxisymCondition2D2N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
*loop elems *onlyincond *onlyingroup
*format "%i"
 *Tcl( getCondId *ElemsNum *CondElemFace )
*end elems
*endif
*endif
*# Point Condition Blocks
*set group *GroupName *nodes
*if(GroupNumEntities > 0)
*set cond point_PointLoadCondition2D1N *nodes
*add cond point_AxisymPointLoadCondition2D1N *nodes
*if(CondNumEntities > 0)	    
*loop nodes *onlyincond *onlyingroup
*set var point_conditions_num=operation(point_conditions_num+1)
*format "%i"
 *Tcl( getCondId *NodesNum 0 )
*end nodes
*endif
*endif
*Set cond group_DeformableBodies *groups
 End MeshConditions

End Mesh

*end groups    
*endif


*# Set the start number for each element and condition:
*set var RigidWallsNum  = RigidWallsStart

*Set cond group_RigidWalls *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Mesh *cond(Group_ID)

 Begin MeshNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End MeshNodes

 Begin MeshElements
 End MeshElements
      
 Begin MeshConditions
*set group *GroupName *elems 
*if(GroupNumEntities) 
*loop elems *onlyingroup 
*set var RigidWallsNum=operation(RigidWallsNum+1)
*format "%i"
 *RigidWallsNum
*end elems
*endif
 End MeshConditions

End Mesh

*end groups        
*endif



*# Note: About elements/conditions: it is important that point elements/conditions are added AFTER regular points/conditions to keep numeration of elemental/conditional data consistent.
*# This is why point elements/conditions get their own blocks.
*#
*Tcl(resetCondId) *\
*# Clear list of condition Ids
