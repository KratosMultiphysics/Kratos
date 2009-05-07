Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

*# Property blocks

Begin Nodes
*#// id	  X	Y	Z
*loop nodes
*nodesnum	*NodesCoord(1)	*NodesCoord(2)	*NodesCoord(3)
*end nodes
End Nodes

*Set cond surface_Fluid2D *elems
*if(CondNumEntities > 0)
Begin Elements Fluid2D
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
*Set cond volume_Fluid3D *elems
*if(CondNumEntities > 0)
Begin Elements Fluid3D
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
*Set cond surface_Fluid2DCoupled *elems
*if(CondNumEntities > 0)
Begin Elements Fluid2DCoupled
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
*Set cond volume_Fluid3DCoupled *elems
*if(CondNumEntities > 0)
Begin Elements Fluid3DCoupled
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
*Set cond surface_ASGS2D *elems
*if(CondNumEntities > 0)
Begin Elements ASGS2D
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
*Set cond volume_ASGS3D *elems
*if(CondNumEntities > 0)
Begin Elements ASGS3D
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
*Set cond surface_ASGSCompressible2D *elems
*if(CondNumEntities > 0)
Begin Elements ASGSCompressible2D
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
*# Element blocks

*Set cond line_Condition2D *elems
*if(CondNumEntities > 0)
Begin Conditions Condition2D
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
*Set cond surface_Condition3D *elems
*if(CondNumEntities > 0)
Begin Conditions Condition3D
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
*Set cond line_Face2DNeumann *elems
*if(CondNumEntities > 0)
Begin Conditions Face2DNeumann
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
*Set cond surface_Face3DNeumann *elems
*if(CondNumEntities > 0)
Begin Conditions Face3DNeumann
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
*Set cond volume_MESH_VELOCITY *nodes
*Add cond surface_MESH_VELOCITY *nodes
*Add cond line_MESH_VELOCITY *nodes
*Add cond point_MESH_VELOCITY *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(MESH_VELOCITY_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData MESH_VELOCITY_X
*loop nodes *OnlyInCond
*if(cond(MESH_VELOCITY_X,int)==1)
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
*if(cond(MESH_VELOCITY_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData MESH_VELOCITY_Y
*loop nodes *OnlyInCond
*if(cond(MESH_VELOCITY_Y,int)==1)
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
*if(cond(MESH_VELOCITY_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData MESH_VELOCITY_Z
*loop nodes *OnlyInCond
*if(cond(MESH_VELOCITY_Z,int)==1)
*format "%i%i%f"
*NodesNum *cond(Fix_Z) *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond volume_BODY_FORCE *nodes
*Add cond surface_BODY_FORCE *nodes
*Add cond line_BODY_FORCE *nodes
*Add cond point_BODY_FORCE *nodes
*if(CondNumEntities > 0)
*# Check if some node has its X value set
*set var Xset=0
*loop nodes *OnlyInCond
*if(cond(BODY_FORCE_X,int)==1)
*set var Xset=1
*endif
*end nodes
*if(Xset == 1)
Begin NodalData BODY_FORCE_X
*loop nodes *OnlyInCond
*if(cond(BODY_FORCE_X,int)==1)
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
*if(cond(BODY_FORCE_Y,int)==1)
*set var Yset=1
*endif
*end nodes
*if(Yset == 1)
Begin NodalData BODY_FORCE_Y
*loop nodes *OnlyInCond
*if(cond(BODY_FORCE_Y,int)==1)
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
*if(cond(BODY_FORCE_Z,int)==1)
*set var Zset=1
*endif
*end nodes
*if(Zset == 1)
Begin NodalData BODY_FORCE_Z
*loop nodes *OnlyInCond
*if(cond(BODY_FORCE_Z,int)==1)
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
*Set cond volume_AIR_PRESSURE *nodes
*Add cond surface_AIR_PRESSURE *nodes
*Add cond line_AIR_PRESSURE *nodes
*Add cond point_AIR_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData AIR_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(AIR_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond volume_WATER_PRESSURE *nodes
*Add cond surface_WATER_PRESSURE *nodes
*Add cond line_WATER_PRESSURE *nodes
*Add cond point_WATER_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData WATER_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(WATER_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond volume_EXTERNAL_PRESSURE *nodes
*Add cond surface_EXTERNAL_PRESSURE *nodes
*Add cond line_EXTERNAL_PRESSURE *nodes
*Add cond point_EXTERNAL_PRESSURE *nodes
*if(CondNumEntities > 0)
Begin NodalData EXTERNAL_PRESSURE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(EXTERNAL_PRESSURE)
*end nodes
End NodalData

*endif
*Set cond volume_VISCOSITY *nodes
*Add cond surface_VISCOSITY *nodes
*Add cond line_VISCOSITY *nodes
*Add cond point_VISCOSITY *nodes
*if(CondNumEntities > 0)
Begin NodalData VISCOSITY
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(VISCOSITY)
*end nodes
End NodalData

*endif
*Set cond volume_DENSITY *nodes
*Add cond surface_DENSITY *nodes
*Add cond line_DENSITY *nodes
*Add cond point_DENSITY *nodes
*if(CondNumEntities > 0)
Begin NodalData DENSITY
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(DENSITY)
*end nodes
End NodalData

*endif
*Set cond volume_DENSITY_AIR *nodes
*Add cond surface_DENSITY_AIR *nodes
*Add cond line_DENSITY_AIR *nodes
*Add cond point_DENSITY_AIR *nodes
*if(CondNumEntities > 0)
Begin NodalData DENSITY_AIR
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(DENSITY_AIR)
*end nodes
End NodalData

*endif
*Set cond volume_POROSITY *nodes
*Add cond surface_POROSITY *nodes
*Add cond line_POROSITY *nodes
*Add cond point_POROSITY *nodes
*if(CondNumEntities > 0)
Begin NodalData POROSITY
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(POROSITY)
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
*Set cond volume_FLAG_VARIABLE *nodes
*Add cond surface_FLAG_VARIABLE *nodes
*Add cond line_FLAG_VARIABLE *nodes
*Add cond point_FLAG_VARIABLE *nodes
*if(CondNumEntities > 0)
Begin NodalData FLAG_VARIABLE
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(FLAG_VARIABLE)
*end nodes
End NodalData

*endif
*Set cond volume_IS_BOUNDARY *nodes
*Add cond surface_IS_BOUNDARY *nodes
*Add cond line_IS_BOUNDARY *nodes
*Add cond point_IS_BOUNDARY *nodes
*if(CondNumEntities > 0)
Begin NodalData IS_BOUNDARY
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(IS_BOUNDARY)
*end nodes
End NodalData

*endif
*Set cond volume_IS_SLIP *nodes
*Add cond surface_IS_SLIP *nodes
*Add cond line_IS_SLIP *nodes
*Add cond point_IS_SLIP *nodes
*if(CondNumEntities > 0)
Begin NodalData IS_SLIP
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Fixed) *cond(IS_SLIP)
*end nodes
End NodalData

*endif
*# Nodal Variable blocks

*# Elemental Variable blocks

*Set cond surface_Slip_Face *elems
*Add cond line_Slip_Face *elems
*if(CondNumEntities > 0)
Begin ConditionalData Slip_Face
*loop elems *OnlyInCond
*format "%i%f"
*ElemsNum *cond(Slip_Face)
*end nodes
End ConditionalData

*endif
*# Conditional Variable blocks
