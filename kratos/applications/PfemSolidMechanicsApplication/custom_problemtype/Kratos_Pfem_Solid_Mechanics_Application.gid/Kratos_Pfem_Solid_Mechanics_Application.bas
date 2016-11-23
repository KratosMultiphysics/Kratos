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

*set cond surface_UpdatedLagrangianElement2D3N *elems
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
*set cond surface_UpdatedLagrangianUwPElement2D3N *elems
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
*set cond surface_UpdatedLagrangianUwPStabElement2D3N *elems
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
*set cond surface_AxisymUpdatedLagrangianUwPStabElement2D3N *elems
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
*set cond surface_UpdatedLagrangianUPElement2D3N *elems
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
*set cond surface_UpdatedLagrangianUPStabElement2D3N *elems
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
*set cond surface_AxisymUpdatedLagrangianUPStabElement2D3N *elems
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
*set cond surface_AxisymUpdatedLagrangianElement2D3N *elems
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
*set cond surface_AxisymUpdatedLagrangianUPElement2D3N *elems
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
*set cond group_RigidBodies *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
*if(strcmp(cond(Parametric_Wall),"False")==0)
*if(strcmp(cond(Body_Surface),"False")==0)
*if(GenData(DOMAIN_SIZE,INT) == 3)
Begin Elements Element3D4N
*else
Begin Elements Element2D3N
*endif
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *elems
*loop elems *onlyingroup
*set var ielem=operation(ielem+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i"
*ielem *ElemsMat *\
*for(i=1;i<=j;i=i+1)*\
 *ElemsConec(*i)*\
*end

*end elems
End Elements

*endif
*endif
*end groups
*endif

*# Condition Blocks
*# start number for each condition type:
*set var RigidWallsstart = icond

*set cond group_RigidBodies *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
*if(strcmp(cond(Body_Surface),"True")==0 || strcmp(cond(Parametric_Wall),"True")==0)
*if(GenData(DOMAIN_SIZE,INT) == 3)
Begin Elements Element3D3N
*else
Begin Elements Element2D2N
*endif
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *elems
*loop elems *onlyingroup
*set var ielem=operation(ielem+1)
*set var i=0
*set var j=ElemsNnode
*format "%i%i%i%i%i"
*ielem *ElemsMat *\
*for(i=1;i<=j;i=i+1)*\
 *ElemsConec(*i)*\
*end

*end elems
End Elements

*endif
*end groups
*endif


*set cond group_POINT_LOAD *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Conditions *cond(ConditionType)
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *nodes
*loop nodes *onlyingroup
*format "%i%i%i"
 *Tcl( setCondId *NodesNum 0 ) 0 *NodesNum
*end nodes
End Conditions

*end groups
*endif

*set cond group_LINE_LOAD *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Conditions *cond(ConditionType)
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *faces
*loop faces *onlyingroup
 *Tcl(Print2DFaceElement *FaceElemsNum *FaceIndex)
*end faces
End Conditions

*end groups
*endif

*set cond group_SURFACE_LOAD *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Conditions *cond(ConditionType)
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *faces
*loop faces *onlyingroup
 *Tcl(Print3DFaceElement *FaceElemsNum *FaceIndex)
*end faces
End Conditions

*end groups
*endif

*set cond group_LINE_PRESSURE *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Conditions *cond(ConditionType)
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *faces
*loop faces *onlyingroup
 *Tcl(Print2DFaceElement *FaceElemsNum *FaceIndex)
*end faces
End Conditions

*end groups
*endif

*set cond group_SURFACE_PRESSURE *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Conditions *cond(ConditionType)
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *faces
*loop faces *onlyingroup
 *Tcl(Print3DFaceElement *FaceElemsNum *FaceIndex)
*end faces
End Conditions

*end groups
*endif

*set cond group_POINT_MOMENT *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin Conditions *cond(ConditionType)
*#// id prop_id	 n1	n2	n3	...
*set group *GroupName *nodes
*loop nodes *onlyingroup
*format "%i%i%i"
 *Tcl( setCondId *NodesNum 0 ) 0 *NodesNum
*end nodes
End Conditions

*end groups
*endif


*# SubModelPart Blocks
*set cond group_DeformableBodies *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities > 0)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
*set group *GroupName *elems 
*if(GroupNumEntities > 0) 
*loop elems *onlyingroup 
*format "%i"
 *ElemsNum
*end elems
*endif
 End SubModelPartElements
      
 Begin SubModelPartConditions
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_LINEAR_MOVEMENT *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_ANGULAR_MOVEMENT *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_POINT_LOAD *groups
*if(CondNumEntities)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
*format "%i"
 *Tcl( getCondId *NodesNum 0 )
*end nodes
*endif
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_LINE_LOAD *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
*set group *GroupName *faces
*if(GroupNumEntities)
*loop faces *onlyingroup
 *Tcl( getCondId *FaceElemsNum *FaceIndex )
*end faces
*endif
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_LINE_PRESSURE *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
*set group *GroupName *faces
*if(GroupNumEntities)
*loop faces *onlyingroup
 *Tcl( getCondId *FaceElemsNum *FaceIndex )
*end faces
*endif
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_SURFACE_LOAD *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
*set group *GroupName *faces
*if(GroupNumEntities)
*loop faces *onlyingroup
 *Tcl( getCondId *FaceElemsNum *FaceIndex )
*end faces
*endif
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_SURFACE_PRESSURE *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
*set group *GroupName *faces
*if(GroupNumEntities)
*loop faces *onlyingroup
 *Tcl( getCondId *FaceElemsNum *FaceIndex )
*end faces
*endif
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_POINT_MOMENT *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
*format "%i"
 *Tcl( getCondId *NodesNum 0 )
*end nodes
*endif
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif
*set cond group_VOLUME_ACCELERATION *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
 End SubModelPartElements
      
 Begin SubModelPartConditions
 End SubModelPartConditions

End SubModelPart
*end groups    
*endif

*set var RigidWallsElemNum = 0
*set var RigidWallsCondNum = 0
*set cond group_RigidBodies *groups
*if(CondNumEntities > 0)
*loop groups *OnlyInCond
Begin SubModelPart *GroupName // *GroupNum

 Begin SubModelPartNodes
*set group *GroupName *nodes
*if(GroupNumEntities)
*loop nodes *onlyingroup
 *NodesNum
*end nodes
*endif
 End SubModelPartNodes

 Begin SubModelPartElements
*set group *GroupName *elems 
*if(GroupNumEntities > 0) 
*loop elems *onlyingroup 
*format "%i"
 *ElemsNum
*end elems
*endif
 End SubModelPartElements
      
 Begin SubModelPartConditions
 End SubModelPartConditions
End SubModelPart
*end groups        
*endif

*# Note: About elements/conditions: it is important that point elements/conditions are added AFTER regular points/conditions to keep numeration of elemental/conditional data consistent.
*# This is why point elements/conditions get their own blocks.
*#
*Tcl(resetCondId) *\
*# Clear list of condition Ids
