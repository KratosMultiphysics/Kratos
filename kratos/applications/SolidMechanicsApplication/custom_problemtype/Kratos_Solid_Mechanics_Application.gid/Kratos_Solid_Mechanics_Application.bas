*# Element and condition indices. We renumber them so each type is numbered from one.
*set var ielem=0
*set var icond=0
*# Define a condition index, which will be used to enforce that condition numbering begins from 1 
Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties


*loop materials
*if(strcmp(MatProp(Type),"Elastic")==0)
*format "%i"
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME *MatProp(CONSTITUTIVE_LAW_NAME)
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
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME *MatProp(CONSTITUTIVE_LAW_NAME)
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
INFINITY_HARDENING_MODULUS *MatProp(INFINITY_HARDENING_MODULUS,real)
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

*#// SMALL DISPLACEMENT
*Set cond surface_SmallDisplacementElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement2D3N
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
*Set cond surface_SmallDisplacementElement2D4N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement2D4N
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
*Set cond surface_SmallDisplacementElement2D6N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement2D6N
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
*Set cond surface_SmallDisplacementElement2D8N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement2D8N
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
*Set cond volume_SmallDisplacementElement3D4N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D4N
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
*Set cond volume_SmallDisplacementElement3D6N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D6N
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
*Set cond volume_SmallDisplacementElement3D8N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D8N
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
*Set cond volume_SmallDisplacementElement3D10N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D10N
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
*Set cond volume_SmallDisplacementElement3D15N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D15N
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
*Set cond volume_SmallDisplacementElement3D20N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D20N
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
*Set cond volume_SmallDisplacementElement3D27N *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementElement3D27N
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
*#// TOTAL LAGRANGIAN
*Set cond surface_TotalLagrangianElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement2D3N
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
*Set cond surface_TotalLagrangianElement2D4N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement2D4N
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
*Set cond surface_TotalLagrangianElement2D6N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement2D6N
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
*Set cond surface_TotalLagrangianElement2D8N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement2D8N
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
*Set cond volume_TotalLagrangianElement3D4N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D4N
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
*Set cond volume_TotalLagrangianElement3D6N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D6N
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
*Set cond volume_TotalLagrangianElement3D8N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D8N
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
*Set cond volume_TotalLagrangianElement3D10N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D10N
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
*Set cond volume_TotalLagrangianElement3D15N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D15N
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
*Set cond volume_TotalLagrangianElement3D20N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D20N
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
*Set cond volume_TotalLagrangianElement3D27N *elems
*if(CondNumEntities > 0)
Begin Elements TotalLagrangianElement3D27N
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
*#// UPDATED LAGRANGIAN
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
*Set cond surface_UpdatedLagrangianElement2D4N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement2D4N
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
*Set cond surface_UpdatedLagrangianElement2D6N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement2D6N
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
*Set cond surface_UpdatedLagrangianElement2D8N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement2D8N
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
*Set cond volume_UpdatedLagrangianElement3D4N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D4N
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
*Set cond volume_UpdatedLagrangianElement3D6N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D6N
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
*Set cond volume_UpdatedLagrangianElement3D8N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D8N
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
*Set cond volume_UpdatedLagrangianElement3D10N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D10N
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
*Set cond volume_UpdatedLagrangianElement3D15N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D15N
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
*Set cond volume_UpdatedLagrangianElement3D20N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D20N
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
*Set cond volume_UpdatedLagrangianElement3D27N *elems
*if(CondNumEntities > 0)
Begin Elements UpdatedLagrangianElement3D27N
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
*#// SPATIAL LAGRANGIAN
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
*Set cond surface_SpatialLagrangianElement2D4N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement2D4N
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
*Set cond surface_SpatialLagrangianElement2D6N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement2D6N
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
*Set cond surface_SpatialLagrangianElement2D8N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement2D8N
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
*Set cond volume_SpatialLagrangianElement3D4N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D4N
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
*Set cond volume_SpatialLagrangianElement3D6N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D6N
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
*Set cond volume_SpatialLagrangianElement3D8N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D8N
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
*Set cond volume_SpatialLagrangianElement3D10N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D10N
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
*Set cond volume_SpatialLagrangianElement3D15N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D15N
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
*Set cond volume_SpatialLagrangianElement3D20N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D20N
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
*Set cond volume_SpatialLagrangianElement3D27N *elems
*if(CondNumEntities > 0)
Begin Elements SpatialLagrangianElement3D27N
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
*#// U-P ELEMENTS
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
*#// AXISYM SMALL DISPLACEMENT
*Set cond surface_AxisymSmallDisplacementElement2D3N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSmallDisplacementElement2D3N
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
*Set cond surface_AxisymSmallDisplacementElement2D4N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSmallDisplacementElement2D4N
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
*Set cond surface_AxisymSmallDisplacementElement2D6N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSmallDisplacementElement2D6N
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
*Set cond surface_AxisymSmallDisplacementElement2D8N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSmallDisplacementElement2D8N
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
*#// AXISYM SPATIAL LAGRANGIAN
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
*Set cond surface_AxisymSpatialLagrangianElement2D4N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSpatialLagrangianElement2D4N
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
*Set cond surface_AxisymSpatialLagrangianElement2D6N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSpatialLagrangianElement2D6N
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
*Set cond surface_AxisymSpatialLagrangianElement2D8N *elems
*if(CondNumEntities > 0)
Begin Elements AxisymSpatialLagrangianElement2D8N
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

*# Point Condition Blocks

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
*Set cond line_LineLoadCondition2D3N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions LineLoadCondition2D3N
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
*Set cond line_LineLoadAxisymCondition2D3N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions LineLoadAxisymCondition2D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond surface_SurfaceLoadCondition3D3N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions SurfaceLoad3DCondition3D3N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond surface_SurfaceLoadCondition3D4N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions SurfaceLoad3DCondition3D4N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond surface_SurfaceLoadCondition3D6N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions SurfaceLoad3DCondition3D6N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond surface_SurfaceLoadCondition3D8N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions SurfaceLoad3DCondition3D8N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif
*Set cond surface_SurfaceLoadCondition3D9N *OverFaceElements *CanRepeat
*if(CondNumEntities > 0)
Begin Conditions SurfaceLoad3DCondition3D9N
*#// id prop_id	 n1	n2	n3	...
*loop elems *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i"
*Tcl( setCondId *ElemsNum *CondElemFace ) *ElemsMat *GlobalNodes*\

*end elems
End Conditions

*endif

*# Condition Blocks

*Set cond point_PointLoad3DCondition *nodes
*if(CondNumEntities > 0)
Begin Conditions PointLoad3DCondition
*loop nodes *OnlyInCond
*set var icond=operation(icond+1)
*format "%i%i%i"
*icond 0 *NodesNum
*end nodes
End Conditions

*endif
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
*# Point Condition Blocks

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

*#Tcl( WriteMeshGroups *FileId ) *\

*# Nodal Variable blocks

*# Elemental Variable blocks

*# Conditional Variable blocks

*# Note: About elements/conditions: it is important that point elements/conditions are added AFTER regular points/conditions to keep numeration of elemental/conditional data consistent.
*# This is why point elements/conditions get their own blocks.
*#
*Tcl(resetCondId) *\
*# Clear list of condition Ids
