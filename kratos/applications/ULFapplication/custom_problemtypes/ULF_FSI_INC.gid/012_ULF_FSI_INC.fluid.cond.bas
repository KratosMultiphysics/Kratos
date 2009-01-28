*Set cond surface_Condition3D *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
CONDITIONS[*ElemsNum] = Condition3D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems
*Set cond line_Condition2D *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
CONDITIONS[*ElemsNum] = Condition2D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems


// Fixing degrees of freedom in nodes

*Set cond volume_BULK_MODULUS_(Fluid) *nodes
*Add cond surface_BULK_MODULUS_(Fluid) *nodes
*Add cond line_BULK_MODULUS_(Fluid) *nodes
*Add cond point_BULK_MODULUS_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_VISCOSITY_(Fluid) *nodes
*Add cond surface_VISCOSITY_(Fluid) *nodes
*Add cond line_VISCOSITY_(Fluid) *nodes
*Add cond point_VISCOSITY_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_DENSITY_(Fluid) *nodes
*Add cond surface_DENSITY_(Fluid) *nodes
*Add cond line_DENSITY_(Fluid) *nodes
*Add cond point_DENSITY_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_IS_STRUCTURE_(Fluid) *nodes
*Add cond surface_IS_STRUCTURE_(Fluid) *nodes
*Add cond line_IS_STRUCTURE_(Fluid) *nodes
*Add cond point_IS_STRUCTURE_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_IS_INTERFACE_(Fluid) *nodes
*Add cond surface_IS_INTERFACE_(Fluid) *nodes
*Add cond line_IS_INTERFACE_(Fluid) *nodes
*Add cond point_IS_INTERFACE_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_IS_BOUNDARY_(Fluid) *nodes
*Add cond surface_IS_BOUNDARY_(Fluid) *nodes
*Add cond line_IS_BOUNDARY_(Fluid) *nodes
*Add cond point_IS_BOUNDARY_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*Add cond surface_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*Add cond line_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*Add cond point_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
*end nodes
*Set cond volume_DISPLACEMENT_(Fluid) *nodes
*Add cond surface_DISPLACEMENT_(Fluid) *nodes
*Add cond line_DISPLACEMENT_(Fluid) *nodes
*Add cond point_DISPLACEMENT_(Fluid) *nodes
*loop nodes  *OnlyInCond
*if(strcmp(cond(DISPLACEMENT_X),"1")==0)
*format "%i%f"
NODES[*NodesNum].Fix(DISPLACEMENT_X);
*endif
*if(strcmp(cond(DISPLACEMENT_Y),"1")==0)
*format "%i%f"
NODES[*NodesNum].Fix(DISPLACEMENT_Y);
*endif
*if(strcmp(cond(DISPLACEMENT_Z),"1")==0)
*format "%i%f"
NODES[*NodesNum].Fix(DISPLACEMENT_Z);
*endif
*end nodes
*Set cond volume_BODY_FORCE_(Fluid) *nodes
*Add cond surface_BODY_FORCE_(Fluid) *nodes
*Add cond line_BODY_FORCE_(Fluid) *nodes
*Add cond point_BODY_FORCE_(Fluid) *nodes
*loop nodes  *OnlyInCond
*if(strcmp(cond(BODY_FORCE_X),"1")==0)
*format "%i%f"
*endif
*if(strcmp(cond(BODY_FORCE_Y),"1")==0)
*format "%i%f"
*endif
*if(strcmp(cond(BODY_FORCE_Z),"1")==0)
*format "%i%f"
*endif
*end nodes
