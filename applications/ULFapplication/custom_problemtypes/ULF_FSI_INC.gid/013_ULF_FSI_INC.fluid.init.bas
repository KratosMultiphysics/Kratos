// Fixing degrees of freedom in nodes

*Set cond volume_BULK_MODULUS_(Fluid) *nodes
*Add cond surface_BULK_MODULUS_(Fluid) *nodes
*Add cond line_BULK_MODULUS_(Fluid) *nodes
*Add cond point_BULK_MODULUS_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](BULK_MODULUS,0) = *cond(Value);
*end nodes
*Set cond volume_VISCOSITY_(Fluid) *nodes
*Add cond surface_VISCOSITY_(Fluid) *nodes
*Add cond line_VISCOSITY_(Fluid) *nodes
*Add cond point_VISCOSITY_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](VISCOSITY,0) = *cond(Value);
*end nodes
*Set cond volume_DENSITY_(Fluid) *nodes
*Add cond surface_DENSITY_(Fluid) *nodes
*Add cond line_DENSITY_(Fluid) *nodes
*Add cond point_DENSITY_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](DENSITY,0) = *cond(Value);
*end nodes
*Set cond volume_IS_STRUCTURE_(Fluid) *nodes
*Add cond surface_IS_STRUCTURE_(Fluid) *nodes
*Add cond line_IS_STRUCTURE_(Fluid) *nodes
*Add cond point_IS_STRUCTURE_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](IS_STRUCTURE,0) = *cond(Value);
*end nodes
*Set cond volume_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*Add cond surface_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*Add cond line_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*Add cond point_IS_LAGRANGIAN_INLET_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](IS_LAGRANGIAN_INLET,0) = *cond(Value);
*end nodes
*Set cond volume_IS_INTERFACE_(Fluid) *nodes
*Add cond surface_IS_INTERFACE_(Fluid) *nodes
*Add cond line_IS_INTERFACE_(Fluid) *nodes
*Add cond point_IS_INTERFACE_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](IS_INTERFACE,0) = *cond(Value);
*end nodes
*Set cond volume_IS_BOUNDARY_(Fluid) *nodes
*Add cond surface_IS_BOUNDARY_(Fluid) *nodes
*Add cond line_IS_BOUNDARY_(Fluid) *nodes
*Add cond point_IS_BOUNDARY_(Fluid) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](IS_BOUNDARY,0) = *cond(Value);
*end nodes
*Set cond volume_DISPLACEMENT_(Fluid) *nodes
*Add cond surface_DISPLACEMENT_(Fluid) *nodes
*Add cond line_DISPLACEMENT_(Fluid) *nodes
*Add cond point_DISPLACEMENT_(Fluid) *nodes
*loop nodes  *OnlyInCond
*if(strcmp(cond(DISPLACEMENT_X),"1")==0)
*format "%i%f"
NODES[*NodesNum](DISPLACEMENT_X,0) = *cond(Value_X);
*endif
*if(strcmp(cond(DISPLACEMENT_Y),"1")==0)
*format "%i%f"
NODES[*NodesNum](DISPLACEMENT_Y,0) = *cond(Value_Y);
*endif
*if(strcmp(cond(DISPLACEMENT_Z),"1")==0)
*format "%i%f"
NODES[*NodesNum](DISPLACEMENT_Z,0) = *cond(Value_Z);
*endif
*end nodes
*Set cond volume_BODY_FORCE_(Fluid) *nodes
*Add cond surface_BODY_FORCE_(Fluid) *nodes
*Add cond line_BODY_FORCE_(Fluid) *nodes
*Add cond point_BODY_FORCE_(Fluid) *nodes
*loop nodes  *OnlyInCond
*if(strcmp(cond(BODY_FORCE_X),"1")==0)
*format "%i%f"
NODES[*NodesNum](BODY_FORCE_X,0) = *cond(Value_X);
*endif
*if(strcmp(cond(BODY_FORCE_Y),"1")==0)
*format "%i%f"
NODES[*NodesNum](BODY_FORCE_Y,0) = *cond(Value_Y);
*endif
*if(strcmp(cond(BODY_FORCE_Z),"1")==0)
*format "%i%f"
NODES[*NodesNum](BODY_FORCE_Z,0) = *cond(Value_Z);
*endif
*end nodes
