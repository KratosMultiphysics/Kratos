// Fixing degrees of freedom in nodes

*Set cond volume_PRESSURE_(Structure) *nodes
*Add cond surface_PRESSURE_(Structure) *nodes
*Add cond line_PRESSURE_(Structure) *nodes
*Add cond point_PRESSURE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](PRESSURE,0) = *cond(Value);
*end nodes
*Set cond volume_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond surface_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond line_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond point_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](NEGATIVE_FACE_PRESSURE,0) = *cond(Value);
*end nodes
*Set cond volume_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond surface_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond line_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond point_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](POSITIVE_FACE_PRESSURE,0) = *cond(Value);
*end nodes
*Set cond volume_IS_INTERFACE_(Structure) *nodes
*Add cond surface_IS_INTERFACE_(Structure) *nodes
*Add cond line_IS_INTERFACE_(Structure) *nodes
*Add cond point_IS_INTERFACE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](IS_INTERFACE,0) = *cond(Value);
*end nodes
*Set cond volume_DISPLACEMENT_(Structure) *nodes
*Add cond surface_DISPLACEMENT_(Structure) *nodes
*Add cond line_DISPLACEMENT_(Structure) *nodes
*Add cond point_DISPLACEMENT_(Structure) *nodes
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
*Set cond volume_VELOCITY_(Structure) *nodes
*Add cond surface_VELOCITY_(Structure) *nodes
*Add cond line_VELOCITY_(Structure) *nodes
*Add cond point_VELOCITY_(Structure) *nodes
*loop nodes  *OnlyInCond
*if(strcmp(cond(VELOCITY_X),"1")==0)
*format "%i%f"
NODES[*NodesNum](VELOCITY_X,0) = *cond(Value_X);
*endif
*if(strcmp(cond(VELOCITY_Y),"1")==0)
*format "%i%f"
NODES[*NodesNum](VELOCITY_Y,0) = *cond(Value_Y);
*endif
*if(strcmp(cond(VELOCITY_Z),"1")==0)
*format "%i%f"
NODES[*NodesNum](VELOCITY_Z,0) = *cond(Value_Z);
*endif
*end nodes
