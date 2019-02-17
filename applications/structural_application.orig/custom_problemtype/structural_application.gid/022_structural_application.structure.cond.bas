*Set cond surface_Face3D *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
CONDITIONS[*ElemsNum] = Face3D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems
*Set cond line_Face2D *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
CONDITIONS[*ElemsNum] = Face2D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems


// Fixing degrees of freedom in nodes

*Set cond volume_PRESSURE_(Structure) *nodes
*Add cond surface_PRESSURE_(Structure) *nodes
*Add cond line_PRESSURE_(Structure) *nodes
*Add cond point_PRESSURE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum].Fix(PRESSURE);
*end nodes
*Set cond volume_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond surface_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond line_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond point_NEGATIVE_FACE_PRESSURE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum].Fix(NEGATIVE_FACE_PRESSURE);
*end nodes
*Set cond volume_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond surface_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond line_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*Add cond point_POSITIVE_FACE_PRESSURE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum].Fix(POSITIVE_FACE_PRESSURE);
*end nodes
*Set cond volume_IS_INTERFACE_(Structure) *nodes
*Add cond surface_IS_INTERFACE_(Structure) *nodes
*Add cond line_IS_INTERFACE_(Structure) *nodes
*Add cond point_IS_INTERFACE_(Structure) *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum].Fix(IS_INTERFACE);
*end nodes
*Set cond volume_DISPLACEMENT_(Structure) *nodes
*Add cond surface_DISPLACEMENT_(Structure) *nodes
*Add cond line_DISPLACEMENT_(Structure) *nodes
*Add cond point_DISPLACEMENT_(Structure) *nodes
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
*Set cond volume_VELOCITY_(Structure) *nodes
*Add cond surface_VELOCITY_(Structure) *nodes
*Add cond line_VELOCITY_(Structure) *nodes
*Add cond point_VELOCITY_(Structure) *nodes
*loop nodes  *OnlyInCond
*if(strcmp(cond(VELOCITY_X),"1")==0)
*format "%i%f"
NODES[*NodesNum].Fix(VELOCITY_X);
*endif
*if(strcmp(cond(VELOCITY_Y),"1")==0)
*format "%i%f"
NODES[*NodesNum].Fix(VELOCITY_Y);
*endif
*if(strcmp(cond(VELOCITY_Z),"1")==0)
*format "%i%f"
NODES[*NodesNum].Fix(VELOCITY_Z);
*endif
*end nodes
