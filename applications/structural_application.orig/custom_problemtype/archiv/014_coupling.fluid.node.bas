*RealFormat "%10.5f"
NODES = NodesList([
*set var i=1
*Set Cond point_Fluid *nodes
*Add Cond line_Fluid *nodes
*Add Cond surface_Fluid *nodes
*Add Cond volume_Fluid *nodes
*loop nodes *OnlyInCond
*if(i>1)
,
*endif
[*nodesnum,*NodesCoord(1),*NodesCoord(2),*NodesCoord(3)]*\
*set var i=i+1
*end

])

