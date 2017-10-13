*RealFormat "%10.5f"
NODES = NodesList([
*set var i=1
*Set Cond point_Structure *nodes
*Add Cond line_Structure *nodes
*Add Cond surface_Structure *nodes
*Add Cond volume_Structure *nodes
*loop nodes *OnlyInCond
*if(i>1)
,
*endif
[*nodesnum,*NodesCoord(1),*NodesCoord(2),*NodesCoord(3)]*\
*set var i=i+1
*end

])
