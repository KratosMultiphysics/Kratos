*RealFormat "%10.5f"
NODES = NodesList([
*set var i=1
*loop nodes
*if(i>1)
,
*endif
[*nodesnum,*NodesCoord(1),*NodesCoord(2),*NodesCoord(3)]*\
*set var i=i+1
*end

])

