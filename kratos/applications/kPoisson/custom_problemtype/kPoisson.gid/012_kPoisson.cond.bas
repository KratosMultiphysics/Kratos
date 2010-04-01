// Fixing degrees of freedom in nodes

*Set cond Point_Fix_Value *nodes
*Add cond Line_Fix_Value *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum].Fix(DUMMY_UNKNOWN);
*end nodes

*Set cond Point_Dummy_Load *nodes
*loop nodes  *OnlyInCond
NODES[*NodesNum].Fix(DUMMY_POINT_SOURCE);
*end nodes