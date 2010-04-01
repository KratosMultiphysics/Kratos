// Fixing degrees of freedom in nodes

*Set cond Point_Fix_Value *nodes
*Add cond Line_Fix_Value *nodes
*loop nodes  *OnlyInCond
*format "%i%f"
NODES[*NodesNum](DUMMY_UNKNOWN,0) = *cond(Value);
*end nodes
*Set cond Point_Dummy_Load *nodes
*loop nodes  *OnlyInCond
NODES[*NodesNum](DUMMY_POINT_SOURCE,0) = *cond(Q);
*end nodes
