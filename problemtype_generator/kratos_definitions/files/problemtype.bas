*# Element and condition indices. We renumber them so each type is numbered from one.
*set var ielem=0
*set var icond=0
*# Define a condition index, which will be used to enforce that condition numbering begins from 1 
Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData

Begin Properties 0
End Properties

*# Property blocks

Begin Nodes
*#// id	  X	Y	Z
*loop nodes
*format "%i%10.5e%10.5e%10.5e"
*nodesnum	*NodesCoord(1)	*NodesCoord(2)	*NodesCoord(3)
*end nodes
End Nodes

*# Element blocks

*# Point Element Blocks

*# Condition Blocks

*# Point Condition Blocks

*# Nodal Variable blocks

*# Elemental Variable blocks

*# Conditional Variable blocks

*# Note: About elements/conditions: it is important that point elements/conditions are added AFTER regular points/conditions to keep numeration of elemental/conditional data consistent.
*# This is why point elements/conditions get their own blocks.
*#
*Tcl(resetCondId) *\
*# Clear list of condition Ids
