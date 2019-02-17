*RealFormat "%10.5f"
*#
*#LOOP TO COUNT NODES
*#count nodes positive fluid side
*set var i=0
*Set Cond point_positive_Fluidboundary *nodes
*Add Cond line_positive_Fluidboundary *nodes
*Add Cond surface_positive_Fluidboundary *nodes
*loop nodes *OnlyInCond
*set var i=i+1   
*end
*#
*#count nodes negative fluid side
*set var j=0
*Set Cond point_negative_Fluidboundary *nodes
*Add Cond line_negative_Fluidboundary *nodes
*Add Cond surface_negative_Fluidboundary *nodes
*loop nodes *OnlyInCond
*set var j=j+1
*end
*#
*#count nodes structure
*set var k=0
*Set Cond point_Structureboundary *nodes
*Add Cond line_Structureboundary *nodes
*Add Cond surface_Structureboundary *nodes
*loop nodes *OnlyInCond
*set var k=k+1
*end
*#
*#
*#OUTPUT:
*#
*#
COMMUN BOUNDARY NODES

FLUID

FLUID_POSITIVE_SIDE
Number_Of_Nodes
*i
Nodes
*Set Cond point_positive_Fluidboundary *nodes
*Add Cond line_positive_Fluidboundary *nodes
*Add Cond surface_positive_Fluidboundary *nodes
*loop nodes *OnlyInCond
*nodesnum*\

*end

FLUID_NEGATIVE_SIDE
Number_Of_Nodes
*j
Nodes
*Set Cond point_negative_Fluidboundary *nodes
*Add Cond line_negative_Fluidboundary *nodes
*Add Cond surface_negative_Fluidboundary *nodes
*loop nodes *OnlyInCond
*nodesnum*\

*end

STRUCTURE_SIDE
Number_Of_Nodes
*k
Nodes
*Set Cond point_Structureboundary *nodes
*Add Cond line_Structureboundary *nodes
*Add Cond surface_Structureboundary *nodes
*loop nodes *OnlyInCond
*nodesnum*\

*end
