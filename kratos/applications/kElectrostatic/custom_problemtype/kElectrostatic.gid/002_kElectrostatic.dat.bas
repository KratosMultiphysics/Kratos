Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

*loop materials
Begin Properties *MatNum 
*format "%i%f%f%f"  
ELECTRICAL_PERMITTIVITY [3] (*MatProp(Electrical_Permittivity), *MatProp(Electrical_Permittivity), 0.0)
End Properties
*end materials

Begin Nodes
*RealFormat "%10.5f"
*loop nodes
*nodesnum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
End Nodes

Begin Elements Electrostatic2D
*loop elems
*ElemsNum *ElemsMat *ElemsConec
*end elems
End Elements

Begin ElementalData ELECTROSTATIC_SURFACE_CHARGE
*Set cond Surface_Electrical_Charge *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Qs)
*end elems
End ElementalData

Begin NodalData ELECTROSTATIC_POTENTIAL
*Set cond Surface_Voltage *nodes
*Add cond Line_Voltage *nodes
*Add cond Point_Voltage *nodes
*loop nodes  *OnlyInCond
*format "%i%i%f"
*NodesNum 1 *cond(Voltage)
*end nodes
End NodalData

Begin NodalData ELECTROSTATIC_POINT_CHARGE
*Set cond Point_Electrical_Charge *nodes
*format "%i%i%f"
*loop nodes  *OnlyInCond
*NodesNum 0 *cond(Q)
*end nodes
End NodalData

Begin Conditions PointCharge2D
*Set cond Point_Electrical_Charge *nodes
*format "%i%i%f"
*loop nodes *OnlyInCond
1 1 *NodesNum
*end nodes
End Conditions

Begin ElementalData ELECTRIC_DISPLACEMENT_FIELD
*Set cond Contour_Electric_Displacement *elems
*loop elems *OnlyInCond
*ElemsNum [3] (*cond(Dn) ,*cond(Dn) ,*cond(Dn)) //vector
*end elems
End ElementalData

Begin ElementalData INFINIT_COEFFICIENT
*Set cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Exponent)
*end elems
End ElementalData

Begin Conditions Efield2D
*Set cond Contour_Electric_Displacement *elems
*Add cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *ElemsMat *LocalNodes
*end elems
End Conditions


