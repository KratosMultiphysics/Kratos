Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

*loop materials
Begin Properties *MatNum 
*format "%i%f%f%f"  
ELECTRICAL_PERMITTIVITY [3] (*MatProp(Electrical_Permittivity), *MatProp(Electrical_Permittivity), *MatProp(Electrical_Permittivity))
End Properties
*end materials

Begin Nodes
*RealFormat "%10.5f"
*loop nodes
*nodesnum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
End Nodes

Begin Elements Electrostatic3D
*loop elems
*ElemsNum *ElemsMat *ElemsConec
*end elems
End Elements

Begin ElementalData ELECTROSTATIC_VOLUME_CHARGE
*Set cond Volume_Electrical_Charge *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Qv)
*end elems
End ElementalData

Begin NodalData ELECTROSTATIC_POTENTIAL
*Set cond Volume_Voltage *nodes
*Add cond Surface_Voltage *nodes
*Add cond Line_Voltage *nodes
*Add cond Point_Voltage *nodes
*loop nodes  *OnlyInCond
*format "%i%i%f"
*NodesNum 1 *cond(Voltage)
*end nodes
End NodalData

Begin Conditions PointCharge3D
*Set cond Point_Electrical_Charge *nodes
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum 1 *NodesNum
*end nodes
End Conditions

Begin Conditions Efield3D
*Set cond Contour_Electric_Displacement *elems
*Add cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *ElemsMat *GlobalNodes
*end elems
End Conditions

Begin ConditionalData ELECTROSTATIC_POINT_CHARGE
*Set cond Point_Electrical_Charge *nodes
*loop nodes  *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(Q)
*end nodes
End ConditionalData

Begin ConditionalData ELECTRIC_DISPLACEMENT_FIELD
*Set cond Contour_Electric_Displacement *elems
*loop elems *OnlyInCond
*ElemsNum [3] (*cond(Dn) ,*cond(Dn) ,*cond(Dn)) //vector
*end elems
End ConditionalData

Begin ConditionalData INFINIT_COEFFICIENT
*Set cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Exponent)
*end elems
End ConditionalData




