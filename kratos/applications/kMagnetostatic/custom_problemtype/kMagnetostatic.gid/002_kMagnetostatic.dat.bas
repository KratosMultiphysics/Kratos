Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

*loop materials
Begin Properties *MatNum 
*format "%i%f%f%f"  
MAGNETIC_PERMEABILITY [3] (*MatProp(X_Magnetic_Permeability), *MatProp(Y_Magnetic_Permeability), 0.0)
COERCIVITY [3] (*MatProp(X_Coercivity), *MatProp(Y_Coercivity), 0.0)
End Properties
*end materials

Begin Nodes
*RealFormat "%10.5f"
*loop nodes
*nodesnum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
End Nodes

Begin Elements Magnetostatic2D
*loop elems
*ElemsNum *ElemsMat *ElemsConec
*end elems
End Elements

Begin ElementalData MAGNETOSTATIC_SURFACE_CURRENT
*Set cond Surface_Current *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Js)
*end elems
End ElementalData

Begin NodalData MAGNETOSTATIC_POTENTIAL
*Set cond Surface_Potential *nodes
*Add cond Line_Potential *nodes
*Add cond Point_Potential *nodes
*loop nodes  *OnlyInCond
*format "%i%i%f"
*NodesNum 1 *cond(Potential)
*end nodes
End NodalData

Begin Conditions PointCurrent2D
*Set cond Point_Current *nodes
*loop nodes *OnlyInCond
*format "%i%i%f"
*NodesNum 1 *NodesNum
*end nodes
End Conditions

Begin Conditions Mfield2D
*Set cond Contour_Magnetic_Intensity *elems
*Add cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *ElemsMat *GlobalNodes
*end elems
End Conditions

Begin ConditionalData MAGNETOSTATIC_POINT_CURRENT
*Set cond Point_Current *nodes
*loop nodes  *OnlyInCond
*format "%i%i%f"
*NodesNum *cond(I)
*end nodes
End ConditionalData

Begin ConditionalData MAGNETIC_FIELD_INTENSITY
*Set cond Contour_Magnetic_Intensity *elems
*loop elems *OnlyInCond
*ElemsNum [3] (*cond(Hn) ,*cond(Hn) ,*cond(Hn)) //vector
*end elems
End ConditionalData

Begin ConditionalData INFINIT_COEFFICIENT
*Set cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Exponent)
*end elems
End ConditionalData




