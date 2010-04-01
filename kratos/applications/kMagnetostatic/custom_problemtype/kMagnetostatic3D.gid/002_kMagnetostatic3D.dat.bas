Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

*loop materials
Begin Properties *MatNum 
*format "%i%f%f%f"  
MAGNETIC_PERMEABILITY [3] (*MatProp(X_Magnetic_Permeability), *MatProp(Y_Magnetic_Permeability), *MatProp(Z_Magnetic_Permeability))
COERCIVITY [3] (*MatProp(X_Coercivity), *MatProp(Y_Coercivity), *MatProp(Z_Coercivity))
End Properties
*end materials

Begin Nodes
*RealFormat "%10.5f"
*loop nodes
*nodesnum *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
End Nodes

Begin Elements Magnetostatic3D
*loop elems
*ElemsNum *ElemsMat *ElemsConec
*end elems
End Elements

Begin ElementalData MAGNETOSTATIC_VOLUME_CURRENT
*Set cond Volume_Current *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Jx) *cond(Jy) *cond(Jz)
*end elems
End ElementalData

Begin ElementalData MAGNETOSTATIC_SURFACE_CURRENT
*Set cond Surface_Current *elems
*loop elems *OnlyInCond
*ElemsNum *cond(Js)
*end elems
End ElementalData

Begin NodalData MAGNETOSTATIC_VECTOR_POTENTIAL
*Set cond Surface_Potential *nodes
*Add cond Line_Potential *nodes
*Add cond Point_Potential *nodes
*loop nodes  *OnlyInCond
*format "%i%i%f"
*NodesNum 1 *cond(X_Potential) *cond(Y_Potential) *cond(Z_Potential)
*end nodes
End NodalData

Begin Conditions Mfield3D
*Set cond Contour_Magnetic_Intensity *elems
*Add cond Infinit_Condition *elems
*loop elems *OnlyInCond
*ElemsNum *ElemsMat *GlobalNodes
*end elems
End Conditions

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




