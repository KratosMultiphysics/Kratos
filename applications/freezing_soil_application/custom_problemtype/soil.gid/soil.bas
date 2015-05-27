//problemtype for KRATOS freezing soil application
//(c) 2013 Meng-Meng Zhou, Ruhr-University Bochum

Begin ModelPartData
//  VARIABLE_NAME value
End ModelPartData
*loop materials
*if(strcmp(MatProp(ConstitutiveLaw),"UserDefined")==0)
Begin Properties *MatNum
*if(strcmp(GenData(Enable_Gravity),"1")==0)
GRAVITY [3] ( *GenData(Gravity_X,real), *GenData(Gravity_Y,real), *GenData(Gravity_Z,real) )
*else
GRAVITY [3] ( 0.0, 0.0, 0.0 )
*endif
End Properties
*else
Begin Properties *MatNum
DENSITY *MatProp(Density)  //scalar
*set var bodyForceX = Operation(MatProp(Density,real)*GenData(Gravity_X,real))
*set var bodyForceY = Operation(MatProp(Density,real)*GenData(Gravity_Y,real))
*set var bodyForceZ = Operation(MatProp(Density,real)*GenData(Gravity_Z,real))
*if(strcmp(GenData(Enable_Gravity),"1")==0)
BODY_FORCE [3] (*bodyForceX, *bodyForceY, *bodyForceZ)
GRAVITY [3] ( *GenData(Gravity_X,real), *GenData(Gravity_Y,real), *GenData(Gravity_Z,real) )
*else
BODY_FORCE [3] (0.0, 0.0, 0.0)
GRAVITY [3] ( 0.0, 0.0, 0.0 )
*endif
*if(strcmp(MatProp(ConstitutiveLaw),"Isotropic3D")==0)
DENSITY *MatProp(Density,real)
YOUNG_MODULUS *MatProp(Young_modulus,real)
POISSON_RATIO *MatProp(Poisson_ratio,real)
*endif
*if(strcmp(MatProp(ConstitutiveLaw),"TutorialDamageModel")==0)
DENSITY *MatProp(Density,real)
YOUNG_MODULUS *MatProp(Young_modulus,real)
POISSON_RATIO *MatProp(Poisson_ratio,real)
DAMAGE_E0 *MatProp(E0,real)
DAMAGE_EF *MatProp(Ef,real)
*endif
*if(strcmp(MatProp(ConstitutiveLaw),"GroutingMortar")==0)
DENSITY *MatProp(Density,real)
YOUNG_MODULUS *MatProp(Young_modulus,real)
POISSON_RATIO *MatProp(Poisson_ratio,real)
PRIMARY_HYDRATION_TIME *MatProp(prim_hyd_time,real)
PRIMARY_HYDRATION_TIME_GRADIENT *MatProp(gradient_prim_hyd_time,real)
STIFFNESS_RATIO *MatProp(E_ratio,real)
*endif
*if(strcmp(MatProp(ConstitutiveLaw),"TrussMaterial")==0)
DENSITY *MatProp(Density,real)
YOUNG_MODULUS *MatProp(Young_modulus,real)
POISSON_RATIO *MatProp(Poisson_ratio,real)
*endif 
End Properties
*endif
*end materials

Begin Nodes
*RealFormat "%10.5f"
*loop nodes
*nodesnum   *NodesCoord(1)  *NodesCoord(2)  *NodesCoord(3)
*end nodes
End Nodes

*include soil_elements.include


//Begin Conditions Condition2D
//1799 0        644        650
//1800 0        650        663
//1801 0        663        673
//...
//1947 0        972        973
//1948 0        973        974
//End Conditions


*set var condID=1
*set var condID2=1

Begin Conditions MasterContactFace3D3
*set cond Contact_Master
*loop elems *OnlyInCond
*if(ElemsNnodeFace==3)
*condID      *ElemsMat    *\
*GlobalNodes(2) *GlobalNodes(1) *GlobalNodes(3) 
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions MasterContactFace3D6
*set cond Contact_Master
*loop elems *OnlyInCond
*if(ElemsNnodeFace==6)
*condID      *ElemsMat    *\
*Globalnodes(2) *Globalnodes(1) *Globalnodes(3) *Globalnodes(4) *Globalnodes(6) *Globalnodes(5)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions MasterContactFace3D
*set cond Contact_Master
*loop elems *OnlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions MasterContactFace3D8
*set cond Contact_Master
*loop elems *OnlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions MasterContactFace3D9
*set cond Contact_Master
*loop elems *OnlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions SlaveContactFace3D3
*set cond Contact_Slave
*loop elems *OnlyInCond
*if(ElemsNnodeFace==3)
*condID      *ElemsMat    *\
*GlobalNodes(2) *GlobalNodes(1) *GlobalNodes(3) 
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions SlaveContactFace3D6
*set cond Contact_Slave
*loop elems *OnlyInCond
*if(ElemsNnodeFace==6)
*condID      *ElemsMat    *\
*Globalnodes(2) *Globalnodes(1) *Globalnodes(3) *Globalnodes(4) *Globalnodes(6) *Globalnodes(5)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions SlaveContactFace3D
*set cond Contact_Slave
*loop elems *OnlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions SlaveContactFace3D8
*set cond Contact_Slave
*loop elems *OnlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions SlaveContactFace3D9
*set cond Contact_Slave
*loop elems *OnlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions PointForce3D
*set cond Point_Force *nodes
*add cond Line_Force *nodes
*add cond Surface_Force *nodes
*add cond Volume_Force *nodes
*loop nodes *OnlyInCond
*condID 1 *NodesNum
*set var condID= condID+1
*end nodes
End Conditions

Begin Conditions FaceForce3D3N
*set cond Distributed_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==3)
*condID      *ElemsMat    *\
*GlobalNodes(2) *GlobalNodes(1) *GlobalNodes(3) 
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceForce3D6N
*set cond Distributed_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==6)
*condID      *ElemsMat    *\
*Globalnodes(2) *Globalnodes(1) *Globalnodes(3) *Globalnodes(4) *Globalnodes(6) *Globalnodes(5)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceForce3D4N
*set cond Distributed_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceForce3D8N
*set cond Distributed_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceForce3D9N
*set cond Distributed_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceWaterFlux3D4N
*set cond Distributed_Surface_Water_Flux *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceWaterFlux3D8N
*set cond Distributed_Surface_Water_Flux *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceWaterFlux3D9N
*set cond Distributed_Surface_Water_Flux *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions


Begin Conditions FaceHeatFlux3D4N
*set cond Distributed_Surface_Heat_Flux *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceHeatFlux3D8N
*set cond Distributed_Surface_Heat_Flux *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceHeatFlux3D9N
*set cond Distributed_Surface_Heat_Flux *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions


Begin Conditions FaceHeatConvection3D4N
*set cond Distributed_Surface_Heat_Convection *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceHeatConvection3D8N
*set cond Distributed_Surface_Heat_Convection *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceHeatConvection3D9N
*set cond Distributed_Surface_Heat_Convection *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions


Begin Conditions FaceHeatRadiation3D4N
*set cond Distributed_Surface_Heat_Radiation *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceHeatRadiation3D8N
*set cond Distributed_Surface_Heat_Radiation *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceHeatRadiation3D9N
*set cond Distributed_Surface_Heat_Radiation *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions


Begin Conditions FaceLoadPressure3D4N
*set cond Surface_Load_Pressure *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceLoadPressure3D8N
*set cond Surface_Load_Pressure *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions FaceLoadPressure3D9N
*set cond Surface_Load_Pressure *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions


Begin Conditions Face3D3N
*set cond Following_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==3)
*condID      *ElemsMat    *\
*GlobalNodes(2) *GlobalNodes(1) *GlobalNodes(3) 
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions Face3D6N
*set cond Following_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==6)
*condID      *ElemsMat    *\
*Globalnodes(2) *Globalnodes(1) *Globalnodes(3) *Globalnodes(4) *Globalnodes(6) *Globalnodes(5)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions Face3D4N
*set cond Following_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions Face3D8N
*set cond Following_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions Face3D9N
*set cond Following_Surface_Load *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions LiningEndFace3D3N
*set cond Lining_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==3)
*condID      *ElemsMat    *\
*GlobalNodes(2) *GlobalNodes(1) *GlobalNodes(3) 
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions LiningEndFace3D6N
*set cond Lining_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==6)
*condID      *ElemsMat    *\
*Globalnodes(2) *Globalnodes(1) *Globalnodes(3) *Globalnodes(4) *Globalnodes(6) *Globalnodes(5)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions LiningEndFace3D4N
*set cond Lining_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions LiningEndFace3D8N
*set cond Lining_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions LiningEndFace3D9N
*set cond Lining_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions TBMEndFace3D3N
*set cond TBM_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==3)
*condID      *ElemsMat    *\
*GlobalNodes(2) *GlobalNodes(1) *GlobalNodes(3) 
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions TBMEndFace3D6N
*set cond TBM_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==6)
*condID      *ElemsMat    *\
*Globalnodes(2) *Globalnodes(1) *Globalnodes(3) *Globalnodes(4) *Globalnodes(6) *Globalnodes(5)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions TBMEndFace3D4N
*set cond TBM_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==4)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions TBMEndFace3D8N
*set cond TBM_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==8)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin Conditions TBMEndFace3D9N
*set cond TBM_End *elems *canRepeat
*loop elems *onlyInCond
*if(ElemsNnodeFace==9)
*condID      *ElemsMat    *\
*Globalnodes(4) *Globalnodes(3) *Globalnodes(2) *Globalnodes(1) *Globalnodes(7) *Globalnodes(6) *Globalnodes(5) *Globalnodes(8) *Globalnodes(9)
//ElementAssignment *condID *ElemsNum
*set var condID= condID+1
*endif
*end elems
End Conditions

Begin NodalData DISPLACEMENT_X
*set cond Nodal_Displacement *nodes *or(1,int)
*add cond Line_Displacement *nodes *or(1,int)
*add cond Surface_Displacement *nodes *or(1,int)
*add cond Volume_Displacement *nodes *or(1,int)
*loop nodes *OnlyInCond
*NodesNum   *cond(1)   *cond(2)
*end nodes
*set cond IsSide *nodes
*loop nodes *OnlyInCond
*NodesNum   1   0.0
*end nodes
End NodalData

Begin NodalData DISPLACEMENT_Y
*set cond Nodal_Displacement *nodes *or(3,int)
*add cond Line_Displacement *nodes *or(3,int)
*add cond Surface_Displacement *nodes *or(3,int)
*add cond Volume_Displacement *nodes *or(3,int)
*loop nodes *OnlyInCond
*NodesNum   *cond(3)   *cond(4)
*end nodes
*set cond IsSide *nodes
*loop nodes *OnlyInCond
*NodesNum   1   0.0
*end nodes
End NodalData

Begin NodalData DISPLACEMENT_Z
*set cond Nodal_Displacement *nodes *or(5,int)
*add cond Line_Displacement *nodes *or(5,int)
*add cond Surface_Displacement *nodes *or(5,int)
*add cond Volume_Displacement *nodes *or(5,int)
*loop nodes *OnlyInCond
*NodesNum   *cond(5)   *cond(6)
*end nodes
*set cond IsBottom *nodes
*loop nodes *OnlyInCond
*NodesNum   1   0.0
*end nodes
End NodalData 

Begin NodalData WATER_PRESSURE
*set cond Initial_Nodal_Water_Pressure *nodes
*add cond Initial_Line_Water_Pressure *nodes
*add cond Initial_Surface_Water_Pressure *nodes
*add cond Initial_Body_Water_Pressure *nodes
*loop nodes *OnlyInCond
*NodesNum   0   *cond(1)  
*end nodes
*set cond Nodal_Water_Pressure *nodes *or(1,int)
*add cond Line_Water_Pressure *nodes *or(1,int)
*add cond Surface_Water_Pressure *nodes *or(1,int)
*add cond Volume_Water_Pressure *nodes *or(1,int)
*loop nodes *OnlyInCond
*NodesNum   *cond(1)   *cond(2)
*end nodes
End NodalData
 

Begin NodalData TEMPERATURE
*set cond Initial_Nodal_Temperature *nodes
*add cond Initial_Line_Temperature *nodes
*add cond Initial_Surface_Temperature *nodes
*add cond Initial_Body_Temperature *nodes
*loop nodes *OnlyInCond
*NodesNum   0   *cond(1)  
*end nodes
*set cond Nodal_Temperature *nodes *or(1,int)
*add cond Line_Temperature *nodes *or(1,int)
*add cond Surface_Temperature *nodes *or(1,int)
*add cond Volume_Temperature *nodes *or(1,int)
*loop nodes *OnlyInCond
*NodesNum   *cond(1)   *cond(2)
*end nodes
End NodalData
 


Begin NodalData NEGATIVE_FACE_PRESSURE
*set cond Following_Surface_Load *nodes
*loop nodes *OnlyInCond
*NodesNum 0 *cond(2)
*end nodes
End NodalData

Begin NodalData FACE_LOAD_PRESSURE_X
*set cond Surface_Load_Pressure *nodes
*loop nodes *OnlyInCond
*NodesNum 0 *cond(1)
*end nodes
End NodalData

Begin NodalData FACE_LOAD_PRESSURE_Y
*set cond Surface_Load_Pressure *nodes
*loop nodes *OnlyInCond
*NodesNum 0 *cond(2)
*end nodes
End NodalData

Begin NodalData FACE_LOAD_PRESSURE_Z
*set cond Surface_Load_Pressure *nodes
*loop nodes *OnlyInCond
*NodesNum 0 *cond(3)
*end nodes
End NodalData
 

Begin NodalData FACE_LOAD
*set cond Distributed_Surface_Load *nodes 
*loop nodes *OnlyInCond
*NodesNum 0 [3] ( *cond(1), *cond(2), *cond(3) )
*end nodes
End NodalData


Begin NodalData FACE_WATER_FLUX
*set cond Distributed_Surface_Water_Flux *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *cond(1)
*end nodes
End NodalData

Begin NodalData FACE_HEAT_FLUX
*set cond Distributed_Surface_Heat_Flux *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *cond(1)
*end nodes
End NodalData

Begin NodalData CONVECTION_COEFFICIENT
*set cond Distributed_Surface_Heat_Convection *nodes 
*loop nodes *OnlyInCond
*NodesNum 1 *cond(1)
*end nodes
End NodalData

Begin NodalData EMISSIVITY
*set cond Distributed_Surface_Heat_Radiation *nodes 
*loop nodes *OnlyInCond
*NodesNum 1 *cond(1)
*end nodes
End NodalData

Begin NodalData AMBIENT_TEMPERATURE
*set cond Distributed_Surface_Heat_Convection *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *cond(2)
*end nodes
*set cond Distributed_Surface_Heat_Radiation *nodes
*loop nodes *OnlyInCond
*NodesNum 1 *cond(2)
*end nodes
End NodalData
 
Begin ElementalData AREA
*set cond LineElementType *elems
*loop elems *OnlyInCond
*if((strcmp(cond(1),"Truss")==0) || (strcmp(cond(1),"Beam")==0))
*ElemsNum *cond(2)
*endif
*end elems
End ElementalData

Begin ElementalData INERTIA
*set cond LineElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"Beam")==0)
*ElemsNum [3,3] (( *cond(5), 0, 0 ),(0.0, *cond(3), 0.0),(0.0, 0.0, *cond(4)))
*endif
*end elems
End ElementalData

Begin ElementalData ACTIVATION_LEVEL
*set cond Line_Activation_Level *elems
*add cond Surface_Activation_Level *elems
*add cond Volume_Activation_Level *elems
*format "%i%i"
*loop elems *OnlyInCond
*ElemsNum *cond(1)
*end elems
End ElementalData

Begin ElementalData DENSITY_WATER
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*ElemsNum *cond(2)
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*ElemsNum *cond(2)
*endif
*if(strcmp(cond(1),"Grouting_Element")==0)
*ElemsNum *cond(2)
*endif
*end elems
End ElementalData

Begin ElementalData DENSITY_AIR
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*ElemsNum *cond(3)
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*ElemsNum *cond(3)
*endif
*end elems
End ElementalData

Begin ElementalData POROSITY
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*ElemsNum *cond(4)
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*ElemsNum *cond(4)
*endif 
*end elems
End ElementalData

Begin ElementalData PERMEABILITY_WATER
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*format "%i%14.12e"
*ElemsNum *cond(5)
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*format "%i%14.12e"
*ElemsNum *cond(5)
*endif 
*end elems
End ElementalData

Begin ElementalData PERMEABILITY_AIR
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*format "%i%14.12e"
*ElemsNum *cond(6)
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*format "%i%14.12e"
*ElemsNum *cond(6)
*endif
*end elems
End ElementalData

Begin ElementalData FIRST_SATURATION_PARAM
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*ElemsNum 2.5
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*ElemsNum 2.5
*endif 
*end elems
End ElementalData

Begin ElementalData SECOND_SATURATION_PARAM
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*ElemsNum 0.4
*endif 
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*ElemsNum 0.4
*endif
*end elems
End ElementalData

Begin ElementalData AIR_ENTRY_VALUE
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"UnsaturatedSoil_2Phase")==0)
*ElemsNum 3000.0
*endif
*if(strcmp(cond(1),"UnsaturatedSoil_3Phase")==0)
*ElemsNum 3000.0
*endif
*end elems
End ElementalData

Begin ElementalData PERMEABILITY_28_DAYS
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"Grouting_Element")==0)
*ElemsNum *cond(7)
*endif
*end elems
End ElementalData

Begin ElementalData PERMEABILITY_1_DAY
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"Grouting_Element")==0)
*ElemsNum *cond(8)
*endif
*end elems
End ElementalData

Begin ElementalData PERMEABILITY_TRANSITION
*set cond VolumeElementType *elems
*loop elems *OnlyInCond
*if(strcmp(cond(1),"Grouting_Element")==0)
*ElemsNum *cond(9)
*endif
*end elems
End ElementalData

Begin ElementalData GRAVITY
*loop elems
*if(strcmp(GenData(Enable_Gravity),"1")==0)
*ElemsNum [3] (*GenData(Gravity_X,real), *GenData(Gravity_Y,real), *GenData(Gravity_Z,real))
*endif
*end elems
End ElementalData
