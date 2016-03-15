*realformat "%10.5e"
*intformat "%i"
Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

Begin Table 1 TIME MONTHS
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+1))
*end for
End Table

Begin Table 2 TIME WATER_LEVEL
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+2))
*end for
End Table

Begin Table 3 TIME OUTER_TEMPERATURE
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+3))
*end for
End Table

Begin Table 4 TIME REFERENCE_TEMEPERATURE
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+4))
*end for
End Table


*# Properties
*loop materials
*if(strcmp(MatProp(Constitutive_Law_Name),"LinearElastic2DPlaneStress")==0 || strcmp(MatProp(Constitutive_Law_Name),"LinearElastic2DPlaneStrain")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  *MatProp(Constitutive_Law_Name)
DENSITY                *MatProp(Density,real)
YOUNG_MODULUS          *MatProp(Young_Modulus,real)
POISSON_RATIO          *MatProp(Poisson_Ratio,real)
THICKNESS              *MatProp(Thickness,real)
End Properties

*elseif(strcmp(MatProp(Constitutive_Law_Name),"LinearElastic3D")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  LinearElastic3D
DENSITY                *MatProp(Density,real)
YOUNG_MODULUS          *MatProp(Young_Modulus,real)
POISSON_RATIO          *MatProp(Poisson_Ratio,real)
End Properties

*elseif(strcmp(MatProp(Constitutive_Law_Name),"ThermalLinearElastic2DPlaneStress")==0 || strcmp(MatProp(Constitutive_Law_Name),"ThermalLinearElastic2DPlaneStrain")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME          *MatProp(Constitutive_Law_Name)
DENSITY                        *MatProp(Density,real)
YOUNG_MODULUS                  *MatProp(Young_Modulus,real)
POISSON_RATIO                  *MatProp(Poisson_Ratio,real)
THICKNESS                      *MatProp(Thickness,real)
THERMAL_EXPANSION_COEFFICIENT  *MatProp(Thermal_Expansion_Coefficient,real)
End Properties

*elseif(strcmp(MatProp(Constitutive_Law_Name),"ThermalLinearElastic3D")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME          ThermalLinearElastic3D
DENSITY                        *MatProp(Density,real)
YOUNG_MODULUS                  *MatProp(Young_Modulus,real)
POISSON_RATIO                  *MatProp(Poisson_Ratio,real)
THERMAL_EXPANSION_COEFFICIENT  *MatProp(Thermal_Expansion_Coefficient,real)
End Properties

*endif
*end materials

*# Nodes
Begin Nodes
*loop nodes
*nodesnum  *NodesCoord(1) *NodesCoord(2) *NodesCoord(3)
*end nodes
End Nodes


*# Elements
*Set cond Triangles *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement2D3N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic>0 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement2D6N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=6;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Quadrilaterals *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==1 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement2D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==2 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement2D9N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=9;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Tetrahedra *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement3D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic>0 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement3D10N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=10;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Hexahedra *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement3D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==1 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement3D20N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=20;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==2 && strcmp(GenData(Deformation_Hypothesis),"Small_Displacements")==0)
Begin Elements SmallDisplacementThermoMechanicElement3D27N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=27;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*# Conditions
*# Note: CondId must be set to 0 before the first condition
*set var CondId=0
*Set cond point_Point_Load *nodes
*if(CondNumEntities > 0) 
*if(GenData(Domain_Size,int)==2)
Begin Conditions PointLoadCondition2D
*loop nodes *OnlyInCond
*set var CondId=CondId+1
*CondId  1  *NodesNum
*end nodes
End Conditions

*else
Begin Conditions PointLoadCondition3D
*loop nodes *OnlyInCond
*set var CondId=CondId+1
*CondId  1  *NodesNum
*end nodes
End Conditions

*endif
*endif
*Set cond line_Line_Load *elems *CanRepeat
*if(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2)
Begin Conditions LineLoadCondition2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic>0 && GenData(Domain_Size,int)==2)
Begin Conditions LineLoadCondition3N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*endif
*Set cond Line_Normal_Load *elems *CanRepeat
*Add cond Line_Uplift_Load *elems *CanRepeat
*if(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2)
Begin Conditions LineNormalLoadCondition2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic>0 && GenData(Domain_Size,int)==2)
Begin Conditions LineNormalLoadCondition3N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*endif
*Set cond surface_Surface_Load *elems *CanRepeat
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceLoadCondition3N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==4)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceLoadCondition4N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==4)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *GlobalNodes(*i)*\delta = self.delta_time
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==6)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceLoadCondition6N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==6)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=6;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==8)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceLoadCondition8N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==8)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==9)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceLoadCondition9N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==9)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=9;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*endif
*Set cond Surface_Normal_Load *elems *CanRepeat
*Add cond Surface_Uplift_Load *elems *CanRepeat
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceNormalLoadCondition3N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==4)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceNormalLoadCondition4N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==4)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==6)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceNormalLoadCondition6N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==6)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=6;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==8)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceNormalLoadCondition8N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==8)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==9)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceNormalLoadCondition9N
*loop elems *OnlyInCond
*if(ElemsNNodeFace==9)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=9;i=i+1)
 *GlobalNodes(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*endif
*# NodalData
*Set cond volume_Solid_Displacement *nodes
*Add cond surface_Solid_Displacement *nodes
*Add cond line_Solid_Displacement *nodes
*Add cond point_Solid_Displacement *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(SOLID_DISPLACEMENT_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData DISPLACEMENT_X
*loop nodes *OnlyInCond
*if(cond(SOLID_DISPLACEMENT_X,int)==1)
*NodesNum  *cond(Fix_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(SOLID_DISPLACEMENT_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData DISPLACEMENT_Y
*loop nodes *OnlyInCond
*if(cond(SOLID_DISPLACEMENT_Y,int)==1)
*NodesNum  *cond(Fix_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(SOLID_DISPLACEMENT_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData DISPLACEMENT_Z
*loop nodes *OnlyInCond
*if(cond(SOLID_DISPLACEMENT_Z,int)==1)
*NodesNum  *cond(Fix_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond volume_Temperature *nodes
*Add cond surface_Temperature *nodes
*Add cond line_Temperature *nodes
*Add cond point_Temperature *nodes
*if(CondNumEntities > 0)
*set var num_values(int)=GenData(Evolution_Data,int)
Begin NodalData TEMPERATURE
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Uniform")==0)
*NodesNum  *cond(Fixed_Temperature)  *cond(Temperature_Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var ts = cond(Surface_Temp,real)
*set var tb = cond(Bottom_Temp,real)
*set var H = cond(Height_Dam,real)
*set var A = cond(Temperature_Amplitude,real)
*set var freq = cond(Angular_Frequency,real)
*set var t = GenData(Evolution_Data,2,real)
*set var t0 = cond(Day_Ambient_Temp,real)
*set var RefCoord = cond(Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var aux = operation(RefCoord-NodesCoord(1))
*set var aux1 = operation((tb-(ts*exp(-0.04*H)))/(1-(exp(-0.04*H))))
*set var Temperature = operation(aux1+((ts-aux1)*(exp(-0.04*aux)))+(A*(exp(-0.018*aux))*(cos(freq*(t-(t0/30.0)-2.15+(1.30*exp(-0.085*aux)))))))
*if(aux >= 0)
*NodesNum   *cond(Fixed_Temperature)  *Temperature 
*else
*NodesNum   *cond(Fixed_Temperature)  *GenData(Evolution_Data,4,real) 
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)	
*set var ts = cond(Surface_Temp,real)
*set var tb = cond(Bottom_Temp,real)
*set var H = cond(Height_Dam,real)
*set var A = cond(Temperature_Amplitude,real)
*set var freq = cond(Angular_Frequency,real)
*set var t = GenData(Evolution_Data,2,real)
*set var t0 = cond(Day_Ambient_Temp,real)
*set var RefCoord = cond(Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var aux = operation(RefCoord-NodesCoord(2))
*set var aux1 = operation((tb-(ts*exp(-0.04*H)))/(1-(exp(-0.04*H))))
*set var Temperature = operation(aux1+((ts-aux1)*(exp(-0.04*aux)))+(A*(exp(-0.018*aux))*(cos(freq*(t-(t0/30.0)-2.15+(1.30*exp(-0.085*aux)))))))
*if(aux >= 0)
*NodesNum   *cond(Fixed_Temperature)  *Temperature 
*else
*NodesNum   *cond(Fixed_Temperature)  *GenData(Evolution_Data,4,real) 
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var ts = cond(Surface_Temp,real)
*set var tb = cond(Bottom_Temp,real)
*set var H = cond(Height_Dam,real)
*set var A = cond(Temperature_Amplitude,real)
*set var freq = cond(Angular_Frequency,real)
*set var t = GenData(Evolution_Data,2,real)
*set var t0 = cond(Day_Ambient_Temp,real)
*set var RefCoord = cond(Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var aux = operation(RefCoord-NodesCoord(3))
*set var aux1 = operation((tb-(ts*exp(-0.04*H)))/(1-(exp(-0.04*H))))
*set var Temperature = operation(aux1+((ts-aux1)*(exp(-0.04*aux)))+(A*(exp(-0.018*aux))*(cos(freq*(t-(t0/30.0)-2.15+(1.30*exp(-0.085*aux)))))))
*if(aux >= 0)
*NodesNum   *cond(Fixed_Temperature)  *Temperature 
*else
*NodesNum   *cond(Fixed_Temperature)  *GenData(Evolution_Data,4,real) 
*endif
*endif
*end nodes
End NodalData

*endif
*Set cond point_Point_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(POINT_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData POINT_LOAD_X
*loop nodes *OnlyInCond
*if(cond(POINT_LOAD_X,int)==1)
*NodesNum  *cond(Fix_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(POINT_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData POINT_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(POINT_LOAD_Y,int)==1)
*NodesNum  *cond(Fix_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(POINT_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData POINT_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(POINT_LOAD_Z,int)==1)
*NodesNum  *cond(Fix_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond line_Line_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData LINE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_X,int)==1)
*NodesNum  *cond(Fix_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData LINE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_Y,int)==1)
*NodesNum  *cond(Fix_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData LINE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(LINE_LOAD_Z,int)==1)
*NodesNum  *cond(Fix_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(LINE_NORMAL_LOAD,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData NORMAL_CONTACT_STRESS
*loop nodes *OnlyInCond
*if(cond(LINE_NORMAL_LOAD,int)==1)
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0)
*NodesNum  *cond(Fix_Normal)  *cond(Normal_Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Fix_Normal)  *Pressure
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Fix_Normal)  *Pressure
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*if(Pressure > 0)
*NodesNum  *cond(Fix_Normal)  *Pressure
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*endif
*endif
*end nodes
End NodalData
*endif
*endif

*Set cond Line_Uplift_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_CONTACT_STRESS
*loop nodes *OnlyInCond
*if(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fix_Normal)  *uplift
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fix_Normal)  *uplift
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"X")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fix_Normal)  *uplift
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fix_Normal)  *uplift
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fix_Normal)  *uplift
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fix_Normal)  *uplift
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
Warning, please change the direction!
*endif
*end nodes
End NodalData
*endif

*Set cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(LINE_TANGENTIAL_LOAD,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData TANGENTIAL_CONTACT_STRESS
*loop nodes *OnlyInCond
*if(cond(LINE_TANGENTIAL_LOAD,int)==1)
*NodesNum  *cond(Fix_Tangential)  *cond(Tangential_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond surface_Surface_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData SURFACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_X,int)==1)
*NodesNum  *cond(Fix_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData SURFACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_Y,int)==1)
*NodesNum  *cond(Fix_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData SURFACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(SURFACE_LOAD_Z,int)==1)
*NodesNum  *cond(Fix_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond Surface_Normal_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_CONTACT_STRESS
*loop nodes *OnlyInCond
*if(strcmp(cond(Load_Distribution),"Uniform")==0)
*NodesNum  *cond(Fixed)  *cond(Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Fixed)  *Pressure
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Fixed)  *Pressure
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*if(Pressure > 0)
*NodesNum  *cond(Fixed)  *Pressure
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*endif
*end nodes
End NodalData

*endif
*Set cond Surface_Uplift_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_CONTACT_STRESS
*loop nodes *OnlyInCond
*if(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fixed)  *uplift
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fixed)  *uplift
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"X")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fixed)  *uplift
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fixed)  *uplift
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fixed)  *uplift
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Fixed)  *uplift
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
Warning, please change the direction!
*endif
*end nodes
End NodalData

*endif
*Set cond volume_Thermal_Conditions *nodes
*Add cond surface_Thermal_Conditions *nodes
*if(CondNumEntities > 0)
Begin NodalData CONDUCTIVITY
*loop nodes *OnlyinCond
*NodesNum  *cond(Fixed_Conductivity)  *cond(Conductivity_Value)
*end nodes
End NodalData

Begin NodalData SPECIFIC_HEAT
*loop nodes *OnlyinCond
*NodesNum  *cond(Specific_Heat_Value)  *cond(Fixed_Specific_Heat)
*end nodes
End NodalData

*endif
*Set cond volume_Mechanical_Conditions *nodes
*Add cond surface_Mechanical_Conditions *nodes
*if(CondNumEntities > 0)
Begin NodalData DENSITY
*loop nodes *OnlyinCond
*NodesNum  *cond(Nodal_Density_Value)  *cond(Fixed_Nodal_Density)
*end nodes
End NodalData

*endif
*Set cond volume_Body_Acceleration *nodes
*Add cond surface_Body_Acceleration *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(BODY_ACCELERATION_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData VOLUME_ACCELERATION_X
*loop nodes *OnlyInCond
*if(cond(BODY_ACCELERATION_X,int)==1)
*NodesNum  *cond(Fix_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(BODY_ACCELERATION_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData VOLUME_ACCELERATION_Y
*loop nodes *OnlyInCond
*if(cond(BODY_ACCELERATION_Y,int)==1)
*NodesNum  *cond(Fix_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(BODY_ACCELERATION_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData VOLUME_ACCELERATION_Z
*loop nodes *OnlyInCond
*if(cond(BODY_ACCELERATION_Z,int)==1)
*NodesNum  *cond(Fix_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData
*endif
*endif

*Set cond surface_Temperature *nodes
*Add cond line_Temperature *nodes
*Add cond point_Temperature *nodes
*if(CondNumEntities > 0)
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Bofang")==0)
Begin Mesh 1 // Related to Bofang Temperature Conditions

 Begin MeshData
  GRAVITY_DIRECTION        *cond(Gravity_Direction)
  COORDINATE_BASE_DAM      *cond(Coordinate_at_Base_Dam_in_Gravity_Direction,real)
  SURFACE_TEMP             *cond(Surface_Temp,real)
  BOTTOM_TEMP              *cond(Bottom_Temp,real)
  HEIGHT_DAM               *cond(Height_Dam,real)
  AMPLITUDE                *cond(Temperature_Amplitude,real)
  FREQUENCY                *cond(Angular_Frequency,real)
  DAY_MAXIMUM              *cond(Day_Ambient_Temp,real)
 End MeshData
 
 Begin MeshNodes
*break 
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Bofang")==0)
 *NodesNum
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Bofang")==0)
 End MeshNodes
End Mesh 
*break
*endif
*end nodes 

*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Uniform")==0)
Begin Mesh 2 // Related to Uniform  Temperature Conditions

 Begin MeshNodes
*break 
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Uniform")==0)
*if(strcmp(cond(Fixed_Temperature),"1")==0)
 *NodesNum
*endif
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Uniform")==0)
*if(strcmp(cond(Fixed_Temperature),"1")==0)
 End MeshNodes
End Mesh
*break 
*endif
*endif
*end nodes
*else
Begin Mesh 1 // Related to Bofang Temperature Conditions

 Begin MeshData
  GRAVITY_DIRECTION        Y
  COORDINATE_BASE_DAM      0.0
  SURFACE_TEMP             0.0
  BOTTOM_TEMP              0.0
  HEIGHT_DAM               0.0
  AMPLITUDE                0.0
  FREQUENCY                0.0
  DAY_MAXIMUM              0.0
 End MeshData
 
 Begin MeshNodes
 End MeshNodes
End Mesh

Begin Mesh 2 // Related to Uniform  Temperature Conditions

 Begin MeshNodes
 End MeshNodes
End Mesh   
*endif

*Set cond Surface_Normal_Load *nodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Load_Distribution),"Hydrostatic")==0)
Begin Mesh 3 // Related to Hydrostatic Conditions
 Begin MeshData

  GRAVITY_DIRECTION        *cond(Gravity_Direction)
  COORDINATE_BASE_DAM      *cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real)
  SPECIFIC_WEIGHT          *cond(Specific_Weight,real)
 End MeshData

*break 
*endif
*end nodes
 Begin MeshNodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Load_Distribution),"Hydrostatic")==0)
 *NodesNum
*endif
*end nodes
 End MeshNodes
End Mesh  
*endif

*Set cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0)
Begin Mesh 3 // Related to Hydrostatic Conditions

 Begin MeshData
  GRAVITY_DIRECTION        *cond(Gravity_Direction)
  COORDINATE_BASE_DAM      *cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real)
  SPECIFIC_WEIGHT          *cond(Specific_Weight,real)
 End MeshData

*break 
*endif
*end nodes
 Begin MeshNodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0)
 *NodesNum
*endif
*end nodes
 End MeshNodes
End Mesh  
*endif

*Set cond Surface_Normal_Load *nodes
*Add cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
*else
Begin Mesh 3 // Related to Hydrostatic Conditions

 Begin MeshData
  GRAVITY_DIRECTION        Y
  COORDINATE_BASE_DAM      0.0
  SPECIFIC_WEIGHT          0.0
 End MeshData
 
 Begin MeshNodes
 End MeshNodes
End Mesh   
*endif

*Set cond Line_Uplift_Load *nodes
*Add cond Surface_Uplift_Load *nodes
*if(CondNumEntities > 0 )
*loop nodes *OnlyInCond 
Begin Mesh 4 // Related to Uplift Pressure Conditions

 Begin MeshData
  GRAVITY_DIRECTION             *cond(Gravity_Direction)
  COORDINATE_BASE_DAM           *cond(Upstream_Coordinate_at_Base_Dam_in_Gravity_Direction,real)
  UPLIFT_DIRECTION              *cond(Uplift_Direction)
  COORDINATE_BASE_DAM_UPLIFT    *cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
  BASE_OF_DAM                   *cond(Base_of_Dam,real)
  SPECIFIC_WEIGHT               *cond(Specific_Weight,real)
 End MeshData
 
*break
*end nodes 
 Begin MeshNodes
*loop nodes *OnlyInCond
 *NodesNum
*end nodes
 End MeshNodes
End Mesh
*else
Begin Mesh 4 // Related to Uplift Pressure Conditions

 Begin MeshData
  GRAVITY_DIRECTION             Y
  COORDINATE_BASE_DAM           0.0
  UPLIFT_DIRECTION              X
  COORDINATE_BASE_DAM_UPLIFT    0.0
  BASE_OF_DAM                   0.0
  SPECIFIC_WEIGHT               0.0
 End MeshData
 
 Begin MeshNodes
 End MeshNodes
End Mesh
*endif



