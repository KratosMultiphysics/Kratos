*realformat "%10.5e"
*intformat "%i"
Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData

*# Tables
Begin Table 1 TIME DISPLACEMENT_X
*set var num_values(int)=GenData(Displacement-X_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Displacement-X_Table,*i,real) *GenData(Displacement-X_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 2 TIME DISPLACEMENT_Y
*set var num_values(int)=GenData(Displacement-Y_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Displacement-Y_Table,*i,real) *GenData(Displacement-Y_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 3 TIME DISPLACEMENT_Z
*set var num_values(int)=GenData(Displacement-Z_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Displacement-Z_Table,*i,real) *GenData(Displacement-Z_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 4 TIME TEMPERATURE
*set var num_values(int)=GenData(Temperature_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Temperature_Table,*i,real) *GenData(Temperature_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 5 TIME POINT_LOAD_X
*set var num_values(int)=GenData(Point-Load-X_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Point-Load-X_Table,*i,real) *GenData(Point-Load-X_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 6 TIME POINT_LOAD_Y
*set var num_values(int)=GenData(Point-Load-Y_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Point-Load-Y_Table,*i,real) *GenData(Point-Load-Y_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 7 TIME POINT_LOAD_Z
*set var num_values(int)=GenData(Point-Load-Z_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Point-Load-Z_Table,*i,real) *GenData(Point-Load-Z_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 8 TIME LINE_LOAD_X
*set var num_values(int)=GenData(Line-Load-X_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Line-Load-X_Table,*i,real) *GenData(Line-Load-X_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 9 TIME LINE_LOAD_Y
*set var num_values(int)=GenData(Line-Load-Y_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Line-Load-Y_Table,*i,real) *GenData(Line-Load-Y_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 10 TIME LINE_LOAD_Z
*set var num_values(int)=GenData(Line-Load-Z_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Line-Load-Z_Table,*i,real) *GenData(Line-Load-Z_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 11 TIME SURFACE_LOAD_X
*set var num_values(int)=GenData(Surface-Load-X_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Surface-Load-X_Table,*i,real) *GenData(Surface-Load-X_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 12 TIME SURFACE_LOAD_Y
*set var num_values(int)=GenData(Surface-Load-Y_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Surface-Load-Y_Table,*i,real) *GenData(Surface-Load-Y_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 13 TIME SURFACE_LOAD_Z
*set var num_values(int)=GenData(Surface-Load-Z_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Surface-Load-Z_Table,*i,real) *GenData(Surface-Load-Z_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 14 TIME NORMAL_LOAD
*set var num_values(int)=GenData(NormalLoad_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(NormalLoad_Table,*i,real) *GenData(NormalLoad_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 15 TIME MONTHS
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+1))
*end for
End Table

Begin Table 16 TIME WATER_LEVEL
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+2))
*end for
End Table

Begin Table 17 TIME OUTER_TEMPERATURE_BOFANG
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+3))
*end for
End Table

Begin Table 18 TIME REFERENCE_TEMEPERATURE
*set var num_values(int)=GenData(Evolution_Data,int)
*for(i=1;i<= num_values(int);i=i+5)
*GenData(Evolution_Data,*i,real) *GenData(Evolution_Data,*Operation(i+4))
*end for
End Table


*# Properties
*loop materials
*if(strcmp(MatProp(Element_Type),"Standard")==0 && (strcmp(MatProp(Standard_Constitutive_Law),"LinearElastic2DPlaneStress")==0 || strcmp(MatProp(Standard_Constitutive_Law),"LinearElastic2DPlaneStrain")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  *MatProp(Standard_Constitutive_Law)
DENSITY                *MatProp(Density,real)
YOUNG_MODULUS          *MatProp(Young_Modulus,real)
POISSON_RATIO          *MatProp(Poisson_Ratio,real)
THICKNESS              *MatProp(Thickness,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Standard")==0 && strcmp(MatProp(Standard_Constitutive_Law),"LinearElastic3D")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  LinearElastic3D
DENSITY                *MatProp(Density,real)
YOUNG_MODULUS          *MatProp(Young_Modulus,real)
POISSON_RATIO          *MatProp(Poisson_Ratio,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Standard")==0 && (strcmp(MatProp(Standard_Constitutive_Law),"ThermalLinearElastic2DPlaneStress")==0 || strcmp(MatProp(Standard_Constitutive_Law),"ThermalLinearElastic2DPlaneStrain")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME          *MatProp(Standard_Constitutive_Law)
DENSITY                        *MatProp(Density,real)
YOUNG_MODULUS                  *MatProp(Young_Modulus,real)
POISSON_RATIO                  *MatProp(Poisson_Ratio,real)
THICKNESS                      *MatProp(Thickness,real)
THERMAL_EXPANSION_COEFFICIENT  *MatProp(Thermal_Expansion_Coefficient,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Standard")==0 && strcmp(MatProp(Standard_Constitutive_Law),"ThermalLinearElastic3D")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME          ThermalLinearElastic3D
DENSITY                        *MatProp(Density,real)
YOUNG_MODULUS                  *MatProp(Young_Modulus,real)
POISSON_RATIO                  *MatProp(Poisson_Ratio,real)
THERMAL_EXPANSION_COEFFICIENT  *MatProp(Thermal_Expansion_Coefficient,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Interface")==0 && (strcmp(MatProp(Interface_Constitutive_Law),"BilinearCohesive2DPlaneStrain")==0 || strcmp(MatProp(Interface_Constitutive_Law),"BilinearCohesive2DPlaneStress")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME     BilinearCohesive2DLaw
DENSITY                   *MatProp(Density,real)
YOUNG_MODULUS             *MatProp(Young_Modulus,real)
POISSON_RATIO             *MatProp(Poisson_Ratio,real)
MINIMUM_JOINT_WIDTH       *MatProp(Minimum_Joint_Width,real)
CRITICAL_DISPLACEMENT     *MatProp(Critical_Displacement,real)
RESIDUAL_STRESS           *MatProp(Residual_Stress,real)
DAMAGE_THRESHOLD          *MatProp(Damage_Threshold,real)
FRICTION_COEFFICIENT      *MatProp(Friction_Coefficient,real)
THICKNESS                 *MatProp(Thickness,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Interface")==0 && (strcmp(MatProp(Interface_Constitutive_Law),"BilinearCohesive3D")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME     BilinearCohesive3DLaw
DENSITY                   *MatProp(Density,real)
YOUNG_MODULUS             *MatProp(Young_Modulus,real)
POISSON_RATIO             *MatProp(Poisson_Ratio,real)
MINIMUM_JOINT_WIDTH       *MatProp(Minimum_Joint_Width,real)
CRITICAL_DISPLACEMENT     *MatProp(Critical_Displacement,real)
RESIDUAL_STRESS           *MatProp(Residual_Stress,real)
DAMAGE_THRESHOLD          *MatProp(Damage_Threshold,real)
FRICTION_COEFFICIENT      *MatProp(Friction_Coefficient,real)
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

*Set cond Surface_Joint_Elements *elems
*if(CondNumEntities > 0)
Begin Elements SmallDisplacementInterfaceElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif

*Set cond Volume_Joint_Elements *elems
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==6)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements SmallDisplacementInterfaceElement3D6N
*loop elems *OnlyInCond
*if(ElemsNnode==6)
*ElemsNum  *ElemsMat *\
*for(i=1;i<=6;i=i+1)
 *ElemsConec(*i)*\
*end for

*endif
*end elems
End Elements

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==8)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements SmallDisplacementInterfaceElement3D8N
*loop elems *OnlyInCond
*if(ElemsNnode==8)
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*endif
*end elems
End Elements

*endif
*endif

*# Conditions
*# Note: CondId must be set to 0 before the first condition
*set var CondId=0
*Set cond point_Point_Load *nodes
*if(CondNumEntities > 0) 
*if(GenData(Domain_Size,int)==2)
Begin Conditions PointLoadCondition2D1N
*loop nodes *OnlyInCond
*set var CondId=CondId+1
*CondId  1  *NodesNum
*end nodes
End Conditions

*else
Begin Conditions PointLoadCondition3D1N
*loop nodes *OnlyInCond
*set var CondId=CondId+1
*CondId  1  *NodesNum
*end nodes
End Conditions

*endif
*endif
*Set cond line_Line_Load *elems *CanRepeat
*Add cond Line_Normal_Load *elems *CanRepeat
*Add cond Line_Uplift_Load *elems *CanRepeat
*if(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2)
Begin Conditions LineLoadCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic>0 && GenData(Domain_Size,int)==2)
Begin Conditions LineLoadCondition2D3N
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
*Add cond Surface_Normal_Load *elems *CanRepeat
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
Begin Conditions SurfaceLoadCondition3D3N
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
Begin Conditions SurfaceLoadCondition3D4N
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
Begin Conditions SurfaceLoadCondition3D6N
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
Begin Conditions SurfaceLoadCondition3D8N
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
Begin Conditions SurfaceLoadCondition3D9N
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
*set var freq = 0.5235987756
*set var t = GenData(Evolution_Data,2,real)
*set var t0 = cond(Day_Ambient_Temp,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
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
*set var freq = 0.5235987756
*set var t = GenData(Evolution_Data,2,real)
*set var t0 = cond(Day_Ambient_Temp,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
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
*set var freq = 0.5235987756
*set var t = GenData(Evolution_Data,2,real)
*set var t0 = cond(Day_Ambient_Temp,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
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
*NodesNum  *cond(Modify_X)  *cond(X_Value)
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
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
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
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
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
*NodesNum  *cond(Modify_X)  *cond(X_Value)
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
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
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
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif

*Set cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NEGATIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0)
*NodesNum  *cond(Modify_Normal)  *cond(Normal_Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*endif
*end nodes
End NodalData
*endif

*Set cond Surface_Normal_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NEGATIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0)
*NodesNum  *cond(Modify_Normal)  *cond(Normal_Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*endif
*end nodes
End NodalData
*endif

*Set cond Line_Uplift_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NEGATIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*if(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"X")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify_Normal)  *uplift
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
Warning, please change the direction!
*endif
*end nodes
End NodalData
*endif

*Set cond Surface_Uplift_Load *nodes
*if(CondNumEntities > 0)
Begin NodalData NEGATIVE_FACE_PRESSURE
*loop nodes *OnlyInCond
*if(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"X")==0 && strcmp(cond(Uplift_Direction),"X")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(3)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify_Normal)  *uplift
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
Warning, please change the direction!
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(1)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real) + GenData(Evolution_Data,3,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*set var Bas = cond(Base_of_Dam,real)
*set var upcoord = cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
*set var uplift = operation(pressure*(1.0-(1.0/Bas)*(abs(NodesCoord(2)-upcoord))))
*if(uplift > 0)
*NodesNum  *cond(Modify)  *uplift
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0 && strcmp(cond(Uplift_Direction),"Z")==0)
Warning, please change the direction!
*endif
*end nodes
End NodalData
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
*NodesNum  *cond(Modify_X)  *cond(X_Value)
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
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
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
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
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
*NodesNum  *cond(Fixed_Specific_Heat) *cond(Specific_Heat_Value)  
*end nodes
End NodalData

*endif
*Set cond volume_Mechanical_Conditions *nodes
*Add cond surface_Mechanical_Conditions *nodes
*if(CondNumEntities > 0)
Begin NodalData DENSITY
*loop nodes *OnlyinCond
*NodesNum  *cond(Fixed_Nodal_Density) *cond(Nodal_Density_Value)
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

*# Meshes
*Set cond volume_Solid_Displacement *nodes
*Add cond surface_Solid_Displacement *nodes
*Add cond line_Solid_Displacement *nodes
*Add cond point_Solid_Displacement *nodes
Begin Mesh 1 //DISPLACEMENT_X mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_X,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 2 //DISPLACEMENT_Y mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Y,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 3 //DISPLACEMENT_Z mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Z,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

*Set cond surface_Temperature *nodes
*Add cond line_Temperature *nodes
*Add cond point_Temperature *nodes
Begin Mesh 4 //TEMPERATURE mesh

	Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Temperature_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
*if(strcmp(cond(Fixed_Temperature),"1")==0)
 *NodesNum
*endif
*endif
*end nodes
*endif
	End MeshNodes

End Mesh

*Set cond point_Point_Load *nodes
Begin Mesh 5 //POINT_LOAD_X mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_X,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 6 //POINT_LOAD_Y mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Y,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 7 //POINT_LOAD_Z mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Z,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

*Set cond line_Line_Load *nodes
Begin Mesh 8 //LINE_LOAD_X mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_X,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 9 //LINE_LOAD_Y mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Y,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 10 //LINE_LOAD_Z mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Z,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

*Set cond surface_Surface_Load *nodes
Begin Mesh 11 //SURFACE_LOAD_X mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_X,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 12 //SURFACE_LOAD_Y mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Y,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

Begin Mesh 13 //SURFACE_LOAD_Z mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Z,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

*Set cond Surface_Normal_Load *nodes
*set var DoWrite=0
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
Begin Mesh 14 // Related to Uniform Normal Load Conditions

    Begin MeshNodes
*set var DoWrite=1
*break
*endif
*end nodes
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
    *NodesNum
*endif
*end nodes
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
    End MeshNodes
    
End Mesh
*break
*endif
*end nodes
*else
*Set cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
Begin Mesh 14 // Related to Uniform Normal Load Conditions

    Begin MeshNodes
*set var DoWrite=1
*break
*endif
*end nodes
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
    *NodesNum
*endif
*end nodes
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Uniform")==0 && (cond(Modify,int)==1))
    End MeshNodes
    
End Mesh
*break
*endif
*end nodes
*endif
*endif
*if(DoWrite==0)
Begin Mesh 14 // Related to Uniform Normal Load Conditions

    Begin MeshNodes
    End MeshNodes
    
End Mesh
*endif

*Set cond surface_Temperature *nodes
*Add cond line_Temperature *nodes
*Add cond point_Temperature *nodes
*if(CondNumEntities > 0)
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Bofang")==0  && (cond(Modify,int)==1)) 
Begin Mesh 15 // Related to Bofang Temperature Conditions

 Begin MeshData
  GRAVITY_DIRECTION        *cond(Gravity_Direction)
  COORDINATE_BASE_DAM      *cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real)
  SURFACE_TEMP             *cond(Surface_Temp,real)
  BOTTOM_TEMP              *cond(Bottom_Temp,real)
  HEIGHT_DAM               *cond(Height_Dam,real)
  AMPLITUDE                *cond(Temperature_Amplitude,real)
  DAY_MAXIMUM              *cond(Day_Ambient_Temp,real)
 End MeshData
 
    Begin MeshNodes
*break 
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Bofang")==0  && (cond(Modify,int)==1))
    *NodesNum
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Temperature_Distribution),"Bofang")==0  && (cond(Modify,int)==1))
    End MeshNodes
	
End Mesh 
*break
*endif
*end nodes 
*else
Begin Mesh 15 // Related to Bofang Temperature Conditions

 Begin MeshData
  GRAVITY_DIRECTION        Y
  COORDINATE_BASE_DAM      0.0
  SURFACE_TEMP             0.0
  BOTTOM_TEMP              0.0
  HEIGHT_DAM               0.0
  AMPLITUDE                0.0
  DAY_MAXIMUM              0.0 
 End MeshData
 
    Begin MeshNodes
    End MeshNodes
	
End Mesh   
*endif

*Set cond Surface_Normal_Load *nodes
*set var DoWrite=0
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0  && (cond(Modify_Normal,int)==1))
Begin Mesh 16 // Related to Hydrostatic Conditions
 
 Begin MeshData
  GRAVITY_DIRECTION        *cond(Gravity_Direction)
  COORDINATE_BASE_DAM      *cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real)
  SPECIFIC_WEIGHT          *cond(Specific_Weight,real)
 End MeshData

    Begin MeshNodes
*set var DoWrite=1
*break 
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0 && (cond(Modify_Normal,int)==1))
    *NodesNum
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0 && (cond(Modify_Normal,int)==1))
    End MeshNodes
    
End Mesh
*break 
*endif
*end nodes 
*else
*Set cond Line_Normal_Load *nodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0  && (cond(Modify_Normal,int)==1))
Begin Mesh 16 // Related to Hydrostatic Conditions
 
 Begin MeshData
  GRAVITY_DIRECTION        *cond(Gravity_Direction)
  COORDINATE_BASE_DAM      *cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real)
  SPECIFIC_WEIGHT          *cond(Specific_Weight,real)
 End MeshData

    Begin MeshNodes
*set var DoWrite=1
*break 
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0 && (cond(Modify_Normal,int)==1))
    *NodesNum
*endif
*end nodes
*loop nodes *OnlyinCond
*if(strcmp(cond(Normal_Load_Distribution),"Hydrostatic")==0 && (cond(Modify_Normal,int)==1))
    End MeshNodes
    
End Mesh
*break 
*endif
*end nodes 
*endif
*endif
*if(DoWrite==0)
Begin Mesh 16 // Related to Hydrostatic Conditions

 Begin MeshData
  GRAVITY_DIRECTION        Y
  COORDINATE_BASE_DAM      0.0
  SPECIFIC_WEIGHT          0.0
 End MeshData

    Begin MeshNodes
    End MeshNodes
 
End Mesh   
*endif

*Set cond Surface_Uplift_Load *nodes
*set var DoWrite=0
*if(CondNumEntities > 0 )
*loop nodes *OnlyInCond 
*if(cond(Modify,int)==1)
Begin Mesh 17 // Related to Uplift Pressure Conditions

 Begin MeshData
  GRAVITY_DIRECTION             *cond(Gravity_Direction)
  COORDINATE_BASE_DAM           *cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real)
  UPLIFT_DIRECTION              *cond(Uplift_Direction)
  COORDINATE_BASE_DAM_UPLIFT    *cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
  BASE_OF_DAM                   *cond(Base_of_Dam,real)
  SPECIFIC_WEIGHT               *cond(Specific_Weight,real)
 End MeshData
 
*set var DoWrite=1 
*break
*endif
*end nodes 
    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes
	
End Mesh
*else
*Set cond Line_Uplift_Load *nodes
*if(CondNumEntities > 0 )
*loop nodes *OnlyInCond 
*if(cond(Modify,int)==1)
Begin Mesh 17 // Related to Uplift Pressure Conditions

 Begin MeshData
  GRAVITY_DIRECTION             *cond(Gravity_Direction)
  COORDINATE_BASE_DAM           *cond(Reservoir_Bottom_Coordinate_in_Gravity_Direction,real)
  UPLIFT_DIRECTION              *cond(Uplift_Direction)
  COORDINATE_BASE_DAM_UPLIFT    *cond(Upstream_Coordinate_at_Base_Dam_in_Uplift_Direction,real)
  BASE_OF_DAM                   *cond(Base_of_Dam,real)
  SPECIFIC_WEIGHT               *cond(Specific_Weight,real)
 End MeshData
 
*set var DoWrite=1 
*break
*endif
*end nodes 
    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify,int)==1)
 *NodesNum
*endif
*end nodes
    End MeshNodes
	
End Mesh
*endif
*endif
*if(DoWrite==0)
Begin Mesh 17 // Related to Uplift Pressure Conditions

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
