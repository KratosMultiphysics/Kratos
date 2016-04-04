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

Begin Table 4 TIME WATER_PRESSURE
*set var num_values(int)=GenData(Pressure_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Pressure_Table,*i,real) *GenData(Pressure_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 5 TIME FORCE_X
*set var num_values(int)=GenData(Force-X_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Force-X_Table,*i,real) *GenData(Force-X_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 6 TIME FORCE_Y
*set var num_values(int)=GenData(Force-Y_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Force-Y_Table,*i,real) *GenData(Force-Y_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 7 TIME FORCE_Z
*set var num_values(int)=GenData(Force-Z_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Force-Z_Table,*i,real) *GenData(Force-Z_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 8 TIME FACE_LOAD_X
*set var num_values(int)=GenData(Face-Load-X_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Face-Load-X_Table,*i,real) *GenData(Face-Load-X_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 9 TIME FACE_LOAD_Y
*set var num_values(int)=GenData(Face-Load-Y_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Face-Load-Y_Table,*i,real) *GenData(Face-Load-Y_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 10 TIME FACE_LOAD_Z
*set var num_values(int)=GenData(Face-Load-Z_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Face-Load-Z_Table,*i,real) *GenData(Face-Load-Z_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 11 TIME NORMAL_CONTACT_STRESS
*set var num_values(int)=GenData(Normal-Load_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Normal-Load_Table,*i,real) *GenData(Normal-Load_Table,*Operation(i+1),real,real)
*end for
End Table

Begin Table 12 TIME TANGENTIAL_CONTACT_STRESS
*set var num_values(int)=GenData(Tangential-Load_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Tangential-Load_Table,*i,real) *GenData(Tangential-Load_Table,*Operation(i+1),real)
*end for
End Table

Begin Table 13 TIME NORMAL_FLUID_FLUX
*set var num_values(int)=GenData(Normal-Fluid-Flux_Table,int)
*for(i=1;i<= num_values(int);i=i+2)
*GenData(Normal-Fluid-Flux_Table,*i,real) *GenData(Normal-Fluid-Flux_Table,*Operation(i+1),real)
*end for
End Table


*# Properties
*loop materials
*if(strcmp(MatProp(Element_Type),"Standard")==0 && (strcmp(MatProp(Standard_Constitutive_Law),"LinearElastic2DPlaneStrain")==0 || strcmp(MatProp(Standard_Constitutive_Law),"LinearElastic2DPlaneStress")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  *MatProp(Standard_Constitutive_Law)
YOUNG_MODULUS          *MatProp(Young_Modulus,real)
POISSON_RATIO          *MatProp(Poisson_Ratio,real)
DENSITY_SOLID          *MatProp(Solid_Density,real)
DENSITY_WATER          *MatProp(Fluid_Density,real)
POROSITY               *MatProp(Porosity,real)
BULK_MODULUS_SOLID     *MatProp(Solid_Bulk_Modulus,real)
BULK_MODULUS_FLUID     *MatProp(Fluid_Bulk_Modulus,real)
PERMEABILITY_XX        *MatProp(Intrinsic_Permeability_XX,real)
PERMEABILITY_YY        *MatProp(Intrinsic_Permeability_YY,real)
PERMEABILITY_XY        *MatProp(Intrinsic_Permeability_XY,real)
DYNAMIC_VISCOSITY      *MatProp(Dynamic_Viscosity,real)
THICKNESS              *MatProp(Thickness,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Standard")==0 && (strcmp(MatProp(Standard_Constitutive_Law),"LinearElastic3D")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  LinearElastic3D
YOUNG_MODULUS          *MatProp(Young_Modulus,real)
POISSON_RATIO          *MatProp(Poisson_Ratio,real)
DENSITY_SOLID          *MatProp(Solid_Density,real)
DENSITY_WATER          *MatProp(Fluid_Density,real)
POROSITY               *MatProp(Porosity,real)
BULK_MODULUS_SOLID     *MatProp(Solid_Bulk_Modulus,real)
BULK_MODULUS_FLUID     *MatProp(Fluid_Bulk_Modulus,real)
PERMEABILITY_XX        *MatProp(Intrinsic_Permeability_XX,real)
PERMEABILITY_YY        *MatProp(Intrinsic_Permeability_YY,real)
PERMEABILITY_ZZ        *MatProp(Intrinsic_Permeability_ZZ,real)
PERMEABILITY_XY        *MatProp(Intrinsic_Permeability_XY,real)
PERMEABILITY_YZ        *MatProp(Intrinsic_Permeability_YZ,real)
PERMEABILITY_ZX        *MatProp(Intrinsic_Permeability_ZX,real)
DYNAMIC_VISCOSITY      *MatProp(Dynamic_Viscosity,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Interface")==0 && (strcmp(MatProp(Interface_Constitutive_Law),"BilinearCohesive2DPlaneStrain")==0 || strcmp(MatProp(Interface_Constitutive_Law),"BilinearCohesive2DPlaneStress")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME     BilinearCohesive2D
YOUNG_MODULUS             *MatProp(Young_Modulus,real)
POISSON_RATIO             *MatProp(Poisson_Ratio,real)
DENSITY_SOLID             *MatProp(Solid_Density,real)
DENSITY_WATER             *MatProp(Fluid_Density,real)
POROSITY                  *MatProp(Porosity,real)
BULK_MODULUS_SOLID        *MatProp(Solid_Bulk_Modulus,real)
BULK_MODULUS_FLUID        *MatProp(Fluid_Bulk_Modulus,real)
TRANSVERSAL_PERMEABILITY  *MatProp(Transversal_Permeability,real)
DYNAMIC_VISCOSITY         *MatProp(Dynamic_Viscosity,real)
MINIMUM_JOINT_WIDTH       *MatProp(Minimum_Joint_Width,real)
CRITICAL_DISPLACEMENT     *MatProp(Critical_Displacement,real)
RESIDUAL_STRESS           *MatProp(Residual_Stress,real)
DAMAGE_THRESHOLD          *MatProp(Damage_Threshold,real)
FRICTION_COEFFICIENT      *MatProp(Friction_Coefficient,real)
THICKNESS                 *MatProp(Thickness,real)
End Properties

*elseif(strcmp(MatProp(Element_Type),"Interface")==0 && (strcmp(MatProp(Interface_Constitutive_Law),"BilinearCohesive3D")==0))
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME     BilinearCohesive3D
YOUNG_MODULUS             *MatProp(Young_Modulus,real)
POISSON_RATIO             *MatProp(Poisson_Ratio,real)
DENSITY_SOLID             *MatProp(Solid_Density,real)
DENSITY_WATER             *MatProp(Fluid_Density,real)
POROSITY                  *MatProp(Porosity,real)
BULK_MODULUS_SOLID        *MatProp(Solid_Bulk_Modulus,real)
BULK_MODULUS_FLUID        *MatProp(Fluid_Bulk_Modulus,real)
TRANSVERSAL_PERMEABILITY  *MatProp(Transversal_Permeability,real)
DYNAMIC_VISCOSITY         *MatProp(Dynamic_Viscosity,real)
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
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Elements UPwSmallStrainElement2D3N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements UPwSmallStrainFICElement2D3N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic>0)
Begin Elements SmallStrainUPwDiffOrderElement2D6N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=6;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Quadrilaterals *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Elements UPwSmallStrainElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements UPwSmallStrainFICElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==1)
Begin Elements SmallStrainUPwDiffOrderElement2D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==2)
Begin Elements SmallStrainUPwDiffOrderElement2D9N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=9;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Surface_Interface_Elements *elems
*if(CondNumEntities > 0)
Begin Elements UPwSmallStrainInterfaceElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Surface_Link_Interface_Elements *elems
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements UPwSmallStrainLinkInterfaceElement2D4N
*loop elems *OnlyInCond
*if(ElemsNnode==3)
*ElemsNum  *ElemsMat  *ElemsConec(1) *ElemsConec(2) *ElemsConec(2) *ElemsConec(3)
*endif
*end elems
End Elements

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements UPwSmallStrainLinkInterfaceElement2D4N
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*endif
*end elems
End Elements

*endif
*endif
*Set cond Tetrahedra *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Elements UPwSmallStrainElement3D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements UPwSmallStrainFICElement3D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic>0)
Begin Elements SmallStrainUPwDiffOrderElement3D10N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=10;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Hexahedra *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Elements UPwSmallStrainElement3D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements UPwSmallStrainFICElement3D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==1)
Begin Elements SmallStrainUPwDiffOrderElement3D20N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=20;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==2)
Begin Elements SmallStrainUPwDiffOrderElement3D27N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=27;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Volume_Interface_Elements *elems
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==6)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements UPwSmallStrainInterfaceElement3D6N
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
Begin Elements UPwSmallStrainInterfaceElement3D8N
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
*Set cond Volume_Link_Interface_Elements *elems
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements UPwSmallStrainLinkInterfaceElement3D6N
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*ElemsNum  *ElemsMat  *ElemsConec(1) *ElemsConec(2) *ElemsConec(3) *ElemsConec(4) *ElemsConec(2) *ElemsConec(3)
*endif
*end elems
End Elements

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==6)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Elements UPwSmallStrainLinkInterfaceElement3D6N
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
Begin Elements UPwSmallStrainLinkInterfaceElement3D8N
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
*# Note: CondId must be set to 0 before the first condition (the first condition must have ID = 1)
*set var CondId=0
*Set cond Point_Force *nodes
*if(CondNumEntities > 0)
*if(GenData(Domain_Size,int)==2)
Begin Conditions UPwForceCondition2D1N
*loop nodes *OnlyInCond
*set var CondId=CondId+1
*CondId  1  *NodesNum
*end nodes
End Conditions

*else
Begin Conditions UPwForceCondition3D1N
*loop nodes *OnlyInCond
*set var CondId=CondId+1
*CondId  1  *NodesNum
*end nodes
End Conditions

*endif
*endif
*Set cond Line_Face_Load *elems *CanRepeat
*if(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2)
Begin Conditions UPwFaceLoadCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic>0 && GenData(Domain_Size,int)==2)
Begin Conditions LineLoadDiffOrderCondition2D3N
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
*if(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2)
Begin Conditions UPwNormalFaceLoadCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic>0 && GenData(Domain_Size,int)==2)
Begin Conditions LineNormalLoadDiffOrderCondition2D3N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*endif
*Set cond Line_Normal_Fluid_Flux *elems *CanRepeat
*if(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Conditions UPwNormalFluxCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Conditions UPwNormalFluxFICCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic>0 && GenData(Domain_Size,int)==2)
Begin Conditions LineNormalFluidFluxDiffOrderCondition2D3N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*endif
*Set cond Line_Interface_Face_Load *elems *CanRepeat
*if(CondNumEntities > 0 && GenData(Domain_Size,int)==2)
Begin Conditions UPwFaceLoadInterfaceCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Conditions

*endif
*Set cond Line_Interface_Normal_Fluid_Flux *elems *CanRepeat
*if(CondNumEntities > 0 && GenData(Domain_Size,int)==2)
Begin Conditions UPwNormalFluxInterfaceCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Conditions

*endif
*Set cond Surface_Face_Load *elems *CanRepeat
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwFaceLoadCondition3D3N
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
Begin Conditions UPwFaceLoadCondition3D4N
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
Begin Conditions SurfaceLoadDiffOrderCondition3D6N
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
Begin Conditions SurfaceLoadDiffOrderCondition3D8N
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
Begin Conditions SurfaceLoadDiffOrderCondition3D9N
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
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwNormalFaceLoadCondition3D3N
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
Begin Conditions UpwNormalFaceLoadCondition3D4N
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
Begin Conditions SurfaceNormalLoadDiffOrderCondition3D6N
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
Begin Conditions SurfaceNormalLoadDiffOrderCondition3D8N
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
Begin Conditions SurfaceNormalLoadDiffOrderCondition3D9N
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
*Set cond Surface_Normal_Fluid_Flux *elems *CanRepeat
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwNormalFluxCondition3D3N
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
Begin Conditions UPwNormalFluxCondition3D4N
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
*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwNormalFluxFICCondition3D3N
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
Begin Conditions UPwNormalFluxFICCondition3D4N
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
*elseif(CondNumEntities > 0 && IsQuadratic>0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNNodeFace==6)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D6N
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
Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D8N
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
Begin Conditions SurfaceNormalFluidFluxDiffOrderCondition3D9N
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
*Set cond Surface_Interface_Face_Load *elems *CanRepeat
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwFaceLoadInterfaceCondition3D4N
*loop elems *OnlyInCond
*if(ElemsNnode==3)
*set var CondId=CondId+1
*CondId  *ElemsMat  *ElemsConec(1) *ElemsConec(2) *ElemsConec(2) *ElemsConec(3)
*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwFaceLoadInterfaceCondition3D4N
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*endif
*Set cond Surface_Interface_Normal_Fluid_Flux *elems *CanRepeat
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==3)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwNormalFluxInterfaceCondition3D4N
*loop elems *OnlyInCond
*if(ElemsNnode==3)
*set var CondId=CondId+1
*CondId  *ElemsMat  *ElemsConec(1) *ElemsConec(2) *ElemsConec(2) *ElemsConec(3)
*endif
*end elems
End Conditions

*endif
*set var DoWrite=0
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*set var DoWrite=1
*break
*endif
*end elems
*if(DoWrite==1)
Begin Conditions UPwNormalFluxInterfaceCondition3D4N
*loop elems *OnlyInCond
*if(ElemsNnode==4)
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*endif
*end elems
End Conditions

*endif
*endif

*# NodalData
*Set cond Volume_Solid_Displacement *nodes
*Add cond Surface_Solid_Displacement *nodes
*Add cond Line_Solid_Displacement *nodes
*Add cond Point_Solid_Displacement *nodes
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
*Set cond Volume_Fluid_Pressure *nodes
*Add cond Surface_Fluid_Pressure *nodes
*Add cond Line_Fluid_Pressure *nodes
*Add cond Point_Fluid_Pressure *nodes
*if(CondNumEntities > 0)
Begin NodalData WATER_PRESSURE
*loop nodes *OnlyInCond
*if(strcmp(cond(Pressure_Distribution),"Uniform")==0)
*NodesNum  *cond(Fixed)  *cond(Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Fixed)  *Pressure
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Fixed)  *Pressure
*else
*NodesNum  *cond(Fixed)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
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
*Set cond Point_Force *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FORCE_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FORCE_X
*loop nodes *OnlyInCond
*if(cond(FORCE_X,int)==1)
*NodesNum  *cond(Modify_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FORCE_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FORCE_Y
*loop nodes *OnlyInCond
*if(cond(FORCE_Y,int)==1)
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FORCE_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FORCE_Z
*loop nodes *OnlyInCond
*if(cond(FORCE_Z,int)==1)
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond Line_Face_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*NodesNum  *cond(Modify_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
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
*NodesNum  *cond(Modify_Normal)  *cond(Normal_Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*if(Pressure > 0)
*NodesNum  *cond(Modify_Normal)  *Pressure
*else
*NodesNum  *cond(Modify_Normal)  0.0
*endif
*endif
*endif
*end nodes
End NodalData

*endif
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
*NodesNum  *cond(Modify_Tangential)  *cond(Tangential_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond Line_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_FLUID_FLUX
*loop nodes *OnlyInCond
*NodesNum  *cond(Modify)  *cond(Value)
*end nodes
End NodalData

*endif
*Set cond Line_Interface_Face_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*NodesNum  *cond(Modify_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond Line_Interface_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_FLUID_FLUX
*loop nodes *OnlyInCond
*NodesNum  *cond(Modify)  *cond(Value)
*end nodes
End NodalData

*endif
*Set cond Surface_Face_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*NodesNum  *cond(Modify_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
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
*NodesNum  *cond(Modify)  *cond(Value)
*elseif(strcmp(cond(Gravity_Direction),"X")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Modify)  *Pressure
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Modify)  *Pressure
*else
*NodesNum  *cond(Modify)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(3)))
*if(Pressure > 0)
*NodesNum  *cond(Modify)  *Pressure
*else
*NodesNum  *cond(Modify)  0.0
*endif
*endif
*end nodes
End NodalData

*endif
*Set cond Surface_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_FLUID_FLUX
*loop nodes *OnlyInCond
*NodesNum  *cond(Modify)  *cond(Value)
*end nodes
End NodalData

*endif
*Set cond Surface_Interface_Face_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_X
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_X,int)==1)
*NodesNum  *cond(Modify_X)  *cond(X_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Y
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Y,int)==1)
*NodesNum  *cond(Modify_Y)  *cond(Y_Value)
*endif
*end nodes
End NodalData

*endif
*set var DoWrite=0
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*set var DoWrite=1
*break
*endif
*end nodes
*if(DoWrite == 1)
Begin NodalData FACE_LOAD_Z
*loop nodes *OnlyInCond
*if(cond(FACE_LOAD_Z,int)==1)
*NodesNum  *cond(Modify_Z)  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif
*Set cond Surface_Interface_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_FLUID_FLUX
*loop nodes *OnlyInCond
*NodesNum  *cond(Modify)  *cond(Value)
*end nodes
End NodalData

*endif
*Set cond Volume_Body_Acceleration *nodes
*Add cond Surface_Body_Acceleration *nodes
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
*NodesNum  0  *cond(X_Value)
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
*NodesNum  0  *cond(Y_Value)
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
*NodesNum  0  *cond(Z_Value)
*endif
*end nodes
End NodalData

*endif
*endif

*# Meshes
*Set cond Volume_Solid_Displacement *nodes
*Add cond Surface_Solid_Displacement *nodes
*Add cond Line_Solid_Displacement *nodes
*Add cond Point_Solid_Displacement *nodes
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

*Set cond Volume_Fluid_Pressure *nodes
*Add cond Surface_Fluid_Pressure *nodes
*Add cond Line_Fluid_Pressure *nodes
*Add cond Point_Fluid_Pressure *nodes
Begin Mesh 4 //WATER_PRESSURE mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

*Set cond Point_Force *nodes
Begin Mesh 5 //FORCE_X mesh

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

Begin Mesh 6 //FORCE_Y mesh

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

Begin Mesh 7 //FORCE_Z mesh

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

*set var DoWrite=0
*Set cond Surface_Face_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=1
Begin Mesh 8 //FACE_LOAD_X mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_X,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh

Begin Mesh 9 //FACE_LOAD_Y mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_Y,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh

Begin Mesh 10 //FACE_LOAD_Z mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_Z,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh
*endif
*Set cond Line_Face_Load *nodes
*if(DoWrite == 0 && CondNumEntities > 0)
*set var DoWrite=1
Begin Mesh 8 //FACE_LOAD_X mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_X,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh

Begin Mesh 9 //FACE_LOAD_Y mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_Y,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh

Begin Mesh 10 //FACE_LOAD_Z mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_Z,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh
*elseif(DoWrite == 0)
*set var DoWrite=1
Begin Mesh 8 //FACE_LOAD_X mesh

    Begin MeshNodes
    End MeshNodes

End Mesh

Begin Mesh 9 //FACE_LOAD_Y mesh

    Begin MeshNodes
    End MeshNodes

End Mesh

Begin Mesh 10 //FACE_LOAD_Z mesh

    Begin MeshNodes
    End MeshNodes

End Mesh
*endif

*set var DoWrite=0
*Set cond Surface_Normal_Load *nodes
*if(CondNumEntities > 0)
*set var DoWrite=1
Begin Mesh 11 //NORMAL_CONTACT_STRESS mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh
*endif
*Set cond Line_Normal_Load *nodes
*if(DoWrite == 0 && CondNumEntities > 0)
*set var DoWrite=1
Begin Mesh 11 //NORMAL_CONTACT_STRESS mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify_Normal,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh
*elseif(DoWrite == 0)
*set var DoWrite=1
Begin Mesh 11 //NORMAL_CONTACT_STRESS mesh

    Begin MeshNodes
    End MeshNodes

End Mesh
*endif

Begin Mesh 12 //TANGENTIAL_CONTACT_STRESS mesh

    Begin MeshNodes
*if(CondNumEntities > 0)
*loop nodes *OnlyInCond
*if(cond(Modify_Tangential,int)==1)
    *NodesNum
*endif
*end nodes
*endif
    End MeshNodes

End Mesh

*set var DoWrite=0
*Set cond Surface_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
*set var DoWrite=1
Begin Mesh 13 //NORMAL_FLUID_FLUX mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh
*endif
*Set cond Line_Normal_Fluid_Flux *nodes
*if(DoWrite == 0 && CondNumEntities > 0)
*set var DoWrite=1
Begin Mesh 13 //NORMAL_FLUID_FLUX mesh

    Begin MeshNodes
*loop nodes *OnlyInCond
*if(cond(Modify,int)==1)
    *NodesNum
*endif
*end nodes
    End MeshNodes

End Mesh
*elseif(DoWrite == 0)
*set var DoWrite=1
Begin Mesh 13 //NORMAL_FLUID_FLUX mesh

    Begin MeshNodes
    End MeshNodes

End Mesh
*endif
