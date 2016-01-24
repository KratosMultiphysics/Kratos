*realformat "%10.5e"
*intformat "%i"
Begin ModelPartData
//VARIABLE_NAME value
End ModelPartData


*# Properties
*loop materials
*if(strcmp(MatProp(Constitutive_Law_Name),"LinearElastic2DPlaneStrain")==0 || strcmp(MatProp(Constitutive_Law_Name),"LinearElastic2DPlaneStress")==0)
Begin Properties  *MatNum
CONSTITUTIVE_LAW_NAME  *MatProp(Constitutive_Law_Name)
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

*elseif(strcmp(MatProp(Constitutive_Law_Name),"LinearElastic3D")==0)
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
Begin Elements SmallStrainUPwElement2D3N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=3;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements SmallStrainUPwFICElement2D3N
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
Begin Elements SmallStrainUPwElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements SmallStrainUPwFICElement2D4N
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
*Set cond Quadrilateral_Interfaces *elems
*if(CondNumEntities > 0 && IsQuadratic==0)
Begin Elements SmallStrainUPwInterfaceElement2D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Tetrahedra *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Elements SmallStrainUPwElement3D4N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=4;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements SmallStrainUPwFICElement3D4N
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
*Set cond Prismatic_Interfaces *elems
*if(CondNumEntities > 0 && IsQuadratic==0)
Begin Elements SmallStrainUPwInterfaceElement3D6N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=6;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*endif
*Set cond Hexahedra *elems
*if(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"False")==0)
Begin Elements SmallStrainUPwElement3D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

*elseif(CondNumEntities > 0 && IsQuadratic==0 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Elements SmallStrainUPwFICElement3D8N
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
*Set cond Hexahedral_Interfaces *elems
*if(CondNumEntities > 0 && IsQuadratic==0)
Begin Elements SmallStrainUPwInterfaceElement3D8N
*loop elems *OnlyInCond
*ElemsNum  *ElemsMat *\
*for(i=1;i<=8;i=i+1)
 *ElemsConec(*i)*\
*end for

*end elems
End Elements

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
*Set cond volume_Fluid_Pressure *nodes
*Add cond surface_Fluid_Pressure *nodes
*Add cond line_Fluid_Pressure *nodes
*Add cond point_Fluid_Pressure *nodes
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
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(1)))
*if(Pressure > 0)
*NodesNum  *cond(Fix_Normal)  *Pressure
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Y")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
*set var Pressure = operation(SpecWei*(RefCoord-NodesCoord(2)))
*if(Pressure > 0)
*NodesNum  *cond(Fix_Normal)  *Pressure
*else
*NodesNum  *cond(Fix_Normal)  0.0
*endif
*elseif(strcmp(cond(Gravity_Direction),"Z")==0)
*set var SpecWei = cond(Specific_Weight,real)
*set var RefCoord = cond(Reference_Coordinate,real)
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
*Set cond Line_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_FLUID_FLUX
*loop nodes *OnlyInCond
*NodesNum  *cond(Fixed)  *cond(Value)
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
*Set cond Surface_Normal_Fluid_Flux *nodes
*if(CondNumEntities > 0)
Begin NodalData NORMAL_FLUID_FLUX
*loop nodes *OnlyInCond
*NodesNum  *cond(Fixed)  *cond(Value)
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
Begin Conditions LineNormalLoadCondition2D2N
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
Begin Conditions LineNormalFluidFluxCondition2D2N
*loop elems *OnlyInCond
*set var CondId=CondId+1
*CondId  *ElemsMat *\
*for(i=1;i<=2;i=i+1)
 *GlobalNodes(*i)*\
*end for

*end elems
End Conditions

*elseif(CondNumEntities > 0 && IsQuadratic==0 && GenData(Domain_Size,int)==2 && strcmp(GenData(FIC_Stabilization),"True")==0)
Begin Conditions LineNormalFluidFluxFICCondition2D2N
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
Begin Conditions SurfaceNormalLoadCondition3D3N
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
Begin Conditions SurfaceNormalLoadCondition3D4N
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
Begin Conditions SurfaceNormalFluidFluxCondition3D3N
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
Begin Conditions SurfaceNormalFluidFluxCondition3D4N
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
Begin Conditions SurfaceNormalFluidFluxFICCondition3D3N
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
Begin Conditions SurfaceNormalFluidFluxFICCondition3D4N
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
