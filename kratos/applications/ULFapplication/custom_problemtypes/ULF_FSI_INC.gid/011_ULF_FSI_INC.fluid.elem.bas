// Reading Elements

ElementsGroup = fluid_group;

*Set cond surface_UpdatedLagrangianFluid2D *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = UpdatedLagrangianFluid2D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems

*Set cond surface_UpdatedLagrangianFluid2Dinc *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = UpdatedLagrangianFluid2Dinc([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems

*Set cond volume_UpdatedLagrangianFluid3D *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = UpdatedLagrangianFluid3D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems

*Set cond volume_UpdatedLagrangianFluid3Dinc *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = UpdatedLagrangianFluid3Dinc([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems


ElementsGroup = solid_group;


*Set cond surface_MembraneElement *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = MembraneElement([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems

*Set cond surface_TotalLagrangian2D3N *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = TotalLagrangian2D3N([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems
*Set cond volume_TotalLagrangian3D4N *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = TotalLagrangian3D4N([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems
