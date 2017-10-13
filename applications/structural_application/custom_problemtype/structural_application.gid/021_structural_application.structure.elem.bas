// Reading Elements

ElementsGroup = structural_group;

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
*Set cond surface_TotalLagrangian2D4N *elems
*loop elems *onlyInCond
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = TotalLagrangian2D4N([*\
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
