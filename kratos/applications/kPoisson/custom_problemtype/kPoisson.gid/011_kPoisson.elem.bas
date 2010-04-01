// Reading Elements

ElementsGroup = dummy_group;

*loop elems
*Set var i=0
*set var j= ElemsNnode
*format "%i%i%i%i%i%i%i%i"
ELEMENTS[*ElemsNum] = Poisson2D([*\
*for(i=1;i<j;i=i+1)*\
*ElemsConec(*i),*\
*end*\
*ElemsConec(*ElemsNnode)],*ElemsMat);
*end elems
