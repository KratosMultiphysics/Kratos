*loop materials
*if(strcmp(MatProp(1),"Solid")==0)
*format "%i%f"
PROPERTIES[*MatNum][DENSITY] = *MatProp(Density,real); 
*format "%i%f"
PROPERTIES[*MatNum][YOUNG_MODULUS] = *MatProp(YoungModulus,real); 
*format "%i%f"
PROPERTIES[*MatNum][POISSON_RATIO] = *MatProp(PoissonRatio,real); 
*format "%i%f"
PROPERTIES[*MatNum][THICKNESS] = *MatProp(Thickness,real); 
*format "%i%f"
PROPERTIES[*MatNum][BODY_FORCE] = [0.000, *operation(MatProp(Gravity,real)*MatProp(Density,real)), 0.000];
*endif
*end materials

