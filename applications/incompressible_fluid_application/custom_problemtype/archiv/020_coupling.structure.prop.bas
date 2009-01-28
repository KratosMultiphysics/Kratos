PROPERTIES[1][END_TIME] =   *GenData(END_TIME);
PROPERTIES[1][DELTA_TIME] = *GenData(DELTA_TIME);

*loop materials
*if(strcmp(MatProp(ID),"Structure")==0)
*format "%i%f"
PROPERTIES[*MatNum][DENSITY] = *MatProp(Density,real); 
*format "%i%f"
PROPERTIES[*MatNum][YOUNG_MODULUS] = *MatProp(Young_Module,real);
*format "%i%f"
PROPERTIES[*MatNum][POISSON_RATIO] = *MatProp(Poisson_Ratio,real);
*format "%i%f"
PROPERTIES[*MatNum][THICKNESS] = *MatProp(Thickness,real);
*format "%i%f"
PROPERTIES[*MatNum][CROSS_SECTION] = *MatProp(Cross_Section_Area,real);
*format "%i%f"

*endif
*end materials


