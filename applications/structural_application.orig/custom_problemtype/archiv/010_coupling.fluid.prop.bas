PROPERTIES[1][END_TIME] =   *GenData(END_TIME);
PROPERTIES[1][DELTA_TIME] = *GenData(DELTA_TIME);

*loop materials
*if(strcmp(MatProp(ID),"Fluid")==0)
*format "%i%f"
PROPERTIES[*MatNum][DENSITY] = *MatProp(Density,real); 
*format "%i%f"
PROPERTIES[*MatNum][VISCOSITY] = *MatProp(Viscosity,real);
*format "%i%f"
PROPERTIES[*MatNum][BODY_FORCE] = [0.000, *MatProp(gravity,real), 0.000];
*format "%i%f"
PROPERTIES[*MatNum][TIME_PARAMETER] = 0.5;
*format "%i%f"
PROPERTIES[*MatNum][STABILIZATION_FACTOR] = 1.000000;

*endif
*end materials
