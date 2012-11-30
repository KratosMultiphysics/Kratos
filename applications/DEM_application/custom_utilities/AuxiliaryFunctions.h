/* 
 * File:   GeometryFunctions.h
 * Author: msantasusana
 *
 * Created on 21 de mayo de 2012, 19:40
 */

#ifndef _KRATOSAUXILIARYFUNCTIONS_H
#define	_KRATOSAUXILIARYFUNCTIONS_H


#include <cmath>


namespace Kratos
{
    namespace AuxiliaryFunctions
    {

	static inline void CalculateAlphaFactor(int n_neighbours, double external_sphere_area, double total_equiv_area , double& alpha)
	{

	    double external_polyhedron_area = 0.0;
	    
	    switch (n_neighbours)
	    {
		case 4:
		    external_polyhedron_area = 3.30797*external_sphere_area;
		    break;
		case 5:
		    external_polyhedron_area = 2.60892*external_sphere_area;
		    break;
		case 6:
		    external_polyhedron_area = 1.90986*external_sphere_area;
		    break;
		case 7:
		    external_polyhedron_area = 1.78192*external_sphere_area;
		    break;
		case 8:
		    external_polyhedron_area = 1.65399*external_sphere_area;
		    break;
		case 9:
		    external_polyhedron_area = 1.57175*external_sphere_area;
		    break;
		case 10:
		    external_polyhedron_area = 1.48951*external_sphere_area;
		    break;
		case 11:
		    external_polyhedron_area = 1.40727*external_sphere_area;
		    break;
		case 12:
		    external_polyhedron_area = 1.32503*external_sphere_area;
		    break;
		case 13:
		    external_polyhedron_area = 1.31023*external_sphere_area;
		    break;
		case 14:
		    external_polyhedron_area = 1.29542*external_sphere_area;
		    break;
		case 15:
		    external_polyhedron_area = 1.28061*external_sphere_area;
		    break;
		case 16:
		    external_polyhedron_area = 1.26580*external_sphere_area;
		    break;
		case 17:
		    external_polyhedron_area = 1.25099*external_sphere_area;
		    break;
		case 18:
		    external_polyhedron_area = 1.23618*external_sphere_area;
		    break;
		case 19:
		    external_polyhedron_area = 1.22138*external_sphere_area;
		    break;
		case 20:
		    external_polyhedron_area = 1.20657*external_sphere_area;
		    break;
		default:
		    external_polyhedron_area = 1.15*external_sphere_area;
		    break;
	    }//switch (n_neighbours)
	      
	    alpha    = external_polyhedron_area/total_equiv_area;
	      
	}//CalculateAlphaFactor
  
    }
}
#endif	/* _KRATOSAUXILIARYFUNCTIONS_H */

