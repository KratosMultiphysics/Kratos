/* 
 * File:   GeometryFunctions.h
 * Author: msantasusana
 *
 * Created on 21 de mayo de 2012, 19:40
 */

#ifndef _DEM_AUXILIARY_FUNCTIONS_H
#define	_DEM_AUXILIARY_FUNCTIONS_H


#include <cmath>
#include "containers/array_1d.h"
#include "../../../kratos/includes/serializer.h"

namespace Kratos {
    
    namespace AuxiliaryFunctions {

	static inline void CalculateAlphaFactor3D(int n_neighbours, double external_sphere_area, double total_equiv_area , double& alpha) {

	    double external_polyhedron_area = 0.0;
	    switch (n_neighbours) {
                  
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
			
		alpha = external_polyhedron_area/total_equiv_area;
			
	}//CalculateAlphaFactor
	  
	static inline void CalculateAlphaFactor2D(int n_neighbours, double external_sphere_perimeter, double total_equiv_perimeter , double& alpha)
        {

            double external_polyhedron_perimeter = 0.0;
        
            switch (n_neighbours) {
              
                case 3:
                    external_polyhedron_perimeter = 1.65399*external_sphere_perimeter;
                    break;
                  
                case 4:
                    external_polyhedron_perimeter = 1.27324*external_sphere_perimeter;
                    break;
              
                case 5:
                    external_polyhedron_perimeter = 1.15633*external_sphere_perimeter;
                    break;
                  
                case 6:
                    external_polyhedron_perimeter = 1.10266*external_sphere_perimeter;
                    break;
              
                case 7:
                    external_polyhedron_perimeter = 1.07303*external_sphere_perimeter;
                    break;
              
                case 8:
                    external_polyhedron_perimeter = 1.05479*external_sphere_perimeter;
                    break;
              
                case 9:
                    external_polyhedron_perimeter = 1.04270*external_sphere_perimeter;
                    break;
              
                case 10:
                    external_polyhedron_perimeter = 1.03425*external_sphere_perimeter;
                    break;
              
                case 11:
                    external_polyhedron_perimeter = 1.02811*external_sphere_perimeter;
                    break;
              
                case 12:
                    external_polyhedron_perimeter = 1.02349*external_sphere_perimeter;
                    break;
              
                case 13:
                    external_polyhedron_perimeter = 1.01993*external_sphere_perimeter;
                    break;
              
                case 14:
                    external_polyhedron_perimeter = 1.01713*external_sphere_perimeter;
                    break;
                     
                default:
                    external_polyhedron_perimeter = 1.0*external_sphere_perimeter;
                    break;
                       
            }//switch (n_neighbours)
            
            alpha = external_polyhedron_perimeter/total_equiv_perimeter;
            
        }//CalculateAlphaFactor
	  
	  
	static inline void SwitchCase(int case_opt, bool& delta_OPTION, bool& continuum_simulation_OPTION) {
		
	    switch (case_opt) {
                
                case 0:
                    delta_OPTION = false;
	            continuum_simulation_OPTION = false;
	            break;
                      
	        case 1:
                    delta_OPTION = true;
		    continuum_simulation_OPTION = false;
	            break;
                      
	        case 2:
		    delta_OPTION = true;
		    continuum_simulation_OPTION = true;
		    break;
                      
		case 3:
		    delta_OPTION = false;
		    continuum_simulation_OPTION = true;
	            break;
                      
		default:
                    delta_OPTION = false;
	            continuum_simulation_OPTION = false;
				  
	    }
			  
	} //SwitchCase
	  
	inline array_1d<double,3> LinearTimeIncreasingFunction(array_1d<double,3> external_total_applied_force, double current_time, double final_time)
	{		
            array_1d<double,3> externally_applied_force_now = external_total_applied_force*current_time/final_time;
            return externally_applied_force_now;
        
        }// inline array_1d<double,3> LinearTimeIncreasingFunction
 	  
    }//namespace AuxiliaryFunctions

  
}//namespace Kratos

#endif	/* _KRATOSAUXILIARYFUNCTIONS_H */

