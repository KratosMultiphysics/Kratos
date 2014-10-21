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

namespace Kratos
{
    
class PropertiesProxy {
            
  public:
      
      int     GetId()                                                          { return mId;                                 }
      void    SetId(int id)                                                    { mId = id;                                   }
       
      double  GetYoung()                                                       { return *mYoung;                             }
      double* pGetYoung()                                                      { return  mYoung;                             }
      void    SetYoungFromProperties(double* young)                            { mYoung = young;                             }  
      
      double  GetPoisson()                                                     { return *mPoisson;                           }
      double* pGetPoisson()                                                    { return  mPoisson;                           }
      void    SetPoissonFromProperties(double* poisson)                        { mPoisson = poisson;                         }                    
            
      double  GetRollingFriction()                                             { return *mRollingFriction;                   }      
      double* pGetRollingFriction()                                            { return  mRollingFriction;                   }  
      void    SetRollingFrictionFromProperties(double* rolling_friction)       { mRollingFriction = rolling_friction;        }      
      
      double  GetTgOfFrictionAngle()                                           { return *mTgOfFrictionAngle;                 } 
      double* pGetTgOfFrictionAngle()                                          { return  mTgOfFrictionAngle;                 } 
      void    SetTgOfFrictionAngleFromProperties(double* tg_of_friction_angle) { mTgOfFrictionAngle = tg_of_friction_angle;  }
      
      double  GetLnOfRestitCoeff()                                             { return *mLnOfRestitCoeff;                   } 
      double* pGetLnOfRestitCoeff()                                            { return  mLnOfRestitCoeff;                   } 
      void    SetLnOfRestitCoeffFromProperties(double* ln_of_restit_coeff)     { mLnOfRestitCoeff = ln_of_restit_coeff;      }  
      
      double  GetDensity()                                                     { return *mDensity;                           }
      double* pGetDensity()                                                    { return  mDensity;                           }
      void    SetDensityFromProperties(double* density)                        { mDensity = density;                         }  
      
      int     GetParticleMaterial()                                            { return *mParticleMaterial;                  }
      int*    pGetParticleMaterial()                                           { return  mParticleMaterial;                  }
      void    SetParticleMaterialFromProperties(int* particle_material)        { mParticleMaterial = particle_material;      }  
      
      PropertiesProxy operator=(PropertiesProxy props){
          
          mId                = props.GetId();
          mYoung             = props.pGetYoung();
          mPoisson           = props.pGetPoisson();
          mRollingFriction   = props.pGetRollingFriction();
          mTgOfFrictionAngle = props.pGetTgOfFrictionAngle();
          mLnOfRestitCoeff   = props.pGetLnOfRestitCoeff();
          mDensity           = props.pGetDensity();
          mParticleMaterial  = props.pGetParticleMaterial();
                  
          return *this;
      }
           
  private:
      int     mId;
      double* mYoung;
      double* mPoisson;
      double* mRollingFriction;
      double* mTgOfFrictionAngle;
      double* mLnOfRestitCoeff;
      double* mDensity;
      int*    mParticleMaterial;
      
      
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          rSerializer.save("mId",mId);
          rSerializer.save("mYoung",mYoung);
          rSerializer.save("mPoisson",mPoisson);
          rSerializer.save("mRollingFriction",mRollingFriction);
          rSerializer.save("mTgOfFrictionAngle",mTgOfFrictionAngle);
          rSerializer.save("mLnOfRestitCoeff",mLnOfRestitCoeff);
          rSerializer.save("mDensity",mDensity);
          rSerializer.save("mParticleMaterial",mParticleMaterial);
      }

      virtual void load(Serializer& rSerializer)
      {
          rSerializer.load("mId",mId);
          rSerializer.load("mYoung",mYoung);
          rSerializer.load("mPoisson",mPoisson);
          rSerializer.load("mRollingFriction",mRollingFriction);
          rSerializer.load("mTgOfFrictionAngle",mTgOfFrictionAngle);
          rSerializer.load("mLnOfRestitCoeff",mLnOfRestitCoeff);
          rSerializer.load("mDensity",mDensity);
          rSerializer.load("mParticleMaterial",mParticleMaterial);          
      }
      
};
  

    namespace AuxiliaryFunctions
    {

	  static inline void CalculateAlphaFactor3D(int n_neighbours, double external_sphere_area, double total_equiv_area , double& alpha)
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
	  
	  static inline void CalculateAlphaFactor2D(int n_neighbours, double external_sphere_perimeter, double total_equiv_perimeter , double& alpha)
      {

          double external_polyhedron_perimeter = 0.0;
          
          switch (n_neighbours)
          {
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
            
          alpha    = external_polyhedron_perimeter/total_equiv_perimeter;
            
      }//CalculateAlphaFactor
	  
	  
	  static inline void SwitchCase(int case_opt, bool& delta_OPTION, bool& continuum_simulation_OPTION)
	  {
		
		  switch (case_opt)
			
			{
			  
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

