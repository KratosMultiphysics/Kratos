/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
/* *********************************************************
*
*   Last Modified by:    $Author: Nelson Lafontaine $
*   Date:                $Date: 26-06-2009 $
*   Revision:            $Revision: 1.0 $
*
* ***********************************************************/


/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

//#include "custom_utilities/tensor_utils.h"
#include "includes/ublas_interface.h"
#include "includes/properties.h"
#include "soft_hard_behavior/friction_softening.h"
#include <cmath>
#include <string>
#include <iostream>


namespace Kratos
{
	    
		   Friction_Softening::Friction_Softening():SofteningHardeningCriteria() {}
		   Friction_Softening::~Friction_Softening(){}
		   double  Friction_Softening::FunctionBehavior(const Vector& Imput_Parameters)
		   {
		     /*
		     const double& he   = Imput_Parameters[0]; /// Longituf del elemento
		     const double& Ep   = Imput_Parameters[1]; /// Deformacion plastica efectiva
		     const double& Ft   = (*mprops)[FT];
	             const double& Ec   = (*mprops)[YOUNG_MODULUS];
	             const double& GE   = (*mprops)[FRACTURE_ENERGY];
	             const double Ef    = (2.00 * GE)/(Ft * he);
		     double result      = 0.00;
		     double frac        = 2.00 * std::sqrt(Ep * Ef ) /(Ep + Ef);
		     double sin         = 0.00;
		     double friction    =  PI * (*mprops)[INTERNAL_FRICTION_ANGLE] / 180.00;
		     if(Ep<Ef){
		        result =  std::asin(frac *std::sin(friction)); 
		        result =  180.00 * result/ PI;}
		     else
		        result = (*mprops)[INTERNAL_FRICTION_ANGLE]; 
		     */
		     const double fric = (*mprops)[INTERNAL_FRICTION_ANGLE];
		     double sinphi     =  std::sin(PI * (fric) / 180.00);
		     const double Kl   =  1.00;
		     const double& Kp  =  Imput_Parameters[3];
		     double frac       =  2.00 * std::sqrt(Kl * Kp ) /(Kl + Kp);
		     double result     =  0.00;
		     
		     return (*mprops)[INTERNAL_FRICTION_ANGLE];
		     
		     if(Kp<Kl)
		         result = std::asin(frac*sinphi) * 180/PI;
		      else
		        result   =  (*mprops)[INTERNAL_FRICTION_ANGLE];
		      
		      if(result!=result)
		      {   KRATOS_WATCH(Kp)
		          KRATOS_WATCH(Imput_Parameters ) 
			  KRATOS_WATCH(frac)
			  KRATOS_WATCH(result)
			  KRATOS_ERROR(std::logic_error,  "FRICTION" , ""); 
		      }
			
		      if(result<1.00)
			 result = 1.00;
		      
		     return result;
		   }
		   
		   double  Friction_Softening::FirstDerivateFunctionBehavior(const Vector& Imput_Parameters)
		   {
		     const double& he   = Imput_Parameters[0]; /// Longituf del elemento
//		     const double& Ep   = Imput_Parameters[1]; /// Deformacion plastica efectiva
		     const double& Ft   = (*mprops)[FT];
//	             const double& Ec   = (*mprops)[YOUNG_MODULUS];
	             const double& GE   = (*mprops)[FRACTURE_ENERGY];
//	             const double Ef    = (2.00 * GE)/(Ft * he); /// parametro a cambiar
		     double result      = 0.00;
//		     double frac        = 2.00 * std::sqrt(Ep * Ef ) /(Ep + Ef);
//		     double sin         = 0.00;
//		     double friction    =  PI * (*mprops)[INTERNAL_FRICTION_ANGLE] / 180.00;
		     /*
		        if(Ep<Ef){
		        double aux         =  2.00 * std::sin(friction) * std::sqrt(Ep * Ef ) /(Ep + Ef);
			const double aux_1 =  std::sqrt(Ef*Ep)*((Ep+Ef)*(Ep+Ef));
			aux                =  std::sqrt(1.00-aux *aux);
			if(aux_1>1E-6){
		           result =  -(Ep-Ef)*Ef*std::sin(friction)/aux_1;
			   result = 0.00; // result/aux;
			  }
 		        }
		     else
		        result = 0.00;
		     */
		     return result;
		   }


                     
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/


