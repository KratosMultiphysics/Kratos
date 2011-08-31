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
#include "soft_hard_behavior/cohesion_softening.h"
#include <cmath>
#include <string>
#include <iostream>


namespace Kratos
{
	    
		   Cohesion_Softening::Cohesion_Softening():SofteningHardeningCriteria() {}
		   Cohesion_Softening::~Cohesion_Softening(){}
		   double  Cohesion_Softening::FunctionBehavior(Vector& Imput_Parameters)
		   {
		     const double PI    = 3.1415926535898; 
		     const double& he   = Imput_Parameters[0]; /// Longituf del elemento
		     const double& Ep   = Imput_Parameters[1]; /// Deformacion plastica efectiva
		     const double& Ft   = (*mprops)[FT];
		     const double& Ec   = (*mprops)[YOUNG_MODULUS];
		     const double& GE   = (*mprops)[FRACTURE_ENERGY];
		     const double Hs    = (Ft*Ft*Imput_Parameters[0]) / (2.00 * GE);
		     const double elev  = -2.00 * Hs * Imput_Parameters[1]/ Ft; 
		     const double fric  = (*mprops)[INTERNAL_FRICTION_ANGLE] * PI/180.00;
                     const double N_phi =  tan(PI/4.00 + fric/2.00); 
		     const double one   =  Ft * exp(elev);
		     const double alfa  = (*mprops)[FC] / (*mprops)[FT];
		     
		     double cohe   = one  * alfa / (N_phi *2.00);
		     double param  = 0.30 *(*mprops)[COHESION];
		     double result = cohe > param ?  cohe:param;    
                     
		     
		     return result;  
		     
		     /*
		     //Imput_Parameters[0] = effectivePlasticStrain
		     double result = 0;
		     const double& cohe = (*mprops)[COHESION];  
  		     if(Imput_Parameters[0]<0.01)
  		         result =   -cohe/0.01 * Imput_Parameters[0] + cohe;
  		     return result;  
		     */
		   }
		   
		   double  Cohesion_Softening::FirstDerivateFunctionBehavior(Vector& Imput_Parameters)
		   {
		     double result      = 0.00;
// 		     const double PI    = 3.1415926535898; 
// 		     const double& he   = Imput_Parameters[0]; /// Longituf del elemento
// 		     const double& Ep   = Imput_Parameters[1]; /// Deformacion plastica efectiva
// 		     const double& Ft   = (*mprops)[FT];
// 		     const double& Ec   = (*mprops)[YOUNG_MODULUS];
// 		     const double& GE   = (*mprops)[FRACTURE_ENERGY];
// 		     const double Hs    = (Ft*Ft*Imput_Parameters[0]) / (2.00 * GE);
// 		     const double elev  = -2.00 * Hs * Imput_Parameters[1]/ Ft; 
// 		     const double fric  = (*mprops)[INTERNAL_FRICTION_ANGLE] * PI/180.00;
// 		     const double N_phi =  tan(PI/4.00 + fric/2.00); 
// 		     const double one   =  Ft * exp(elev);
// 		     const double alfa  = (*mprops)[FC] / (*mprops)[FT];
// 
// 		    result =  (-2.00 * Hs / Ft)  * one * alfa / (N_phi *2.00);  
		     
		     
		     
//  		     const double& cohe = (*mprops)[COHESION];  
//   		     if(Imput_Parameters[0]<0.01)
//   		         result =   -cohe/0.01;
  		     return result; 
		   }


                     
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/


