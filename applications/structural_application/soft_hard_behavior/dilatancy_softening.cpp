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
#include "soft_hard_behavior/dilatancy_softening.h"
#include <cmath>
#include <string>
#include <iostream>


namespace Kratos
{
	    
		   Dilatancy_Softening::Dilatancy_Softening():SofteningHardeningCriteria() {}
		   Dilatancy_Softening::~Dilatancy_Softening(){}
		   double  Dilatancy_Softening::FunctionBehavior(const Vector& Imput_Parameters)
		   {
		     //const double PI    = 3.1415926535898; 
		     const double& he    =  Imput_Parameters[0];   /// Longituf del elemento
		     const double& Ep    =  Imput_Parameters[1];   /// Deformacion plastica efectiva
		     double angle        =  PI * Imput_Parameters[2]/180.00;
		     double result       =  0.00;
		     double friction     =  PI * (*mprops)[INTERNAL_FRICTION_ANGLE] / 180.00;
		     double dilatancy    =  PI * (*mprops)[DILATANCY_ANGLE] / 180.00;
		     double sinTv        =  (sin(friction)-sin(dilatancy))/(1.00 - (sin(friction)*sin(dilatancy)));
		     double sinangle     =  sin(angle);
		     result              =  (sinangle-sinTv)/(1.00-sinangle*sinTv); 
		     result              =  std::asin(result);
		     
		     return  (*mprops)[DILATANCY_ANGLE];
		     
		     if(angle< std::asin(sinTv))
                         result =  0.00;  
		     else
		         result =  180.00 * result/PI;
		     
		     if(result!=result)
			  KRATOS_ERROR(std::logic_error,  "DILATANCY" , "");
		     
		     if(result > (*mprops)[DILATANCY_ANGLE])
		        result= (*mprops)[DILATANCY_ANGLE];
		     
		     
		     if(result < 1.00)
		        result = 1.00; /// Initial value no must be 0: Cercano a cero es mas razonable 
		     
		     
		     return result; 
		     
		   }
		   
		   double  Dilatancy_Softening::FirstDerivateFunctionBehavior(const Vector& Imput_Parameters)
		   {
		     double result = 0; 
		     return result; 
		   }


                     
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/


