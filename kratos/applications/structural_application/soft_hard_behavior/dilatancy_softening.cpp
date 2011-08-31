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
		   double  Dilatancy_Softening::FunctionBehavior(Vector& Imput_Parameters)
		   {
		     //Imput_Parameters[0] = effectivePlasticStrain
		     //Imput_Parameters[1] = currentcohesion
		     double result = 0;
		     result =   (*mprops)[DILATANCY_ANGLE];
// 		     if(Imput_Parameters[0]<0.1)
// 		         result =  -2000 * Imput_Parameters[0] +  (*mprops)[COHESION];
// 		     else
// 		       result = 0;
		     return result;
		       
		   }
		   double  Dilatancy_Softening::FirstDerivateFunctionBehavior(Vector& Imput_Parameters)
		   {
		     double result = 0; 
// 		     if(Imput_Parameters[0]<0.1)
// 		         result =  -2000;
// 		     else
// 		       result = 0;
		     
		     return result; 
		   }


                     
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/


