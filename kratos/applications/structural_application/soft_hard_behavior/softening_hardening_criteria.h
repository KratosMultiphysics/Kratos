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

#if !defined( KRATOS_SOFTENING_HARDENING_CRITERIA)
#define KRATOS_SOFTENING_HARDENING_CRITERIA

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

#include "custom_utilities/tensor_utils.h"
#include "includes/ublas_interface.h"
#include "includes/properties.h"
#include <cmath>
#include <string>
#include <iostream>



namespace Kratos
{


      class SofteningHardeningCriteria
	  {
	    public:
           
            SofteningHardeningCriteria(){}
            virtual ~SofteningHardeningCriteria(){}
            KRATOS_CLASS_POINTER_DEFINITION( SofteningHardeningCriteria );
	    
	    virtual double FunctionSofteningHardeningBehavior(const double& A, const double& r_o, const double& r) 
	      { 
		    return 0.00;
	      }
         
            // Calcula la funcion de abalandamiento o endurecimiento junto con su derivada respecto kp_punto (0 <= kp_punto <= 1.00) 
            virtual void FunctionSofteningHardeningBehavior(const double& capap, const double& sigma, double& Result, double& der_Result) 
	      { 
		    return;
	      }

            virtual void Linear_Strain_Softening(Vector& Principal_Stress)
             {
               return;
             } 
  
	    protected:

	    private:
                 
          };
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/
#endif /* FLUENCY_CRITERIA defined */

