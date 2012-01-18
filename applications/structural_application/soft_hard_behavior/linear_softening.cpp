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
#include "soft_hard_behavior/linear_softening.h"
#include <cmath>
#include <string>
#include <iostream>


namespace Kratos
{
	    
		   Linear_Softening::Linear_Softening():SofteningHardeningCriteria() {}
		   Linear_Softening::~Linear_Softening(){}
               

                  double  Linear_Softening::FunctionBehavior(const Vector& Imput_Parameters)
                  {
		     const double& he         = Imput_Parameters[0]; /// Longituf del elemento
		     const double& Ep         = Imput_Parameters[1]; /// Deformacion plastica efectiva
		     const double& Ft         = (*mprops)[FT];
//		     const double& Ec         = (*mprops)[YOUNG_MODULUS];
		     const double& GE         = (*mprops)[FRACTURE_ENERGY];
		     const double Hs          = Ft*Ft*he/(2.00 * GE);
		     const double fact        = 1.00 - Hs * Ep/Ft;
		     const double ft          = Ft * fact; 
		     //const double param     = 0.30 *(*mprops)[FT];
		     double result            = ft;  //> param ?  ft:param;    
                     return result;       
		  }
                  
                  double Linear_Softening::FirstDerivateFunctionBehavior(const Vector& Imput_Parameters)
                  {     
		     const double& he         = Imput_Parameters[0]; /// Longituf del elemento
//		     const double& Ep         = Imput_Parameters[1]; /// Deformacion plastica efectiva
		     const double& Ft         = (*mprops)[FT];
//		     const double& Ec         = (*mprops)[YOUNG_MODULUS];
		     const double& GE         = (*mprops)[FRACTURE_ENERGY];
		     const double Hs          = Ft*Ft*he/(2.00 * GE);
		     //if(ft>param)
		     double result            = -Hs;
		     return result; 
		  }

                   // for geomaterials model
                   void Linear_Softening::FunctionSofteningHardeningBehavior(const double& capap, const double& sigma, double& Result, double& der_Result)

                     {
                         //Result     = 0.00;
                         //der_Result = -sigma;
			 if (capap<=1.00)
                               {
                                  Result     = sigma*(1.00 - capap);
                                  der_Result = -sigma;
                               }
                      //KRATOS_WATCH(Result)
                      //KRATOS_WATCH(der_Result)     
                      return;
                     }
                     
                  double Linear_Softening::Calculate( const Vector& Imput_Parameters)
		   {
		     const double& Ft   = (*mprops)[FT];
		     const double& Ec   = (*mprops)[YOUNG_MODULUS];
		     const double& GE   = (*mprops)[FRACTURE_ENERGY];
		     const double& l    = Imput_Parameters[0];
		     const double& r_o  = Imput_Parameters[1];  // Ft/sqrt(Ec);
		     const double& r    = Imput_Parameters[2];
		     
		     const double Hs_barra = Ft * Ft /(  2.00 * Ec * GE );
		     const double ls_barra = 1.00 / Hs_barra; 
 		     double Hs       = l / (ls_barra - l);
		     if (Hs < 0.00)
 		      {    
 		        Hs  = 0.00;
                        std::cout<<"Warning: Softening Parameters is less than zero. Please refine more."<<std::endl;
 		      }
 		     
 		     const double ru =  r_o * (1.00  + 1.00/Hs); 
		     double q_r      = 0.00;
		     if ( r>=r_o  && r <= ru)
		     {
		       q_r = r_o - Hs * (r - r_o );
		     }
		     double a    =  1.00 - (q_r/r);
		     return a;
		    }
                     
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/


