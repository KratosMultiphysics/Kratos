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
		   double  Cohesion_Softening::FunctionBehavior(const Vector& Imput_Parameters)
		   {
		     const double c     = (*mprops)[COHESION];
//		     const double& Ft   = (*mprops)[FT];
//		     const double& Fc   = (*mprops)[FC];
//		     const double Ro    = Fc/Ft;
		     const double& he   = Imput_Parameters[0];       /// Longituf del elemento
		     const double& Ep   = Imput_Parameters[1];       /// Deformacion plastica efectiva
//	             const double& Em   = (*mprops)[YOUNG_MODULUS];
	             const double& GE   = (*mprops)[FRACTURE_ENERGY];
	             const double Ec    = 2.00 * GE /(std::sqrt(PI) * c * he); 
                     const double elev  = (Ep/Ec)*(Ep/Ec); 
//		     double cohe        = (*mprops)[COHESION]*std::exp(-elev);
		     double result      = c;  //cohe;   //> param ?  cohe:param;    
		     return result;          //= (*mprops)[COHESION];
		   }
		   
		   double  Cohesion_Softening::FirstDerivateFunctionBehavior(const Vector& Imput_Parameters)
		   {
		     
		     const double c     = (*mprops)[COHESION];
//		     const double& Ft   = (*mprops)[FT];
//		     const double& Fc   = (*mprops)[FC];
//		     const double Ro    = Fc/Ft;
		     const double& he   = Imput_Parameters[0];       /// Longituf del elemento
//		     const double& Ep   = Imput_Parameters[1];       /// Deformacion plastica efectiva
//	             const double& Em   = (*mprops)[YOUNG_MODULUS];
	             const double& GE   = (*mprops)[FRACTURE_ENERGY];
	             const double Ec    = 2.00 * GE /(std::sqrt(PI) * c * he); 
//                     const double elev  = (Ep/Ec)*(Ep/Ec); 
		     
		     double cohe   =  0.00;   //2.00 * (Ep/Ec) * (*mprops)[COHESION]*std::exp(-elev);
		     double result =  cohe;   //> param ?  cohe:param;    
		     return result;           //= (*mprops)[COHESION];
		   }

                   //[0]->kptn1
		   //[1]->kpn
                   double Cohesion_Softening::EvolucionLaws(const Vector& Imput_Parameters, const array_1d<double,3>& Sigma)
                   {
		     
		      double Kpn1  =  Imput_Parameters[0];
		      double Kpn   =  Imput_Parameters[1];
		      double cn    =  Imput_Parameters[2];
		      double co    =  (*mprops)[COHESION];
		      double result      =  0.00;
//		      const double toler =  1E-3;
		      
		      if(Kpn1<1.00){
		      
			
		      result             = co;
	              double teta_a      =  Tensor_Utils<double>::Mc_aully(Sigma);
		      double teta_b      =  std::fabs(Sigma[0]) + std::fabs(Sigma[1]) + std::fabs(Sigma[2]);
		      double teta        =  0.00;        
		      if (fabs(teta_b) < 1E-10)
		        {teta = 0.00;}
		      else
		         {teta = teta_a/teta_b;}
		         
		     double fact1 = 0.00; 
		     double fact2 = 0.00; 
		     double value    =  0.00; 
		     double derivate =  0.00;
                     FuntionCCAndDerivate(Kpn1, co, value, derivate);   
		     if(value!=0.00)
		         fact1 = teta * derivate / value;
		     FuntionCCAndDerivate(Kpn1, co, value, derivate);
		     if(value!=0.00)
		         fact2  = (1.00 - teta) * derivate / value;
		     
		     double fact   = (1.00 - (fact1 + fact2)*(Kpn1-Kpn));
		     //double fact   = (1.00 + (fact1 + fact2)*(Kpn1-Kpn)); 
		     if(std::fabs(fact)!=0.00)
		         result = cn / fact;
		     
		     if(result>co)
		        KRATOS_ERROR(std::logic_error,  "ERROR COHESION" , "");
		     
		      }
		      else
			result = 0.00;
		      
		      return  result;
		     }
                     
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/


