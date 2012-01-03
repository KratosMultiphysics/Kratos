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

#if !defined( KRATOS_COHESION_SOFTENING_CRITERIA)
#define KRATOS_COHESION_SOFTENING_CRITERIA

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "soft_hard_behavior/softening_hardening_criteria.h"



namespace Kratos
{

	  class Cohesion_Softening: public SofteningHardeningCriteria
	    {
	    
		public: 
		   Cohesion_Softening();
		  ~Cohesion_Softening();
                   KRATOS_CLASS_POINTER_DEFINITION(Cohesion_Softening);
		   virtual boost::shared_ptr<SofteningHardeningCriteria> Clone() const
                    {
                      boost::shared_ptr<SofteningHardeningCriteria> p_clone(new Cohesion_Softening());
                      return p_clone;
                    }
		   
		   double   FunctionBehavior(const Vector& Imput_Parameters);
		   double   FirstDerivateFunctionBehavior(const  Vector& Imput_Parameters);
		   double   EvolucionLaws(const Vector& Imput_Parameters,const array_1d<double,3>& Sigma);
		   
	      private:
		   /// Esta funcion son con respecto a KPunto
		   ///  exponencial para la chesion en tracción, para kp es lineal
		   /// Sigma = Cohesion inicial
		   void FuntionCCAndDerivate(const double& Kp, const double& Sigma, double& value, double& derivate)
		   {
		     value    = 0.00;
		     derivate = 0.00;
		     if(Kp<1.00){
		        value    = (1.00 - Kp) * Sigma; 
			derivate = -Sigma;
		     }
		   }
		   
		   void FuntionCTAndDerivate(const double& Kp, const double& Sigma, double& value, double& derivate)
		   {
		     value    = 0.00;
		     derivate = 0.00;
		     if(Kp<1.00){
		        const double raiz = std::sqrt((1.00 - Kp));
		        value             = raiz * Sigma; 
			derivate          = -Sigma/(2.00 * raiz);
		     }
		     
		   }
		   
		   
		   
           };    
    

    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/
#endif /* FLUENCY_CRITERIA defined */

