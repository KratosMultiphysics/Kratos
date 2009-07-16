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

#if !defined( KRATOS_FLUENCY_CRITERIA)
#define KRATOS_FLUENCY_CRITERIA

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
//#include "includes/define.h"
//#include "includes/variables.h"
//#include "utilities/math_utils.h"
//#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
#include "includes/ublas_interface.h"
#include <cmath>


namespace Kratos
{

     class FluencyCriteria
    {
        public:


		    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector		  
		    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;

		    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;

		    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz.

		    typedef FluencyCriteria FluencyCriteriaType; 

		    FluencyCriteria(){}

		    virtual ~FluencyCriteria(){}

		    virtual void InitializeMaterial() {}


		    virtual void CalculateEquivalentUniaxialStress(const Vector& StrainVector, 
		    const Vector& StressVector,
		    const Matrix& ConstitutiveMatrix,
		    const double& ro, const double& n, 
		    double& Result) 
		    {}

		    virtual void CalculateEquivalentUniaxialStressViaPrincipalStress(const Vector& StrainVector, 
		    const Vector& StressVector,double& Result)
		    {}


		    virtual void CalculateEquivalentUniaxialStressViaInvariants(const Vector& StrainVector, 
		    const Vector& StressVector,double& Result)

		    {}

		    virtual void CalculateEquivalentUniaxialStressViaCilindricalCoordinate(const Vector& StrainVector, 
		    const Vector& StressVector,double& Result)

		    {}


		    virtual void CalculateDerivateFluencyCriteria(Vector DerivateFluencyCriteria)
		    {}

		    protected:

		    static void  Comprobate_State_Tensor(Matrix& StressTensor, const Vector& StressVector)
		    {
		  // Necesario para calcular eigen valores con subrutina de Jacobi, NO ACEPTA TERMINOS NULOS. 
		  if (fabs(StressTensor(0,0))<1E-10){StressTensor(0,0) = 1E-10; }
		  if (fabs(StressTensor(0,1))<1E-10){StressTensor(0,1) = 1E-10; }   
		  if (fabs(StressTensor(1,0))<1E-10){StressTensor(1,0) = 1E-10; }
		  if (fabs(StressTensor(1,1))<1E-10){StressTensor(1,1) = 1E-10; }

		  if (StressVector.size()== 6) 
		    {
		  if (fabs(StressTensor(0,2))<1E-10){StressTensor(0,2) = 1E-10; }
		  if (fabs(StressTensor(1,2))<1E-10){StressTensor(1,2) = 1E-10; }
		  if (fabs(StressTensor(2,0))<1E-10){StressTensor(2,0) = 1E-10; }
		  if (fabs(StressTensor(2,1))<1E-10){StressTensor(2,1) = 1E-10; }
		  if (fabs(StressTensor(2,2))<1E-10){StressTensor(2,2) = 1E-10; }  
		    }

		  }

		    private:



    }; /* Class FluencyCriteria */
    
    /**
     * definition of CONSTITUTIVE_LAW variable
     */
}  /* namespace Kratos.*/
#endif /* FLUENCY_CRITERIA defined */
