/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
/* *********************************************************   
*          
*   Last Modified by:    $Author: Nelson Lafontaine //inglafontaine@gmail.com    $ 
*   Date:                $Date: 2009-7-13$
*   Revision:            $Revision: 1.00   $
*
* ***********************************************************/

#if !defined(MAPED_SPACE)
#define MAPED_SPACE


#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "custom_utilities/tensor_utils.h"
//#include "custom_strategies/strategies/residualbased_newton_raphson_line_search_strategy.h"
#include <cmath>

namespace Kratos
{
    template<class TDataType> class Maped_Space
    {
        public:
            /** 
             * @name type definitions
             * @{
             */
            //typedef Matrix MatrixType;
		
            //typedef Vector VectorType;
		
            typedef unsigned int IndexType;
            
            typedef unsigned int SizeType;

	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector
			  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz.

            typedef MathUtils<TDataType> MathUtilsType;
	    
	    //typedef typename ResidualBasedNewtonRaphsonLineSearchesStrategy  ResidualBasedNewtonRaphsonLineSearchesStrategy; 


           static inline void Calculate_Omega_Tensor(const Vector& Orthotropic_Elastic_Limit, const double& Isotropic_Elastic_Limit, Matrix& W_Matrix)
	    {
	     
             //Orthotropic_Limit_Elastic = [fxx, fyy, fz, fxy, fyz, fxz] 3D // umbrales de resistencia limites en direccion de ortotropia
             //Orthotropic_Limit_Elastic = [fxx, fyy, fxy] 2D
	      unsigned int dim  = Orthotropic_Elastic_Limit.size();
	      unsigned int size = 3;
	      if (dim==3){size = 2;}
	     
	      Second_Order_Tensor omega(size);
	      Fourth_Order_Tensor W_Tensor(size);
	      if (size==2)
		{

                    
		    W_Matrix.resize(3,3,false);	  

		    omega[0].resize(2,false);
                    omega[1].resize(2,false);

		    W_Tensor[0].resize(2, false);
                    W_Tensor[1].resize(2, false);
		    W_Tensor[0][0].resize(2,2, false); W_Tensor[0][1].resize(2,2, false);
		    W_Tensor[1][0].resize(2,2, false); W_Tensor[1][1].resize(2,2, false);

		    omega[0][0] = sqrt(Orthotropic_Elastic_Limit(0)/Isotropic_Elastic_Limit);
 		    omega[0][1] = sqrt(Orthotropic_Elastic_Limit(2)/Isotropic_Elastic_Limit); 
 		    omega[1][0] = omega[0][1];
 		    omega[1][1] = sqrt(Orthotropic_Elastic_Limit(1)/Isotropic_Elastic_Limit); 
 		   
		    Tensor_Utils<double>::Prod_Second_Order_Tensor(omega, omega, W_Tensor);
                    SD_MathUtils<double>::TensorToMatrix(W_Tensor, W_Matrix);

		}
    
	    }

 /*          static inline void Calculate_Tensor_Ajuste()
	    {
	    }
*/	    
                   
        private:

	    //Second_Order_Tensor momega;
	    //Fourth_Order_Tensor W_Tensor;
	    //Matrix W_Matrix;

    
	    //ResidualBasedNewtonRaphsonLineSearchesStrategy mSolver;  
	    

    };// class SD_MathUtils
}
#endif /* SD_MATH_UTILS defined */
