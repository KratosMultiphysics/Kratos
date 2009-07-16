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
*/

#if !defined(TENSOR_UTILS)
#define TENSOR_UTILS

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include <cmath>


namespace Kratos
{
    template<class TDataType> class Tensor_Utils
    {
        public:
            /** 
             * @name type definitions
             * @{
             */
            typedef Matrix MatrixType;
		
            typedef Vector VectorType;
		
            typedef unsigned int IndexType;
            
            typedef unsigned int SizeType;
            
            typedef  Tensor_Utils<TDataType> Tensor_Utils_Type;

	    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor; // dos opciones: un tensor de segundo orden y/o un vector que almacena un vector
			  
	    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;
			  
            typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;
			  
	    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor; // Acumulo un tensor de 2 orden en una matriz.


            
//***********************************************************************
//***********************************************************************

static inline double Sign(const double &Def)
	{	

	KRATOS_TRY
	if (Def>0.00)
	{
	return 1.00;
	} 
	else if(Def<0)
	{
	return -1.00;
	}
	else
	{
	return  0.00;
	}
	KRATOS_CATCH("")
	}


 //***********************************************************************
//*********************************************************************** 
 
static inline double Mc_aully(const vector<double> & v)

      {
      KRATOS_TRY
      double Result = 0.00;
      for(unsigned int i=0; i<v.size(); ++i)
      {
      Result += fabs(v[i])+v[i];
      }
      return 0.50*Result;
      KRATOS_CATCH("")
      }

//***********************************************************************
//*********************************************************************** 

static inline Vector Sign(const vector<double> &v)

      {	
      KRATOS_TRY
      Vector Result (v.size());
      for(unsigned int i=0; i<v.size(); ++i)
      {
      Result[i] = Sign(v[i]); 
      }
      return Result; 
      KRATOS_CATCH("")
      }

//***********************************************************************
//***********************************************************************

static inline double Trace(const Matrix &A)

      {
      KRATOS_TRY
      double Result = 0.00;
      for(unsigned int i=0; i<A.size1(); ++i) 
      {
      Result += A(i,i); 		
      }
      return Result;
      KRATOS_CATCH("")
      } 


 
//***********************************************************************
//***********************************************************************

static inline void SphericandDesviatoricTensor(const Vector & StressVector, Matrix &SphericComponent, Matrix &DesviatoricComponent)

	{
	  KRATOS_TRY
	  double crit      = 1.0E-15;
	  double zero      = 1.0E-15;	
	  double Sigma_m   = 0.00;          // Tension Media
	  unsigned int dim = StressVector.size();            
	  matrix<double> Tensor(0,0);
	  vector<double> PrincipalStress(0);

	  if(dim==3)
	  {
	  Tensor= zero_matrix<double>(2,2);
	  SphericComponent = zero_matrix<double>(2,2);    
	  DesviatoricComponent = zero_matrix<double>(2,2);
	  PrincipalStress = zero_vector<double>(2);  
	  }
	  else
	  {
	  Tensor= zero_matrix<double>(3,3);
	  SphericComponent = zero_matrix<double>(3,3);    
	  DesviatoricComponent= zero_matrix<double>(3,3);
	  PrincipalStress = zero_vector<double>(3);      
	  }


	  Tensor  = MathUtils<double>::StressVectorToTensor(StressVector);
	  PrincipalStress  = SD_MathUtils<double>::EigenValues(Tensor,crit, zero);

	  if(dim==3)			    				
	  Sigma_m = (sum(PrincipalStress))/2.00;
	  else
	  Sigma_m = (sum(PrincipalStress))/3.00;

	  for(unsigned int i=0; i<SphericComponent.size1(); ++i)
	  {
	  SphericComponent(i,i) = Sigma_m;
	  }

	  noalias(DesviatoricComponent) = Tensor - SphericComponent;
	  KRATOS_CATCH("") 

	  }

//***********************************************************************
//*********************************************************************** 

// Calculo de los invariantes 
// I: Invariante I1, I2, I3
// J: Invariantes J;
// Jdes: Invariantes del tensor Desviador
// Ver Capitulo de Mecanica de medios Continuos para ingenieros de Oliver. pag 233. 

static inline void TensorialInvariants(const Vector & StressVector, Vector &I, Vector &J, Vector &J_des)

{
	KRATOS_TRY
	double crit      = 1.0E-15;
	double zero      = 1.0E-15;
	unsigned int dim  = StressVector.size();
	unsigned int size = 3;
	matrix<double> Tensor(0,0);
	matrix<double> Aux_Tensor(0,0);
	matrix<double> Aux_Matrix(0,0);
	matrix<double> SphericComponent(0,0);
	matrix<double> DesviatoricComponent(0,0);
	vector<double> PrincipalStress(0);  

	if(dim==3)
	{  
	size       = 2; 
	Tensor     = zero_matrix<double>(2,2);
	Aux_Tensor = zero_matrix<double>(2,2);    
	Aux_Matrix = zero_matrix<double>(2,2);
	PrincipalStress = zero_vector<double>(2);
	}
	else
	{
	Tensor     = zero_matrix<double>(3,3);
	Aux_Tensor = zero_matrix<double>(3,3);    
	Aux_Matrix = zero_matrix<double>(3,3);
	PrincipalStress = zero_vector<double>(3);
	}

	// Los invariantes seran representados como vectores

	I     = zero_vector<double>(3);
	J     = zero_vector<double>(3);
	J_des = zero_vector<double>(3);


	Tensor  = MathUtils<double>::StressVectorToTensor(StressVector);
	PrincipalStress  = SD_MathUtils<double>::EigenValues(Tensor,crit, zero);
	for(unsigned int i = 0; i<PrincipalStress.size(); ++i)
	{
		      Aux_Tensor(i,i) = PrincipalStress(i);
	}

	// Invariantes I	
	I[0] = Trace(Aux_Tensor);
	noalias(Aux_Matrix) = prod(trans(Aux_Tensor),Aux_Tensor);
	I[1] = 0.5*(Trace(Aux_Matrix)-I[0]*I[0]);
	I[2] = MathUtils<double>::Det(Tensor);

	// Invariantes J
	J[0] =  I[0];
	J[1] =  0.50*(I[0]*I[0]+2.00*I[1]);
	J[2] =  (I[0]*I[0]*I[0] + 3.00*I[0]*I[1]+3.00*I[2])/3.00;


	SphericandDesviatoricTensor(StressVector, SphericComponent, DesviatoricComponent);

	PrincipalStress  = SD_MathUtils<double>::EigenValues(DesviatoricComponent,crit, zero);
	for(unsigned int i = 0; i<PrincipalStress.size(); ++i)
	{
		      Aux_Tensor(i,i) = PrincipalStress(i);
	}

	J_des[0] = Trace(Aux_Tensor);
	noalias(Aux_Matrix) = prod(trans(Aux_Tensor),Aux_Tensor);
	J_des[1] = 0.5*(Trace(Aux_Matrix));
	J_des[2] = MathUtils<double>::Det(Aux_Tensor);

	KRATOS_CATCH("")
}

static inline void Prod_Second_Order_Tensor(const Second_Order_Tensor& A,const Second_Order_Tensor& B, Fourth_Order_Tensor& Result)
{
      
      
    unsigned int size = A[0].size();  
    
    if (size==3)
    {
    Result.resize(3);
    
    Result[0].resize(3);
    Result[1].resize(3);
    Result[2].resize(3);
    
    Result[0][0].resize(3,3, false); Result[0][1].resize(3,3, false); Result[0][2].resize(3,3, false);
    Result[1][0].resize(3,3, false); Result[1][1].resize(3,3, false); Result[1][2].resize(3,3, false);
    Result[2][0].resize(3,3, false); Result[2][1].resize(3,3, false); Result[2][2].resize(3,3, false);
    }
    else
    {
    Result.resize(2);
    
    Result[0].resize(2);
    Result[1].resize(2);
    
    Result[0][0].resize(2,2, false); Result[0][1].resize(2,2, false); 
    Result[1][0].resize(2,2, false); Result[1][1].resize(2,2, false); 
    }

      for(unsigned int i=0;i<size; i++ ){
	  for(unsigned int j=0;j<size; j++){
		  Result[i][j] = outer_prod(A[i], B[j]);
		                                        }
	          }
}


static inline void Sum_Fourth_Order_Tensor(const Fourth_Order_Tensor& A,const Fourth_Order_Tensor& B, Fourth_Order_Tensor& Result)
{
      unsigned int size = A[0].size();	    
      for (unsigned int i=0; i<size;i++){
	  for (unsigned int j=0; j<size;j++){
		noalias(Result[i][j]) = A[i][j] + B[i][j];}
		  } 
}

static inline void Rest_Fourth_Order_Tensor(const Fourth_Order_Tensor& A,const Fourth_Order_Tensor& B, Fourth_Order_Tensor& Result)
{
      unsigned int size = A[0].size();	    
      for (unsigned int i=0; i<size;i++){
	  for (unsigned int j=0; j<size;j++){
		noalias(Result[i][j]) = A[i][j] - B[i][j];}
		  } 
}


    private:

    };//tensor_utils
}
#endif /* TENSOR_UTILS defined */
