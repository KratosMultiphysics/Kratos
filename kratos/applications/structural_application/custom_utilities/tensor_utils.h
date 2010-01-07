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

static inline double double_product(const Matrix &A,const Matrix &B)  //A:B = Tr(A^tB)
      {
      KRATOS_TRY
      
      unsigned int dim = A.size1();
      double result       = 0.00;
      Matrix A_Trans      = ZeroMatrix(dim,dim);
      Matrix Aux          = ZeroMatrix(dim,dim);
      noalias(A_Trans)    = trans(A);
      noalias(Aux)        = prod(trans(A_Trans),B);
      result              = Trace(Aux);
      return result; 
      KRATOS_CATCH("")
      } 


//***********************************************************************
//*********************************************************************** 

// Calculo de los invariantes 
// I: Invariante I1, I2, I3
// J: Invariantes J;
// Jdes: Invariantes del tensor Desviador
// Ver Capitulo de Mecanica de medios Continuos para ingenieros de Oliver. pag 233. 

//NOTA: Los invariantes se calculan para un estado tridimensional.

static inline void TensorialInvariants(const Matrix& Tensor, Vector& I, Vector& J, Vector& J_des)

{
	KRATOS_TRY
	//int iter                            = 50;
	//double zero                         = 1.0E-15;
        int dim                             = Tensor.size1();
	matrix<double> Aux_Tensor           = ZeroMatrix(3,3);
	matrix<double> SphericComponent     = IdentityMatrix(3,3);
	matrix<double> DesviatoricComponent = ZeroMatrix(3,3);
	vector<double> PrincipalStress      = ZeroVector(3);
        matrix<double> EigenVectors         = ZeroMatrix(3,3);  

	// Los invariantes seran representados como vectores
	I     = zero_vector<double>(3);
	J     = zero_vector<double>(3);
	J_des = zero_vector<double>(3);
 
        if(dim==3)
	    {
              Aux_Tensor = Tensor;
	    } 
        else if(dim==2)
            {
             Aux_Tensor(0,0) = Tensor(0,0); Aux_Tensor(0,1) = Tensor(0,1);
             Aux_Tensor(1,0) = Tensor(1,0); Aux_Tensor(1,1) = Tensor(1,1);
             //Comprobate_State_Tensor(Aux_Tensor);
             }
    
	//SD_MathUtils<double>::EigenVectors(Tensor, EigenVectors, PrincipalStress, zero, iter); 
//      I[0] = PrincipalStress(0) + PrincipalStress(1) + PrincipalStress(2);
// 	I[1] = (PrincipalStress(0)*PrincipalStress(1) + PrincipalStress(0)*PrincipalStress(2) + PrincipalStress(1)*PrincipalStress(2)); //Pag 39 javier Bonet
// 	I[2] = MathUtils<double>::Det(Tensor);

	I[0] = Aux_Tensor(0,0) + Aux_Tensor(1,1) + Aux_Tensor(2,2);
	I[1] = 0.50*(double_product(Aux_Tensor,Aux_Tensor) - I[0]*I[0]); 
	I[2] = MathUtils<double>::Det(Aux_Tensor);
	
  
	

	// Invariantes J
	J[0] =  I[0];
	J[1] =  0.50*(I[0]*I[0]+2.00*I[1]);
	J[2] =  (I[0]*I[0]*I[0] + 3.00*I[0]*I[1]+3.00*I[2])/3.00;

	noalias(SphericComponent)     =  (I(0)/3.00)*SphericComponent;
	noalias(DesviatoricComponent) =  Aux_Tensor - SphericComponent;

	J_des[0] = 0.00;
	J_des[1] = 0.50*double_product(DesviatoricComponent,DesviatoricComponent);
	J_des[2] = MathUtils<double>::Det(DesviatoricComponent);
    
	//KRATOS_WATCH(Tensor)
        //KRATOS_WATCH(DesviatoricComponent)
        //KRATOS_WATCH(I)
        //KRATOS_WATCH(J)
        //KRATOS_WATCH(J_des)


	KRATOS_CATCH("")
}

//***********************************************************************
//*********************************************************************** 

static inline void  Comprobate_State_Tensor(Matrix& StressTensor)
    { 
  if (fabs(StressTensor(0,0))<1E-15){StressTensor(0,0) = 1E-15; }
  if (fabs(StressTensor(1,1))<1E-15){StressTensor(1,1) = 1E-15; }
  if (fabs(StressTensor(2,2))<1E-15){StressTensor(2,2) = 1E-15; }  
    }


//***********************************************************************
//*********************************************************************** 

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
	      for(unsigned int k=0;k<size; k++){
                  for(unsigned int l=0;l<size; l++){
		           Result[i][j](k,l) = A[i](j)*B[k](l);
                                    }
                                 } 
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

static inline void Identy_Fourth_Order_Tensor(unsigned int dim, Fourth_Order_Tensor& Result)
{
      
    if (dim==3)
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

      for(unsigned int i=0;i<dim; i++ ){
	  for(unsigned int j=0;j<dim; j++){
	      for(unsigned int k=0;k<dim; k++){
                  for(unsigned int l=0;l<dim; l++){
                               
			      if((i==k) && (j==l))
		                  {
                                     Result[i][j](k,l) = 1.00;
                                  }
                              else {
                                     Result[i][j](k,l) = 0.00;
                                    }
                                 } 
                             }
                         }     
                     }

}

static inline void Contraction_Double(Fourth_Order_Tensor& Fourth_Tensor, Second_Order_Tensor Second_Tensor, Second_Order_Tensor Result_Tensor)
{ 
  unsigned int size = Second_Tensor[0].size(); 
  for (unsigned int i=0;i<size; i++){
     for(unsigned int j=0; j<size; j++){
        for (unsigned int k=0; k<size; k++){
            for (unsigned int l =0; j<size; l++){ 
                Result_Tensor[i](j)= Result_Tensor[i](j) + Fourth_Tensor[i][j](k,l)*Second_Tensor[k](l);}
                          }
                }
         } 
}

static inline void Contraction_Double(Fourth_Order_Tensor& Fourth_Tensor_A, Fourth_Order_Tensor& Fourth_Tensor_B , Fourth_Order_Tensor& Result)
{ 
    unsigned int size = Fourth_Tensor_B[0].size(); 

    
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


  for (unsigned int i=0;i<size; i++){
     for(unsigned int j=0; j<size; j++){
        for (unsigned int k=0; k<size; k++){
            for (unsigned int l =0; l<size; l++){
               for (unsigned int m=0; m<size; m++){
                  for (unsigned int n =0; n<size; n++){
                       Result[i][j](k,l) = Result[i][j](k,l) + Fourth_Tensor_A[i][j](m,n)*Fourth_Tensor_B[m][n](k,l);}
                          }
                      }
                } 
           }
      }
}


    private:

    };//tensor_utils
}
#endif /* TENSOR_UTILS defined */
