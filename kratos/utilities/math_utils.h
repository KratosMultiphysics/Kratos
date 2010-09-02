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
*   Last Modified by:    $Author: anonymous $
*   Date:                $Date: 2007-12-13 15:29:11 $
*   Revision:            $Revision: 1.8 $
*
* ***********************************************************/


#if !defined(KRATOS_MATH_UTILS )
#define  KRATOS_MATH_UTILS


/* System includes */


/* External includes */


/* Project includes */


namespace Kratos
{
	
    /**@name Kratos Globals */
    /*@{ */
	
	
    /*@} */
    /**@name Type Definitions */       
    /*@{ */
	
	
    /*@} */
	
	
    /**@name  Enum's */       
    /*@{ */
	
	
    /*@} */
    /**@name  Functions */       
    /*@{ */
	
	
	
    /*@} */
    /**@name Kratos Classes */
    /*@{ */
	
	///Various mathematical utilitiy functions.
        /**
     * Defines several utility functions
         */
    template<class TDataType>
            class MathUtils
    {
        public:
            /**@name Type Definitions */       
            /*@{ */
            typedef Matrix MatrixType;
		
            typedef Vector VectorType;
		
            typedef unsigned int IndexType;
		
            /*@} */
		/**@name Life Cycle 
                 */    
            /*@{ */
		
		/** Constructor.
                 */
		
		
		/** Destructor.
                 */
		
            /*@} */
		/**@name Operators 
                 */  
            /*@{ */
		
		
            /*@} */
            /**@name Operations */
            /*@{ */
		
		//***********************************************************************
		//***********************************************************************
            static TDataType Distance(const TDataType& rFirstData, const TDataType& rSecondData)
            {return rFirstData.Distance(rSecondData);}
		
		//***********************************************************************
		//***********************************************************************
            static TDataType Abs(const TDataType& rData)
            {return rData > TDataType(0) ? rData : -rData;}

            static TDataType Min(const TDataType& rValue1, const TDataType& rValue2)
            {return rValue1 > rValue2 ? rValue2 : rValue1;}

            static TDataType Max(const TDataType& rValue1, const TDataType& rValue2)
            {return rValue1 > rValue2 ? rValue1 : rValue2;}
		
		//***********************************************************************
		/** 
            inverts matrices of order 2 and 3
                 */
		//***********************************************************************
            static inline void InvertMatrix(const MatrixType& InputMatrix, 
                                            MatrixType& InvertedMatrix,
                                            TDataType& InputMatrixDet)
            {
                unsigned int size = InputMatrix.size2();
			
                if (size==2) InvertMatrix2(InputMatrix,InvertedMatrix,InputMatrixDet);
                else InvertMatrix3(InputMatrix,InvertedMatrix,InputMatrixDet);			
            }
		
		//***********************************************************************
		/**
		Inversion of a 3*3 matrix (no check is performed on matrix size):
            - InputMatrix is the input matrix (unchanged at output)
            - InvertedMatrix is the inverse of the input matrix
            - InputMatrixDet is the determinant of the input matrix
                 */
		//***********************************************************************
            static void InvertMatrix3(const MatrixType& InputMatrix, 
                                      MatrixType& InvertedMatrix,
                                      TDataType& InputMatrixDet) //VERIFIED!!!
            {
                KRATOS_TRY
                        if(InvertedMatrix.size1() != 3 || InvertedMatrix.size2() != 3)
                        InvertedMatrix.resize(3,3);
			
			//filling the inverted matrix with the algebraic complements
			//first column
                InvertedMatrix(0,0) = InputMatrix(1,1)*InputMatrix(2,2) - InputMatrix(1,2)*InputMatrix(2,1);
                InvertedMatrix(1,0) = -InputMatrix(1,0)*InputMatrix(2,2) + InputMatrix(1,2)*InputMatrix(2,0);
                InvertedMatrix(2,0) = InputMatrix(1,0)*InputMatrix(2,1) - InputMatrix(1,1)*InputMatrix(2,0);
			
			//second column
                InvertedMatrix(0,1) = -InputMatrix(0,1)*InputMatrix(2,2) + InputMatrix(0,2)*InputMatrix(2,1);
                InvertedMatrix(1,1) = InputMatrix(0,0)*InputMatrix(2,2) - InputMatrix(0,2)*InputMatrix(2,0);
                InvertedMatrix(2,1) = -InputMatrix(0,0)*InputMatrix(2,1) + InputMatrix(0,1)*InputMatrix(2,0);
			
			//third column
                InvertedMatrix(0,2) = InputMatrix(0,1)*InputMatrix(1,2) - InputMatrix(0,2)*InputMatrix(1,1);
                InvertedMatrix(1,2) = -InputMatrix(0,0)*InputMatrix(1,2) + InputMatrix(0,2)*InputMatrix(1,0);
                InvertedMatrix(2,2) = InputMatrix(0,0)*InputMatrix(1,1) - InputMatrix(0,1)*InputMatrix(1,0);
			
			//calculation of determinant (of the input matrix)
                InputMatrixDet = InputMatrix(0,0)*InvertedMatrix(0,0) + InputMatrix(0,1)*InvertedMatrix(1,0) + InputMatrix(0,2)*InvertedMatrix(2,0);
			
			//finalizing the calculation of the inverted matrix
                InvertedMatrix /= InputMatrixDet;
                KRATOS_CATCH("")
            }
		
		//***********************************************************************
		/**
		Inversion of a 2*2 matrix (no check is performed on matrix size):
            - InputMatrix is the input matrix (unchanged at output)
            - InvertedMatrix is the inverse of the input matrix
            - InputMatrixDet is the determinant of the input matrix
                 */
		//***********************************************************************
            static void InvertMatrix2(const MatrixType& InputMatrix, 
                                      MatrixType& InvertedMatrix,
                                      TDataType& InputMatrixDet) 
            {
                KRATOS_TRY
                        InvertedMatrix.resize(2,2);
			
                InputMatrixDet = InputMatrix(0,0)*InputMatrix(1,1)-InputMatrix(0,1)*InputMatrix(1,0);
			
                InvertedMatrix(0,0) =  InputMatrix(1,1); InvertedMatrix(0,1) = -InputMatrix(0,1);
                InvertedMatrix(1,0) = -InputMatrix(1,0); InvertedMatrix(1,1) =  InputMatrix(0,0);
			
                InvertedMatrix/=InputMatrixDet;
                KRATOS_CATCH("")
            }
		
		//***********************************************************************
		/**calculates the determinant of a matrix of dimension 2*2
            (no check performed on matrix size)*/
		//***********************************************************************
            static inline TDataType Det2(const MatrixType& A)
            {
                return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
            }
		
		//***********************************************************************
		/**calculates the determinant of a matrix of dimension 3*3
            (no check performed on matrix size)*/
		//***********************************************************************
            static inline TDataType Det3(const MatrixType& A) 
            {			
			//calculating the algebraic complements to the first line
                double a = A(1,1)*A(2,2) - A(1,2)*A(2,1);
                double b = A(1,0)*A(2,2) - A(1,2)*A(2,0);
                double c = A(1,0)*A(2,1) - A(1,1)*A(2,0);
				
                return A(0,0)*a - A(0,1)*b + A(0,2)*c;
            }

		//***********************************************************************
		/**calculates the determinant of a matrix of dimension 2*2 or 3*3
            (no check performed on matrix size)*/
		//***********************************************************************
            static inline TDataType Det(const boost::numeric::ublas::bounded_matrix<double,2,2>& A)
            {
                return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
            }
            static inline TDataType Det(const boost::numeric::ublas::bounded_matrix<double,3,3>& A) 
            {			
			//calculating the algebraic complements to the first line
                double a = A(1,1)*A(2,2) - A(1,2)*A(2,1);
                double b = A(1,0)*A(2,2) - A(1,2)*A(2,0);
                double c = A(1,0)*A(2,1) - A(1,1)*A(2,0);
				
                return A(0,0)*a - A(0,1)*b + A(0,2)*c;
            }
            static inline TDataType Det(const MatrixType& A) 
            {	
                TDataType Det;
                if (A.size1()==2)
                    Det = Det2(A);
                else
                    Det = Det3(A);

                return Det;

            }
		
		//***********************************************************************
		/** performs the dot product of two vectors of dimension 3 
            (no check performed on vector sizes)*/
		//***********************************************************************
            static inline TDataType Dot3(Vector& a, Vector& b)
            {
                return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
            }

		//***********************************************************************
		/** performs the dot product of two vectors of arbitrary size
            (no check performed on vector sizes)*/
		//***********************************************************************
            static inline TDataType Dot(const Vector& FirstVector, const Vector& SecondVector)
            {
                Vector::const_iterator i = FirstVector.begin();
                Vector::const_iterator j = SecondVector.begin();
                TDataType temp();
                while(i != FirstVector.end())
                    temp += *i++ * *j++;
                return temp;
			//return std::inner_product(FirstVector.begin(), FirstVector.end(), SecondVector.begin(), TDataType());
            }	

		//***********************************************************************
		/* 
            calculates the norm of vector "a" which is assumed to be of size 3
            (no check is performed on the vector's size)
                */
		//***********************************************************************
            static inline TDataType Norm3(Vector& a)
            {
                TDataType temp = pow(a[0],2) + pow(a[1],2) + pow(a[2],2);
                temp = sqrt(temp);
                return temp;
            }

		//***********************************************************************
		/* 
            calculates the norm of vector "a" 
                */
		//***********************************************************************

            static inline TDataType Norm(const Vector& a)
            {
                Vector::const_iterator i = a.begin();
                TDataType temp = TDataType();
                while(i != a.end())
                {
                    temp += (*i) * (*i);
                    i++;
                }
                return sqrt(temp);
            }		
		
		//***********************************************************************
		/* 
            performs the vector product of the two input vectors a,b
            a,b are assumed to be of size 3 (no check is performed on vector sizes)
                */
		//***********************************************************************
            static inline Vector CrossProduct(Vector& a, Vector& b)
            {
                Vector c(3);
			
                c[0] = a[1]*b[2] - a[2]*b[1];
                c[1] = a[2]*b[0] - a[0]*b[2];
                c[2] = a[0]*b[1] - a[1]*b[0];
			
                return c;
            }	

		//identical but it assumes that the output vector is given already the correct size
            static inline void CrossProduct(array_1d<double, 3>& c, const array_1d<double, 3>& a, const array_1d<double, 3>& b)
            {
                c[0] = a[1]*b[2] - a[2]*b[1];
                c[1] = a[2]*b[0] - a[0]*b[2];
                c[2] = a[0]*b[1] - a[1]*b[0];
            }

		//***********************************************************************
		/** 
		returns a matrix :
            A = a.tensorproduct.b
            a,b are assumed to be of order 3, no check is performed on the size of the vectors
                 */
		//***********************************************************************
            static inline MatrixType TensorProduct3(Vector& a, Vector& b)
            {
                MatrixType A(3,3);
                A(0,0)=a[0]*b[0]; A(0,1)=a[0]*b[1]; A(0,2)=a[0]*b[2];
                A(1,0)=a[1]*b[0]; A(1,1)=a[1]*b[1]; A(1,2)=a[1]*b[2];
                A(2,0)=a[2]*b[0]; A(2,1)=a[2]*b[1]; A(2,2)=a[2]*b[2];
                return A;
            }
		
		//***********************************************************************
		/** 
            "InputMatrix" is ADDED to "Destination" matrix starting from 
            InitialRow and InitialCol of the destination matrix
            "Destination" is assumed to be able to contain the "input matrix"
            (no check is performed on the bounds)
                 */
		//***********************************************************************
            static inline void  AddMatrix(
                                          MatrixType& Destination,
                                          MatrixType& InputMatrix, 
                                          int InitialRow,
                                          int InitialCol)
            {
                KRATOS_TRY
                        for(unsigned int i=0; i<InputMatrix.size1(); i++)
                        for(unsigned int j=0; j<InputMatrix.size2(); j++) 
                        Destination(InitialRow+i, InitialCol+j) += InputMatrix(i,j);
                KRATOS_CATCH("")
            }
		//***********************************************************************
		/** 
            "InputMatrix" is SUBTRACTED to "Destination" matrix starting from 
            InitialRow and InitialCol of the destination matrix
            "Destination" is assumed to be able to contain the "input matrix"
            (no check is performed on the bounds)
                 */
		//***********************************************************************
            static inline void  SubtractMatrix(
                                               MatrixType& Destination,
                                               MatrixType& InputMatrix, 
                                               int InitialRow,
                                               int InitialCol)
            {
                KRATOS_TRY
                        for(unsigned int i=0; i<InputMatrix.size1(); i++)
                        for(unsigned int j=0; j<InputMatrix.size2(); j++) 
                        Destination(InitialRow+i, InitialCol+j) -= InputMatrix(i,j);
                KRATOS_CATCH("")
            }
		
		//***********************************************************************
		/** 
            "InputMatrix" is WRITTEN on "Destination" matrix starting from 
            InitialRow and InitialCol of the destination matrix
            "Destination" is assumed to be able to contain the "input matrix"
            (no check is performed on the bounds)
            ATTENTION: Destination is overwritten!!
                 */
		//***********************************************************************
            static inline void  WriteMatrix(
                                            MatrixType& Destination,
                                            MatrixType& InputMatrix, 
                                            int InitialRow,
                                            int InitialCol)
            {
                KRATOS_TRY
                        for(unsigned int i=0; i<InputMatrix.size1(); i++)
                        for(unsigned int j=0; j<InputMatrix.size2(); j++) 
                        Destination(InitialRow+i, InitialCol+j) = InputMatrix(i,j);
                KRATOS_CATCH("")
            }
		
		//***********************************************************************
		//***********************************************************************
		//performs the Kroneker product of the Reduced Matrix with the identity matrix of 
		//size "dimension"
            static inline void  ExpandReducedMatrix(
                    MatrixType& Destination,
                    MatrixType& ReducedMatrix,
                    unsigned int dimension)
            {
                KRATOS_TRY
                        unsigned int size=ReducedMatrix.size2();
			
                for (unsigned int i=0;i<size;i++)
                {
                    unsigned int rowindex = i*dimension;
                    for (unsigned int j=0;j<size;j++)
                    {
                        unsigned int colindex = j*dimension;
                        for(unsigned int ii=0;ii<dimension;ii++)
                            Destination(rowindex+ii,colindex+ii) = ReducedMatrix(i,j);
                    }
                }
                KRATOS_CATCH("")
            }
		//***********************************************************************
		//***********************************************************************
		//performs the Kroneker product of the Reduced Matrix with the identity matrix of 
		//size "dimension" ADDING to the destination matrix
            static inline void  ExpandAndAddReducedMatrix(
                    MatrixType& Destination,
                    MatrixType& ReducedMatrix,
                    unsigned int dimension)
            {
                KRATOS_TRY
                        unsigned int size=ReducedMatrix.size2();
			
                for (unsigned int i=0;i<size;i++)
                {
                    int rowindex = i*dimension;
                    for (unsigned int j=0;j<size;j++)
                    {
                        unsigned int colindex = j*dimension;
                        for(unsigned int ii=0;ii<dimension;ii++)
                            Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
                    }
                }
                KRATOS_CATCH("")
            }
		
		//***********************************************************************
		/** 
            performs x += coeff*y;
            no check on bounds is performed
                 */
		//***********************************************************************
            static inline void  VecAdd(
                                       Vector& x,
                                       TDataType coeff,
                                       Vector& y)
            {
                KRATOS_TRY
                        unsigned int size=x.size();
			
                for (unsigned int i=0;i<size;i++)
                {
                    x[i] += coeff * y[i];
                }
                KRATOS_CATCH("")
            }
		
		//***********************************************************************
		//***********************************************************************
            static inline MatrixType StressVectorToTensor(const Vector& StressVector)
            {
                KRATOS_TRY
                        Matrix StressTensor;
			
                if (StressVector.size()==3)
                {
                    StressTensor.resize(2,2);
                    StressTensor(0,0) = StressVector[0]; StressTensor(0,1) = StressVector[2];
                    StressTensor(1,0) = StressVector[2]; StressTensor(1,1) = StressVector[1];
                }
                else if (StressVector.size()==6)
                {
                    StressTensor.resize(3,3);
                    StressTensor(0,0) = StressVector[0]; StressTensor(0,1) = StressVector[3]; StressTensor(0,2) = StressVector[5];
                    StressTensor(1,0) = StressVector[3]; StressTensor(1,1) = StressVector[1]; StressTensor(1,2) = StressVector[4];
                    StressTensor(2,0) = StressVector[5]; StressTensor(2,1) = StressVector[4]; StressTensor(2,2) = StressVector[2];
                }
			
                return StressTensor;
                KRATOS_CATCH("")
				
            }
            /*@} */  
            /**@name Acces */
            /*@{ */
		
		
            /*@} */
            /**@name Inquiry */
            /*@{ */
		
		
            /*@} */      
            /**@name Friends */
            /*@{ */
		
		
            /*@} */
		
        private:
            /**@name Static Member Variables */
            /*@{ */
		
		
            /*@} */
            /**@name Member Variables */
            /*@{ */
		
		
            /*@} */
            /**@name Private Operators*/
            /*@{ */
		
		
            /*@} */
            /**@name Private Operations*/
            /*@{ */
		
		
            /*@} */
            /**@name Private  Acces */
            /*@{ */
		
		
            /*@} */     
            /**@name Private Inquiry */
            /*@{ */
		
		
            /*@} */   
            /**@name Un accessible methods */
            /*@{ */
		
            MathUtils(void);
		
            MathUtils(MathUtils& rSource);
		
		
            /*@} */   
		
    }; /* Class ClassName */
	
    /*@} */
	
    /**@name Type Definitions */       
    /*@{ */
	
	
    /*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_MATH_UTILS  defined */

