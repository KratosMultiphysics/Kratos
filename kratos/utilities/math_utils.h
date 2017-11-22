//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//  Collaborator:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MATH_UTILS )
#define  KRATOS_MATH_UTILS


/* System includes */


/* External includes */
#include <cmath>
#include "includes/ublas_interface.h"
#include "containers/array_1d.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

///Various mathematical utilitiy functions.
/**
 * Defines several utility functions
 */
template<class TDataType>
class MathUtils
{
public:

    ///@name Type Definitions
    ///@{
    
    typedef Matrix MatrixType;

    typedef Vector VectorType;

    typedef std::size_t SizeType;
    
    typedef unsigned int IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /* Constructor */


    /** Destructor */

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This function calculates the number of elements between first and last.
     * @param rFirstData: First element
     * @param rSecondData: Second element
     * @return Distance: Number of elements
     */
    
    static TDataType Distance(
        const TDataType& rFirstData, 
        const TDataType& rSecondData
        )
    {
        return rFirstData.Distance(rSecondData);
    }

    /**
     * In geometry, Heron's formula (sometimes called Hero's formula), named after Hero of Alexandria, gives the area of a triangle by requiring no arbitrary choice of side as base or vertex as origin, contrary to other formulas for the area of a triangle, such as half the base times the height or half the norm of a cross product of two sides.
     * @param a: First length
     * @param b: Second length
     * @param c: Third length
     * @return Heron solution: Heron's formula states that the area of a triangle whose sides have lengths a, b, and c
     */
    
    template<bool check>// = false>
    static inline double Heron(
        double a, 
        double b, 
        double c
        )
    {
        const double s = 0.5 * (a + b + c);
        const double A2 = s * (s - a) * (s - b) * (s - c);
        if(check == true)
        {
            if(A2 < 0.0)
            {
                KRATOS_ERROR << "The square of area is negative, probably the triangle is in bad shape:" << A2 << std::endl;
            }
            else
            {
                return std::sqrt(A2);
            }
        }
        else
        {
            return std::sqrt(std::abs(A2));
        }
    }

    /**
     * It gives you the absolute value of a given value
     * @param rData: The value to compute the absolute value
     * @return The absolute value of rData
     */
    
    static TDataType Abs(const TDataType& rData)
    {
        return rData > TDataType(0) ? rData : -rData;
    }

    /**
     * It gives you the minimum value between two values
     * @param rValue1: The first value
     * @param rValue2: The second value
     * @return The minimum value
     */
        
    static TDataType Min(
        const TDataType& rValue1, 
        const TDataType& rValue2
        )
    {
        return rValue1 > rValue2 ? rValue2 : rValue1;
    }

    /**
     * It gives you the maximum value between two values
     * @param rValue1: The first value
     * @param rValue2: The second value
     * @return The maximum value
     */
        
    static TDataType Max(
        const TDataType& rValue1, 
        const TDataType& rValue2
        )
    {
        return rValue1 > rValue2 ? rValue1 : rValue2;
    }

    /**
     * Calculates the determinant of a 2x2, 3x3 and 4x4 matrix (using bounded matrix for performance) 
     * @param A: The matrix to calculate
     * @return DetA: The determinant of the matrix
     */

    template<unsigned int TDim>
    static inline TDataType DetMat(const boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim>& InputMatrix)
    {
        TDataType InputMatrixDet;
        
        if (TDim == 1)
        {
            InputMatrixDet = InputMatrix(0, 0);
        }
        else if (TDim == 2)
        {
            InputMatrixDet = InputMatrix(0, 0) * InputMatrix(1, 1) - InputMatrix(0, 1) * InputMatrix(1, 0);
        }
        else if (TDim == 3)
        {
            InputMatrixDet = InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 2)
                           + InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(0, 2)
                           + InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 0)
                           - InputMatrix(2, 0) * InputMatrix(1, 1) * InputMatrix(0, 2)
                           - InputMatrix(2, 1) * InputMatrix(1, 2) * InputMatrix(0, 0)
                           - InputMatrix(1, 0) * InputMatrix(0, 1) * InputMatrix(2,2);
        }
        else
        {
            InputMatrixDet = InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 0) - InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 1) + InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 2) + InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 2) + InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 2) - InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 2) + InputMatrix(0, 3) * (InputMatrix(1, 2) * InputMatrix(2, 1) * InputMatrix(3, 0) - InputMatrix(1, 1) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(1, 2) * InputMatrix(2, 0) * InputMatrix(3, 1) + InputMatrix(1, 0) * InputMatrix(2, 2) * InputMatrix(3, 1) + InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 2) - InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 2)) + (InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 0) - InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 1) - InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 2) + InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 2)) * InputMatrix(3, 3) + InputMatrix(0, 2) * (-(InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 0)) + InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 0) + InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 1) - InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 3) + InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 3));
        }
        
        return InputMatrixDet;
    }
    
    /**
     * Calculates the inverse of a 2x2, 3x3 and 4x4 matrices (using bounded matrix for performance) 
     * @param InputMatrix: The matrix to invert
     * @return InvertMatrix: The inverted matrix
     */

    template<unsigned int TDim>
    static inline boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> InvertMatrix(
            const boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim>& InputMatrix,
            TDataType& InputMatrixDet,
            const TDataType Tolerance = std::numeric_limits<double>::epsilon()
            )
    {
        boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> InvertedMatrix;
        
        /* Compute Determinant of the matrix */
        InputMatrixDet = DetMat<TDim>(InputMatrix);
        
        if (std::abs(InputMatrixDet) < Tolerance)
        {
            KRATOS_WATCH(InputMatrix);
            KRATOS_ERROR << " Determinant of the matrix is zero or almost zero!!!, InputMatrixDet = " << InputMatrixDet << std::endl;
        }
        
        if (TDim == 1)
        {
            InvertedMatrix(0, 0) = 1.0 / InputMatrixDet;
        }
        else if (TDim == 2)
        {            
            /* Compute inverse of the Matrix */
            InvertedMatrix(0, 0) =   InputMatrix(1, 1) / InputMatrixDet;
            InvertedMatrix(0, 1) = - InputMatrix(0, 1) / InputMatrixDet;
            InvertedMatrix(1, 0) = - InputMatrix(1, 0) / InputMatrixDet;
            InvertedMatrix(1, 1) =   InputMatrix(0, 0) / InputMatrixDet;
        }
        else if (TDim == 3)
        {
            /* Compute inverse of the Matrix */
            // First column
            InvertedMatrix(0, 0) =   (InputMatrix(1, 1) * InputMatrix(2, 2) - InputMatrix(1, 2) * InputMatrix(2, 1)) / InputMatrixDet;
            InvertedMatrix(1, 0) = - (InputMatrix(1, 0) * InputMatrix(2, 2) - InputMatrix(2, 0) * InputMatrix(1, 2)) / InputMatrixDet;
            InvertedMatrix(2, 0) =   (InputMatrix(1, 0) * InputMatrix(2, 1) - InputMatrix(1, 1) * InputMatrix(2, 0)) / InputMatrixDet;

            // Second column
            InvertedMatrix(0, 1) = - (InputMatrix(0, 1) * InputMatrix(2, 2) - InputMatrix(0, 2) * InputMatrix(2, 1)) / InputMatrixDet;
            InvertedMatrix(1, 1) =   (InputMatrix(0, 0) * InputMatrix(2, 2) - InputMatrix(0, 2) * InputMatrix(2, 0)) / InputMatrixDet;
            InvertedMatrix(2, 1) = - (InputMatrix(0, 0) * InputMatrix(2, 1) - InputMatrix(0, 1) * InputMatrix(2, 0)) / InputMatrixDet;

            // Third column
            InvertedMatrix(0, 2) =   (InputMatrix(0, 1) * InputMatrix(1, 2) - InputMatrix(0, 2) * InputMatrix(1, 1)) / InputMatrixDet;
            InvertedMatrix(1, 2) = - (InputMatrix(0, 0) * InputMatrix(1, 2) - InputMatrix(0, 2) * InputMatrix(1, 0)) / InputMatrixDet;
            InvertedMatrix(2, 2) =   (InputMatrix(0, 0) * InputMatrix(1, 1) - InputMatrix(1, 0) * InputMatrix(0, 1)) / InputMatrixDet;
        }
        else if (TDim == 4)
        {
            /* Compute inverse of the Matrix */
            // First column
            InvertedMatrix(0, 0) = -(InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 1)) + InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 1) + InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 2) - InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 2) - InputMatrix(1, 2) * InputMatrix(2, 1) * InputMatrix(3, 3) + InputMatrix(1, 1) * InputMatrix(2, 2) * InputMatrix(3, 3);
            InvertedMatrix(0, 1) = InputMatrix(0, 3) * InputMatrix(2, 2) * InputMatrix(3, 1) - InputMatrix(0, 2) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(0, 3) * InputMatrix(2, 1) * InputMatrix(3, 2) + InputMatrix(0, 1) * InputMatrix(2, 3) * InputMatrix(3, 2) + InputMatrix(0, 2) * InputMatrix(2, 1) * InputMatrix(3, 3) - InputMatrix(0, 1) * InputMatrix(2, 2) * InputMatrix(3, 3);
            InvertedMatrix(0, 2) = -(InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(3, 1)) + InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(3, 1) + InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(3, 2) - InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(3, 2) - InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(3, 3) + InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(3, 3);
            InvertedMatrix(0, 3) = InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(2, 1) - InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(2, 1) - InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(2, 2) + InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 2) + InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(2, 3) - InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 3);

            // Second column
            InvertedMatrix(1, 0) = InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 0) - InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 2) + InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 2) + InputMatrix(1, 2) * InputMatrix(2, 0) * InputMatrix(3, 3) - InputMatrix(1, 0) * InputMatrix(2, 2) * InputMatrix(3, 3);
            InvertedMatrix(1, 1) = -(InputMatrix(0, 3) * InputMatrix(2, 2) * InputMatrix(3, 0)) + InputMatrix(0, 2) * InputMatrix(2, 3) * InputMatrix(3, 0) + InputMatrix(0, 3) * InputMatrix(2, 0) * InputMatrix(3, 2) - InputMatrix(0, 0) * InputMatrix(2, 3) * InputMatrix(3, 2) - InputMatrix(0, 2) * InputMatrix(2, 0) * InputMatrix(3, 3) + InputMatrix(0, 0) * InputMatrix(2, 2) * InputMatrix(3, 3);
            InvertedMatrix(1, 2) = InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(3, 0) - InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(3, 0) - InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(3, 2) + InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(3, 2) + InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(3, 3) - InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(3, 3);
            InvertedMatrix(1, 3) = -(InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(2, 0)) + InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(2, 0) + InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(2, 2) - InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 2) - InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(2, 3) + InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 3);
            
            // Third column
            InvertedMatrix(2, 0) = -(InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 0)) + InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 0) + InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 1) - InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 3) + InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 3);
            InvertedMatrix(2, 1) = InputMatrix(0, 3) * InputMatrix(2, 1) * InputMatrix(3, 0) - InputMatrix(0, 1) * InputMatrix(2, 3) * InputMatrix(3, 0) - InputMatrix(0, 3) * InputMatrix(2, 0) * InputMatrix(3, 1) + InputMatrix(0, 0) * InputMatrix(2, 3) * InputMatrix(3, 1) + InputMatrix(0, 1) * InputMatrix(2, 0) * InputMatrix(3, 3) - InputMatrix(0, 0) * InputMatrix(2, 1) * InputMatrix(3, 3);
            InvertedMatrix(2, 2) = -(InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(3, 0)) + InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(3, 0) + InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(3, 1) - InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(3, 1) - InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(3, 3) + InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(3, 3);
            InvertedMatrix(2, 3) = InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(2, 0) - InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 0) - InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(2, 1) + InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 1) + InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 3) - InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 3);
            
            // Fourth column
            InvertedMatrix(3, 0) = InputMatrix(1, 2) * InputMatrix(2, 1) * InputMatrix(3, 0) - InputMatrix(1, 1) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(1, 2) * InputMatrix(2, 0) * InputMatrix(3, 1) + InputMatrix(1, 0) * InputMatrix(2, 2) * InputMatrix(3, 1) + InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 2) - InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 2);
            InvertedMatrix(3, 1) = -(InputMatrix(0, 2) * InputMatrix(2, 1) * InputMatrix(3, 0)) + InputMatrix(0, 1) * InputMatrix(2, 2) * InputMatrix(3, 0) + InputMatrix(0, 2) * InputMatrix(2, 0) * InputMatrix(3, 1) - InputMatrix(0, 0) * InputMatrix(2, 2) * InputMatrix(3, 1) - InputMatrix(0, 1) * InputMatrix(2, 0) * InputMatrix(3, 2) + InputMatrix(0, 0) * InputMatrix(2, 1) * InputMatrix(3, 2);
            InvertedMatrix(3, 2) = InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(3, 0) - InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(3, 0) - InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(3, 1) + InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(3, 1) + InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(3, 2) - InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(3, 2);
            InvertedMatrix(3, 3) = -(InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(2, 0)) + InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 0) + InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(2, 1) - InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 1) - InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 2) + InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 2);
            
            InvertedMatrix /= InputMatrixDet;
        }
        else
        {
            KRATOS_ERROR << "::WARNING: Size not implemented. Size: " << TDim << std::endl;
        }
        
        return InvertedMatrix;
    }
    
    /**
     * It inverts non square matrices (https://en.wikipedia.org/wiki/Inverse_element#Matrices)
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */

    static void GeneralizedInvertMatrix(
        const MatrixType& InputMatrix,
        MatrixType& InvertedMatrix,
        TDataType& InputMatrixDet
        )
    {
        const unsigned int size_1 = InputMatrix.size1();
        const unsigned int size_2 = InputMatrix.size2();
        
        if (size_1 == size_2)
        {
            InvertMatrix(InputMatrix, InvertedMatrix, InputMatrixDet);
        }
        else if (size_1 < size_2) // Right inverse
        {
            if (InvertedMatrix.size1() != size_2 || InvertedMatrix.size2() != size_1)
            {
                InvertedMatrix.resize(size_2, size_1, false);
            }
            const Matrix aux = prod(InputMatrix, trans(InputMatrix));
            Matrix auxInv;
            InvertMatrix(aux, auxInv, InputMatrixDet);
	    InputMatrixDet = std::sqrt(InputMatrixDet);
            noalias(InvertedMatrix) = prod(trans(InputMatrix), auxInv);
        }
        else // Left inverse
        {
            if (InvertedMatrix.size1() != size_2 || InvertedMatrix.size2() != size_1)
            {
                InvertedMatrix.resize(size_2, size_1, false);
            }
            const Matrix aux = prod(trans(InputMatrix), InputMatrix);
            Matrix auxInv;
            InvertMatrix(aux, auxInv, InputMatrixDet);
	    InputMatrixDet = std::sqrt(InputMatrixDet);
            noalias(InvertedMatrix) = prod(auxInv, trans(InputMatrix));
        }
    }
    
    /**
     * It inverts matrices of order 2, 3 and 4
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */

    static void InvertMatrix(
        const MatrixType& InputMatrix,
        MatrixType& InvertedMatrix,
        TDataType& InputMatrixDet
        )
    {
        const unsigned int size = InputMatrix.size2();

        if(size == 1)
        {
            if(InvertedMatrix.size1() != 1 || InvertedMatrix.size2() != 1)
            {
                InvertedMatrix.resize(1,1,false);
            }
            InvertedMatrix(0,0) = 1.0/InputMatrix(0,0);
            InputMatrixDet = InputMatrix(0,0);
        }
        else if (size == 2) 
        {
            InvertMatrix2(InputMatrix, InvertedMatrix, InputMatrixDet);
        }
        else if (size == 3) 
        {
            InvertMatrix3(InputMatrix, InvertedMatrix, InputMatrixDet);
        }
        else if (size == 4)
        {
            InvertMatrix4(InputMatrix, InvertedMatrix, InputMatrixDet);
        }
        else
        {   
            const unsigned int size1 = InputMatrix.size1();
            const unsigned int size2 = InputMatrix.size2();
            if(InvertedMatrix.size1() != size1 || InvertedMatrix.size2() != size2)
            {
                InvertedMatrix.resize(size1, size2,false);
            }
            
            int singular = 0;
            using namespace boost::numeric::ublas;
            typedef permutation_matrix<SizeType> pmatrix;
            Matrix A(InputMatrix);
            pmatrix pm(A.size1());
            singular = lu_factorize(A,pm);
            InvertedMatrix.assign( IdentityMatrix(A.size1()));
            lu_substitute(A, pm, InvertedMatrix);
             
            // Calculating determinant
            InputMatrixDet = 1.0;
            
            for (unsigned int i = 0; i < A.size1();i++)
            {
                unsigned int ki = pm[i] == i ? 0 : 1;
                InputMatrixDet *= std::pow(-1.0, ki) * A(i,i);
            }
            
            if (singular == 1)
            {
                KRATOS_ERROR << "::WARNING: Matrix is singular: " << InputMatrix << std::endl;
            }
        }
    }

    /**
     * It inverts matrices of order 2 //VERIFIED!!!
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */
    
    static void InvertMatrix2(
        const MatrixType& InputMatrix,
        MatrixType& InvertedMatrix,
        TDataType& InputMatrixDet
        )
    {
        KRATOS_TRY;
        
        if(InvertedMatrix.size1() != 2 || InvertedMatrix.size2() != 2)
        {
            InvertedMatrix.resize(2,2,false);
        }

        InputMatrixDet = InputMatrix(0,0)*InputMatrix(1,1)-InputMatrix(0,1)*InputMatrix(1,0);

        InvertedMatrix(0,0) =  InputMatrix(1,1);
        InvertedMatrix(0,1) = -InputMatrix(0,1);
        InvertedMatrix(1,0) = -InputMatrix(1,0);
        InvertedMatrix(1,1) =  InputMatrix(0,0);

        InvertedMatrix/=InputMatrixDet;
        
        KRATOS_CATCH("");
    }
    
    /**
     * It inverts matrices of order 3 //VERIFIED!!!
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */
    
    static void InvertMatrix3(
        const MatrixType& InputMatrix,
        MatrixType& InvertedMatrix,
        TDataType& InputMatrixDet
        ) 
    {
        KRATOS_TRY;
        
        if(InvertedMatrix.size1() != 3 || InvertedMatrix.size2() != 3)
        {
            InvertedMatrix.resize(3,3,false);
        }

        // Filling the inverted matrix with the algebraic complements
        // First column
        InvertedMatrix(0,0) = InputMatrix(1,1)*InputMatrix(2,2) - InputMatrix(1,2)*InputMatrix(2,1);
        InvertedMatrix(1,0) = -InputMatrix(1,0)*InputMatrix(2,2) + InputMatrix(1,2)*InputMatrix(2,0);
        InvertedMatrix(2,0) = InputMatrix(1,0)*InputMatrix(2,1) - InputMatrix(1,1)*InputMatrix(2,0);

        // Second column
        InvertedMatrix(0,1) = -InputMatrix(0,1)*InputMatrix(2,2) + InputMatrix(0,2)*InputMatrix(2,1);
        InvertedMatrix(1,1) = InputMatrix(0,0)*InputMatrix(2,2) - InputMatrix(0,2)*InputMatrix(2,0);
        InvertedMatrix(2,1) = -InputMatrix(0,0)*InputMatrix(2,1) + InputMatrix(0,1)*InputMatrix(2,0);

        // Third column
        InvertedMatrix(0,2) = InputMatrix(0,1)*InputMatrix(1,2) - InputMatrix(0,2)*InputMatrix(1,1);
        InvertedMatrix(1,2) = -InputMatrix(0,0)*InputMatrix(1,2) + InputMatrix(0,2)*InputMatrix(1,0);
        InvertedMatrix(2,2) = InputMatrix(0,0)*InputMatrix(1,1) - InputMatrix(0,1)*InputMatrix(1,0);

        // Calculation of determinant (of the input matrix)
        InputMatrixDet = InputMatrix(0,0)*InvertedMatrix(0,0) + InputMatrix(0,1)*InvertedMatrix(1,0) + InputMatrix(0,2)*InvertedMatrix(2,0);

        // Finalizing the calculation of the inverted matrix
        InvertedMatrix /= InputMatrixDet;
        
        KRATOS_CATCH("")
    }
    
    /**
     * It inverts matrices of order 4
     * @param InputMatrix: Is the input matrix (unchanged at output)
     * @return InvertedMatrix: Is the inverse of the input matrix
     * @return InputMatrixDet: Is the determinant of the input matrix
     */
    
    static void InvertMatrix4(
        const MatrixType& InputMatrix,
        MatrixType& InvertedMatrix,
        TDataType& InputMatrixDet
        ) 
    {
        KRATOS_TRY;
        
    if (InvertedMatrix.size1() != 4 || InvertedMatrix.size2() != 4) 
    {
        InvertedMatrix.resize(4, 4, false);
    }

    /* Compute inverse of the Matrix */
    // First column
    InvertedMatrix(0, 0) = -(InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 1)) + InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 1) + InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 2) - InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 2) - InputMatrix(1, 2) * InputMatrix(2, 1) * InputMatrix(3, 3) + InputMatrix(1, 1) * InputMatrix(2, 2) * InputMatrix(3, 3);
    InvertedMatrix(0, 1) = InputMatrix(0, 3) * InputMatrix(2, 2) * InputMatrix(3, 1) - InputMatrix(0, 2) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(0, 3) * InputMatrix(2, 1) * InputMatrix(3, 2) + InputMatrix(0, 1) * InputMatrix(2, 3) * InputMatrix(3, 2) + InputMatrix(0, 2) * InputMatrix(2, 1) * InputMatrix(3, 3) - InputMatrix(0, 1) * InputMatrix(2, 2) * InputMatrix(3, 3);
    InvertedMatrix(0, 2) = -(InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(3, 1)) + InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(3, 1) + InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(3, 2) - InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(3, 2) - InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(3, 3) + InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(3, 3);
    InvertedMatrix(0, 3) = InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(2, 1) - InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(2, 1) - InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(2, 2) + InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 2) + InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(2, 3) - InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 3);

    // Second column
    InvertedMatrix(1, 0) = InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 0) - InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 2) + InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 2) + InputMatrix(1, 2) * InputMatrix(2, 0) * InputMatrix(3, 3) - InputMatrix(1, 0) * InputMatrix(2, 2) * InputMatrix(3, 3);
    InvertedMatrix(1, 1) = -(InputMatrix(0, 3) * InputMatrix(2, 2) * InputMatrix(3, 0)) + InputMatrix(0, 2) * InputMatrix(2, 3) * InputMatrix(3, 0) + InputMatrix(0, 3) * InputMatrix(2, 0) * InputMatrix(3, 2) - InputMatrix(0, 0) * InputMatrix(2, 3) * InputMatrix(3, 2) - InputMatrix(0, 2) * InputMatrix(2, 0) * InputMatrix(3, 3) + InputMatrix(0, 0) * InputMatrix(2, 2) * InputMatrix(3, 3);
    InvertedMatrix(1, 2) = InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(3, 0) - InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(3, 0) - InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(3, 2) + InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(3, 2) + InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(3, 3) - InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(3, 3);
    InvertedMatrix(1, 3) = -(InputMatrix(0, 3) * InputMatrix(1, 2) * InputMatrix(2, 0)) + InputMatrix(0, 2) * InputMatrix(1, 3) * InputMatrix(2, 0) + InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(2, 2) - InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 2) - InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(2, 3) + InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 3);
    
    // Third column
    InvertedMatrix(2, 0) = -(InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 0)) + InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 0) + InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 1) - InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 3) + InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 3);
    InvertedMatrix(2, 1) = InputMatrix(0, 3) * InputMatrix(2, 1) * InputMatrix(3, 0) - InputMatrix(0, 1) * InputMatrix(2, 3) * InputMatrix(3, 0) - InputMatrix(0, 3) * InputMatrix(2, 0) * InputMatrix(3, 1) + InputMatrix(0, 0) * InputMatrix(2, 3) * InputMatrix(3, 1) + InputMatrix(0, 1) * InputMatrix(2, 0) * InputMatrix(3, 3) - InputMatrix(0, 0) * InputMatrix(2, 1) * InputMatrix(3, 3);
    InvertedMatrix(2, 2) = -(InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(3, 0)) + InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(3, 0) + InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(3, 1) - InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(3, 1) - InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(3, 3) + InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(3, 3);
    InvertedMatrix(2, 3) = InputMatrix(0, 3) * InputMatrix(1, 1) * InputMatrix(2, 0) - InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 0) - InputMatrix(0, 3) * InputMatrix(1, 0) * InputMatrix(2, 1) + InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 1) + InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 3) - InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 3);
    
    // Fourth column
    InvertedMatrix(3, 0) = InputMatrix(1, 2) * InputMatrix(2, 1) * InputMatrix(3, 0) - InputMatrix(1, 1) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(1, 2) * InputMatrix(2, 0) * InputMatrix(3, 1) + InputMatrix(1, 0) * InputMatrix(2, 2) * InputMatrix(3, 1) + InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 2) - InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 2);
    InvertedMatrix(3, 1) = -(InputMatrix(0, 2) * InputMatrix(2, 1) * InputMatrix(3, 0)) + InputMatrix(0, 1) * InputMatrix(2, 2) * InputMatrix(3, 0) + InputMatrix(0, 2) * InputMatrix(2, 0) * InputMatrix(3, 1) - InputMatrix(0, 0) * InputMatrix(2, 2) * InputMatrix(3, 1) - InputMatrix(0, 1) * InputMatrix(2, 0) * InputMatrix(3, 2) + InputMatrix(0, 0) * InputMatrix(2, 1) * InputMatrix(3, 2);
    InvertedMatrix(3, 2) = InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(3, 0) - InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(3, 0) - InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(3, 1) + InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(3, 1) + InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(3, 2) - InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(3, 2);
    InvertedMatrix(3, 3) = -(InputMatrix(0, 2) * InputMatrix(1, 1) * InputMatrix(2, 0)) + InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 0) + InputMatrix(0, 2) * InputMatrix(1, 0) * InputMatrix(2, 1) - InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 1) - InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 2) + InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 2);

    // Calculation of determinant (of the input matrix)
    InputMatrixDet = InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 0) - InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 2) * InputMatrix(3, 1) + InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(0, 1) * InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 2) + InputMatrix(0, 0) * InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 2) + InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 2) - InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 2) + InputMatrix(0, 3) * (InputMatrix(1, 2) * InputMatrix(2, 1) * InputMatrix(3, 0) - InputMatrix(1, 1) * InputMatrix(2, 2) * InputMatrix(3, 0) - InputMatrix(1, 2) * InputMatrix(2, 0) * InputMatrix(3, 1) + InputMatrix(1, 0) * InputMatrix(2, 2) * InputMatrix(3, 1) + InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 2) - InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 2)) + (InputMatrix(0, 1) * InputMatrix(1, 2) * InputMatrix(2, 0) - InputMatrix(0, 0) * InputMatrix(1, 2) * InputMatrix(2, 1) - InputMatrix(0, 1) * InputMatrix(1, 0) * InputMatrix(2, 2) + InputMatrix(0, 0) * InputMatrix(1, 1) * InputMatrix(2, 2)) * InputMatrix(3, 3) + InputMatrix(0, 2) * (-(InputMatrix(1, 3) * InputMatrix(2, 1) * InputMatrix(3, 0)) + InputMatrix(1, 1) * InputMatrix(2, 3) * InputMatrix(3, 0) + InputMatrix(1, 3) * InputMatrix(2, 0) * InputMatrix(3, 1) - InputMatrix(1, 0) * InputMatrix(2, 3) * InputMatrix(3, 1) - InputMatrix(1, 1) * InputMatrix(2, 0) * InputMatrix(3, 3) + InputMatrix(1, 0) * InputMatrix(2, 1) * InputMatrix(3, 3));

    // Finalizing the calculation of the inverted matrix
    InvertedMatrix /= InputMatrixDet;
        
        KRATOS_CATCH("");
    }

    /**
     * Calculates the determinant of a matrix of dimension 2x2 or 3x3 (no check performed on matrix size)
     * @param A: Is the input matrix
     * @return The determinant of the 2x2 matrix
     */
        
    static inline TDataType Det(const MatrixType& A)
    {
        TDataType Det;
        
        if (A.size1() == 2)
        {
            Det = Det2(A);
        }
        else if (A.size1() == 3)
        {
            Det = Det3(A);
        }
        else if (A.size1() == 4)
        {
            Det = Det4(A);
        }
        else
        {
            using namespace boost::numeric::ublas;
            typedef permutation_matrix<SizeType> pmatrix;
            Matrix Aux(A);
            pmatrix pm(Aux.size1());
            bool singular = lu_factorize(Aux,pm);
            
            if (singular == true)
            {
                return 0.0;
            }
            
            Det = 1.0;
            
            for (unsigned int i = 0; i < Aux.size1();i++)
            {
                unsigned int ki = pm[i] == i ? 0 : 1;
                Det *= std::pow(-1.0, ki) * Aux(i,i);
            }
        }

        return Det;
    }
    
    /**
     * Calculates the determinant of a matrix of dimension 2x2 or 3x3 (no check performed on matrix size)
     * @param A: Is the input matrix
     * @return The determinant of the 2x2 matrix
     */
        
    static inline TDataType GeneralizedDet(const MatrixType& A)
    {
        TDataType determinant;
        
        if (A.size1() == A.size2())
        {
            determinant = Det(A);
        }
        else if (A.size1() < A.size2()) // Right determinant
        {
            Matrix AAT = prod( A, trans(A) );
            determinant = std::sqrt(Det(AAT));
        }
        else // Left determinant
        {
            Matrix ATA = prod( trans(A), A );
            determinant = std::sqrt(Det(ATA));
        }
        
        return determinant;
    }
    
    /**
     * Calculates the determinant of a matrix of dimension 2x2 (no check performed on matrix size)
     * @param A: Is the input matrix
     * @return The determinant of the 2x2 matrix
     */
    
    static inline TDataType Det2(const MatrixType& A)
    {
        return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
    }

    /**
     * Calculates the determinant of a matrix of dimension 3*3 (no check performed on matrix size)
     * @param A: Is the input matrix
     * @return The determinant of the 3x3 matrix
     */
    
    static inline TDataType Det3(const MatrixType& A)
    {
        // Calculating the algebraic complements to the first line
        const double a = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        const double b = A(1,0)*A(2,2) - A(1,2)*A(2,0);
        const double c = A(1,0)*A(2,1) - A(1,1)*A(2,0);

        return A(0,0)*a - A(0,1)*b + A(0,2)*c;
    }
    
    /**
     * Calculates the determinant of a matrix of dimension 4*4 (no check performed on matrix size)
     * @param A: Is the input matrix
     * @return The determinant of the 4x4 matrix
     */
    
    static inline TDataType Det4(const MatrixType& A)
    {
        const double Det = A(0,1)*A(1,3)*A(2,2)*A(3,0)-A(0,1)*A(1,2)*A(2,3)*A(3,0)-A(0,0)*A(1,3)*A(2,2)*A(3,1)+A(0,0)*A(1,2)*A(2,3)*A(3,1)   
                          -A(0,1)*A(1,3)*A(2,0)*A(3,2)+A(0,0)*A(1,3)*A(2,1)*A(3,2)+A(0,1)*A(1,0)*A(2,3)*A(3,2)-A(0,0)*A(1,1)*A(2,3)*A(3,2)+A(0,3)*(A(1,2)*A(2,1)*A(3,0)-A(1,1)*A(2,2)*A(3,0)-A(1,2)*A(2,0)*A(3,1)+A(1,0)*A(2,2)*A(3,1)+A(1,1)*A(2,0)*A(3,2)
                          -A(1,0)*A(2,1)*A(3,2))+(A(0,1)*A(1,2)*A(2,0)-A(0,0)*A(1,2)*A(2,1)-A(0,1)*A(1,0)*A(2,2)+A(0,0)*A(1,1)*A(2,2))*A(3,3)+A(0,2)*(-(A(1,3)*A(2,1)*A(3,0))+A(1,1)*A(2,3)*A(3,0)+A(1,3)*A(2,0)*A(3,1)-A(1,0)*A(2,3)*A(3,1)-A(1,1)*A(2,0)*A(3,3)+A(1,0)*A(2,1)*A(3,3));
        return Det;
    }

    /**
     * Calculates the determinant of a matrix of dimension 2x2 (in this case for a bounded matrix)
     * @param A: Is the input matrix
     * @return The determinant of the matrix
     */
    
    static inline TDataType Det(const boost::numeric::ublas::bounded_matrix<double,2,2>& A)
    {
        return (A(0,0)*A(1,1)-A(0,1)*A(1,0));
    }
    
    /**
     * Calculates the determinant of a matrix of dimension 3x3 (in this case for a bounded matrix)
     * @param A: Is the input matrix
     * @return The determinant of the matrix
     */
    
    static inline TDataType Det(const boost::numeric::ublas::bounded_matrix<double,3,3>& A)
    {
        // Calculating the algebraic complements to the first line
        const double a = A(1,1)*A(2,2) - A(1,2)*A(2,1);
        const double b = A(1,0)*A(2,2) - A(1,2)*A(2,0);
        const double c = A(1,0)*A(2,1) - A(1,1)*A(2,0);

        return A(0,0)*a - A(0,1)*b + A(0,2)*c;
    }

    /**
     * Performs the dot product of two vectors of dimension 3
     * (no check performed on vector sizes)
     * @param a: First input vector
     * @param b: Second input vector
     * @return The resulting norm
     */
    
    static inline TDataType Dot3(
        Vector& a, 
        Vector& b
        )
    {
        return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    }

    /**
     * Performs the dot product of two vectors of arbitrary size
     * (no check performed on vector sizes)
     * @param FirstVector: First input vector
     * @param SecondVector: Second input vector
     * @return The resulting norm
     */
    
    static inline TDataType Dot(
        const Vector& FirstVector, 
        const Vector& SecondVector
        )
    {
        Vector::const_iterator i = FirstVector.begin();
        Vector::const_iterator j = SecondVector.begin();
        TDataType temp = TDataType();
        while(i != FirstVector.end())
        {
            temp += *i++ * *j++;
        }
        return temp;
        //return std::inner_product(FirstVector.begin(), FirstVector.end(), SecondVector.begin(), TDataType());
    }

    /**
     * Calculates the norm of vector "a" which is assumed to be of size 3
     * (no check is performed on the vector's size)
     * @param a: Input vector
     * @return The resulting norm
     */
    
    static inline TDataType Norm3(Vector& a)
    {
        TDataType temp = std::pow(a[0],2) + std::pow(a[1],2) + std::pow(a[2],2);
        temp = std::sqrt(temp);
        return temp;
    }

    static inline double Norm3(const array_1d<double, 3>& a)
    {
        double temp = std::pow(a[0],2) + std::pow(a[1],2) + std::pow(a[2],2);
        temp = std::sqrt(temp);
        return temp;
    }

    /**
     * Calculates the norm of vector "a"
     * @param a: Input vector
     * @return The resulting norm
     */
    
    static inline TDataType Norm(const Vector& a)
    {
        Vector::const_iterator i = a.begin();
        TDataType temp = TDataType();
        while(i != a.end())
        {
            temp += (*i) * (*i);
            i++;
        }
        return std::sqrt(temp);
    }

    /**
     * Performs the vector product of the two input vectors a,b
     * a,b are assumed to be of size 3 (no check is performed on vector sizes)
     * @param a: First input vector
     * @param b: Second input vector
     * @return The resulting vector
     */
    
    static inline Vector CrossProduct(
        Vector& a, 
        Vector& b
        )
    {
        Vector c(3);

        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];

        return c;
    }

    //this function is deprecated since instead of giving back vec x Tuple it gives back Tuple x Vec ( = -Vec x Tuple)
    //THAT IS -- THIS FUNCTION GIVES BACK THE OPPOSITE SIGN OF THE PRODUCT
    //please use instead CrossProd(c,a,b) which is also in general more optimal since it never creates tmp on return
    KRATOS_DEPRECATED static inline array_1d<double, 3> UnitCrossProduct(
        const array_1d<double, 3>& vec, 
        const array_1d<double, 3>& Tuple
        )
    {
        array_1d<double, 3> cross;

        cross[0] =  Tuple[1]*vec[2] - Tuple[2]*vec[1];
        cross[1] =  Tuple[2]*vec[0] - Tuple[0]*vec[2];
        cross[2] =  Tuple[0]*vec[1] - Tuple[1]*vec[0];

        const double length = std::sqrt(inner_prod(cross, cross));
        cross = (1.00/length) * cross;
        return cross;
    }

    //this function is deprecated since instead of giving back vec x Tuple it gives back Tuple x Vec ( = -Vec x Tuple)
    //THAT IS -- THIS FUNCTION GIVES BACK THE OPPOSITE SIGN OF THE PRODUCT
    //please use instead CrossProd(c,a,b) which is also in general more optimal since it never creates tmp on return
    KRATOS_DEPRECATED static inline array_1d<double, 3> CrossProduct(
        const array_1d<double, 3>& vec, 
        const array_1d<double, 3>& Tuple
        )
    {
        array_1d<double, 3> cross;

        cross[0] =  Tuple[1]*vec[2] - Tuple[2]*vec[1];
        cross[1] =  Tuple[2]*vec[0] - Tuple[0]*vec[2];
        cross[2] =  Tuple[0]*vec[1] - Tuple[1]*vec[0];
 
	return cross;
    }
    
    /**
     * Identical but it assumes that the output vector is given already the correct size
     * @param a: First input vector
     * @param b: Second input vector
     * @return c: The resulting vector
     */
    
    static inline void CrossProduct(
        array_1d<double, 3>& c, 
        const array_1d<double, 3>& a, 
        const array_1d<double, 3>& b
        )
    {
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }
    
    //general version    
    template< class T1, class T2, class T3 >
    static inline void CrossProduct(T1& c, const T2& a, const T3& b ){
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];
    }

    template< class T1, class T2, class T3 >
    static inline void UnitCrossProduct(T1& c, const T2& a, const T3& b ){
        CrossProduct(c,a,b);
        double norm = norm_2(c);
#ifdef KRATOS_DEBUG
        if(norm < 1000.0*std::numeric_limits<double>::epsilon())
            KRATOS_ERROR << "norm is 0 when making the UnitCrossProduct of the vectors " << a << " and " << b << std::endl;
#endif
        c/=norm;
    }
    
    /**
     * Returns a matrix :
     * A = a.tensorproduct.b
     * a,b are assumed to be of order 3, no check is performed on the size of the vectors
     * @param a: First input vector
     * @param b: Second input vector
     */
    
    static inline MatrixType TensorProduct3(
        Vector& a, 
        Vector& b
        )
    {
        MatrixType A(3,3);
        
        A(0,0) = a[0]*b[0];
        A(0,1) = a[0]*b[1];
        A(0,2) = a[0]*b[2];
        A(1,0) = a[1]*b[0];
        A(1,1) = a[1]*b[1];
        A(1,2) = a[1]*b[2];
        A(2,0) = a[2]*b[0];
        A(2,1) = a[2]*b[1];
        A(2,2) = a[2]*b[2];
        
        return A;
    }
    
    /**
     * "InputMatrix" is ADDED to "Destination" matrix starting from
     * InitialRow and InitialCol of the destination matrix
     * "Destination" is assumed to be able to contain the "input matrix"
     * (no check is performed on the bounds)
     * @return Destination: The matric destination
     * @param InputMatrix: The input matrix to be computed
     * @param InitialRow: The initial row to compute
     * @param InitialCol: The initial column to compute
     */
    
    static inline void  AddMatrix(
        MatrixType& Destination,
        MatrixType& InputMatrix,
        int InitialRow,
        int InitialCol
        )
    {
        KRATOS_TRY
        for(unsigned int i = 0; i < InputMatrix.size1(); i++)
        {
            for(unsigned int j = 0; j < InputMatrix.size2(); j++)
            {
                Destination(InitialRow+i, InitialCol+j) += InputMatrix(i,j);
            }
        }
        KRATOS_CATCH("")
    }
    
    /**
     *  "InputMatrix" is SUBTRACTED to "Destination" matrix starting from
     * InitialRow and InitialCol of the destination matrix
     * "Destination" is assumed to be able to contain the "input matrix"
     * (no check is performed on the bounds)
     * @return Destination: The matric destination
     * @param InputMatrix: The input matrix to be computed
     * @param InitialRow: The initial row to compute
     * @param InitialCol: The initial column to compute
     */
        
    static inline void  SubtractMatrix(
        MatrixType& Destination,
        MatrixType& InputMatrix,
        int InitialRow,
        int InitialCol
        )
    {
        KRATOS_TRY;
        
        for(unsigned int i = 0; i<InputMatrix.size1(); i++)
        {
            for(unsigned int j = 0; j<InputMatrix.size2(); j++)
            {
                Destination(InitialRow+i, InitialCol+j) -= InputMatrix(i,j);
            }
        }
        
        KRATOS_CATCH("");
    }

    /**
     * "InputMatrix" is WRITTEN on "Destination" matrix starting from
     * InitialRow and InitialCol of the destination matrix
     * "Destination" is assumed to be able to contain the "input matrix"
     * (no check is performed on the bounds)
     * ATTENTION: Destination is overwritten!!
     * @return Destination: The matric destination
     * @param InputMatrix: The input matrix to be computed
     * @param InitialRow: The initial row to compute
     * @param InitialCol: The initial column to compute
     */
        
    static inline void  WriteMatrix(
        MatrixType& Destination,
        MatrixType& InputMatrix,
        int InitialRow,
        int InitialCol
        )
    {
        KRATOS_TRY;
        
        for(unsigned int i = 0; i < InputMatrix.size1(); i++)
        {
            for(unsigned int j = 0; j < InputMatrix.size2(); j++)
            {
                Destination(InitialRow+i, InitialCol+j) = InputMatrix(i,j);
            }
        }
            
        KRATOS_CATCH("");
    }

    /**
     * Performs the Kroneker product of the Reduced Matrix with the identity matrix of size "dimension"
     * @return Destination: The matric destination
     * @param ReducedMatrix: The reduced matrix to be computed
     * @param dimension: The dimension where we work
     */
    
    static inline void  ExpandReducedMatrix(
        MatrixType& Destination,
        MatrixType& ReducedMatrix,
        unsigned int dimension
        )
    {
        KRATOS_TRY;
        
        const unsigned int size = ReducedMatrix.size2();

        for (unsigned int i = 0; i < size; i++)
        {
            unsigned int rowindex = i*dimension;
            for (unsigned int j = 0; j < size; j++)
            {
                unsigned int colindex = j*dimension;
                for(unsigned int ii=0; ii<dimension; ii++)
                {
                    Destination(rowindex+ii,colindex+ii) = ReducedMatrix(i,j);
                }
            }
        }
        
        KRATOS_CATCH("");
    }

    /**
     * Performs the Kroneker product of the Reduced Matrix with the identity matrix of size "dimension" ADDING to the destination matrix
     * @return Destination: The matric destination
     * @param ReducedMatrix: The reduced matrix to be added
     * @param dimension: The dimension where we work
     */
        
    static inline void  ExpandAndAddReducedMatrix(
        MatrixType& Destination,
        MatrixType& ReducedMatrix,
        const unsigned int dimension
        )
    {
        KRATOS_TRY;
        
        const unsigned int size = ReducedMatrix.size2();
        unsigned int rowindex = 0;
	unsigned int colindex = 0;
        
        for (unsigned int i = 0; i < size; i++)
        {
            rowindex = i * dimension;
            for (unsigned int j = 0; j < size; j++)
            {
                colindex = j * dimension;
                for(unsigned int ii = 0; ii < dimension; ii++)
                {
                    Destination(rowindex+ii,colindex+ii) += ReducedMatrix(i,j);
                }
            } 
        }
        
        KRATOS_CATCH("");
    }

    /**
     * Performs x += coeff*y. no check on bounds is performed
     * @return x: The vector destination
     * @param y: The vector to be added
     * @param coeff: The proportion to be added
     */
        
    static inline void  VecAdd(
        Vector& x,
        TDataType coeff,
        Vector& y)
    {
        KRATOS_TRY
        unsigned int size=x.size();

        for (unsigned int i=0; i<size; i++)
        {
            x[i] += coeff * y[i];
        }
        KRATOS_CATCH("")
    }

   /**
     * Transforms a stess vector into a matrix. Stresses are assumed to be stored
     * in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * @param rStressVector the given stress vector
     * @return the corresponding stress tensor in matrix form
     */
   
    static inline MatrixType StressVectorToTensor(const Vector& rStressVector)
    {
        KRATOS_TRY
        Matrix StressTensor;

        if (rStressVector.size()==3)
        {
            StressTensor.resize(2,2,false);
            StressTensor(0,0) = rStressVector[0];
            StressTensor(0,1) = rStressVector[2];
            StressTensor(1,0) = rStressVector[2];
            StressTensor(1,1) = rStressVector[1];
        }
        else if (rStressVector.size()==4)
        {
            StressTensor.resize(3,3,false);
            StressTensor(0,0) = rStressVector[0];
            StressTensor(0,1) = rStressVector[3];
            StressTensor(0,2) = 0.0;
            StressTensor(1,0) = rStressVector[3];
            StressTensor(1,1) = rStressVector[1];
            StressTensor(1,2) = 0.0;
            StressTensor(2,0) = 0.0;
            StressTensor(2,1) = 0.0;
            StressTensor(2,2) = rStressVector[2];
        }
        else if (rStressVector.size()==6)
        {
            StressTensor.resize(3,3,false);
            StressTensor(0,0) = rStressVector[0];
            StressTensor(0,1) = rStressVector[3];
            StressTensor(0,2) = rStressVector[5];
            StressTensor(1,0) = rStressVector[3];
            StressTensor(1,1) = rStressVector[1];
            StressTensor(1,2) = rStressVector[4];
            StressTensor(2,0) = rStressVector[5];
            StressTensor(2,1) = rStressVector[4];
            StressTensor(2,2) = rStressVector[2];
        }

        return StressTensor;
        
        KRATOS_CATCH("");
    }

   /**
     * Transforms a  vector into a symmetric matrix.  Components are assumed to be stored
     * in the following way:
     * \f$ [ s11, s22, s33, s12, s23, s13 ] \f$ for 3D case and
     * \f$ [ s11, s22, s33, s12 ] \f$ for 2D case.
     * \f$ [ s11, s22, s12 ] \f$ for 2D case.
     * @param rVector the given stress vector
     * @return the corresponding Tensor in matrix form
     */
   
    static inline MatrixType VectorToSymmetricTensor(const Vector& rVector)
    {
        KRATOS_TRY;
        
        Matrix Tensor;

        if (rVector.size() == 3)
        {
            Tensor.resize(2,2,false);
            Tensor(0,0) = rVector[0];
            Tensor(0,1) = rVector[2];
            Tensor(1,0) = rVector[2];
            Tensor(1,1) = rVector[1];
        }
        else if (rVector.size() == 4)
        {
            Tensor.resize(3,3,false);
            Tensor(0,0) = rVector[0];
            Tensor(0,1) = rVector[3];
            Tensor(0,2) = 0.0;
            Tensor(1,0) = rVector[3];
            Tensor(1,1) = rVector[1];
            Tensor(1,2) = 0.0;
            Tensor(2,0) = 0.0;
            Tensor(2,1) = 0.0;
            Tensor(2,2) = rVector[2];
        }
        else if (rVector.size() == 6)
        {
            Tensor.resize(3,3,false);
            Tensor(0,0) = rVector[0];
            Tensor(0,1) = rVector[3];
            Tensor(0,2) = rVector[5];
            Tensor(1,0) = rVector[3];
            Tensor(1,1) = rVector[1];
            Tensor(1,2) = rVector[4];
            Tensor(2,0) = rVector[5];
            Tensor(2,1) = rVector[4];
            Tensor(2,2) = rVector[2];
        }

        return Tensor;
        
        KRATOS_CATCH("");
    }
    
    /**
     * Sign function
     * @param ThisDataType: The value to extract the sign
     * @return The sign of the value
     */
    static inline int Sign(const TDataType& ThisDataType)
    {
        KRATOS_TRY;
        const TDataType& x = ThisDataType;
        return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
        KRATOS_CATCH("");
    }
    

    /**
     * Transforms a strain vector into a matrix. Strains are assumed to be stored
     * in the following way:
     * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     * \f$ [ e11, e22, e33, 2*e12 ] \f$ for 2D case.
     * \f$ [ e11, e22, 2*e12 ] \f$ for 2D case.
     * Hence the deviatoric components of the strain vector are divided by 2
     * while they are stored into the matrix
     * @param rStrainVector the given strain vector
     * @return the corresponding strain tensor in matrix form
     */
    
    static inline MatrixType StrainVectorToTensor( const VectorType& rStrainVector)
    {
        KRATOS_TRY
        Matrix StrainTensor;

        if (rStrainVector.size()==3)
        {
            StrainTensor.resize(2,2, false);
            StrainTensor(0,0) = rStrainVector[0];
            StrainTensor(0,1) = 0.5*rStrainVector[2];
            StrainTensor(1,0) = 0.5*rStrainVector[2];
            StrainTensor(1,1) = rStrainVector[1];
        }
        else if (rStrainVector.size()==4)
        {
            StrainTensor.resize(3,3, false);
            StrainTensor(0,0) = rStrainVector[0];
            StrainTensor(0,1) = 0.5*rStrainVector[3];
            StrainTensor(0,2) = 0;
            StrainTensor(1,0) = 0.5*rStrainVector[3];
            StrainTensor(1,1) = rStrainVector[1];
            StrainTensor(1,2) = 0;
            StrainTensor(2,0) = 0;
            StrainTensor(2,1) = 0;
            StrainTensor(2,2) = rStrainVector[2];
        }
        else if (rStrainVector.size()==6)
        {
            StrainTensor.resize(3,3, false);
            StrainTensor(0,0) = rStrainVector[0];
            StrainTensor(0,1) = 0.5*rStrainVector[3];
            StrainTensor(0,2) = 0.5*rStrainVector[5];
            StrainTensor(1,0) = 0.5*rStrainVector[3];
            StrainTensor(1,1) = rStrainVector[1];
            StrainTensor(1,2) = 0.5*rStrainVector[4];
            StrainTensor(2,0) = 0.5*rStrainVector[5];
            StrainTensor(2,1) = 0.5*rStrainVector[4];
            StrainTensor(2,2) = rStrainVector[2];
        }

        return StrainTensor;
        
        KRATOS_CATCH("");
    }

    /**
     * Transforms a given symmetric Strain Tensor to Voigt Notation:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     * in the 2D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * \f$ [ e11, e22, e33, 2*e12 ] \f$ fir 2D case.
     * in the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (3*1) Vector
     * \f$ [ e11, e22, 2*e12 ] \f$ fir 2D case.
     * @param rStrainTensor the given symmetric second order strain tensor
     * @return the corresponding strain tensor in vector form
     */

    static inline Vector StrainTensorToVector( 
        const Matrix& rStrainTensor, 
        unsigned int rSize = 0 
        )
    {
        KRATOS_TRY;

        Vector StrainVector;

        if(rSize == 0)
        {
            if(rStrainTensor.size1() == 2)
            {
                rSize = 3;
            }
            else if(rStrainTensor.size1() == 3)
            {
                rSize = 6;
            }
        }

        if (rSize == 3)
        {
            StrainVector.resize(3,false);
            StrainVector[0] = rStrainTensor(0,0);
            StrainVector[1] = rStrainTensor(1,1);
            StrainVector[2] = 2.00*rStrainTensor(0,1);
        }
        else if (rSize == 4)
        {
            StrainVector.resize(4,false);
            StrainVector[0] = rStrainTensor(0,0);
            StrainVector[1] = rStrainTensor(1,1);
            StrainVector[2] = rStrainTensor(2,2);
            StrainVector[3] = 2.00*rStrainTensor(0,1);
        }
        else if (rSize == 6)
        {
            StrainVector.resize(6,false);
            StrainVector[0] = rStrainTensor(0,0);
            StrainVector[1] = rStrainTensor(1,1);
            StrainVector[2] = rStrainTensor(2,2);
            StrainVector[3] = 2.00*rStrainTensor(0,1);
            StrainVector[4] = 2.00*rStrainTensor(1,2);
            StrainVector[5] = 2.00*rStrainTensor(0,2);
        }

        return StrainVector;
        
        KRATOS_CATCH("");
     }

    /**
     * Transforms a given symmetric Stress Tensor to Voigt Notation:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (3*1) Vector
     * @param rStressTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */
    
    static inline Vector StressTensorToVector(
        const Matrix& rStressTensor, 
        unsigned int rSize = 0
        )
    {
        KRATOS_TRY;
        
        Vector StressVector;
      
        if(rSize == 0)
        {
            if(rStressTensor.size1() == 2)
            {
                rSize = 3;
            }
            else if(rStressTensor.size1() == 3)
            {
                rSize = 6;
            }
        }

        if (rSize == 3)
        {
            StressVector.resize(3,false);
            StressVector[0] = rStressTensor(0,0);
            StressVector[1] = rStressTensor(1,1);
            StressVector[2] = rStressTensor(0,1);
        }
        else if (rSize == 4)
        {
            StressVector.resize(4,false);
            StressVector[0] = rStressTensor(0,0);
            StressVector[1] = rStressTensor(1,1);
            StressVector[2] = rStressTensor(2,2);
            StressVector[3] = rStressTensor(0,1);
            
        }
        else if (rSize == 6)
        {
            StressVector.resize(6,false);
            StressVector[0] = rStressTensor(0,0);
            StressVector[1] = rStressTensor(1,1);
            StressVector[2] = rStressTensor(2,2);
            StressVector[3] = rStressTensor(0,1);
            StressVector[4] = rStressTensor(1,2);
            StressVector[5] = rStressTensor(0,2);
        }

        return StressVector;
        
        KRATOS_CATCH("");
     }

    /**
     * Transforms a given symmetric Tensor to Voigt Notation:
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (6*1) Vector
     * in the 3D case: from a second order tensor (3*3) Matrix  to a corresponing (4*1) Vector
     * in the 2D case: from a second order tensor (2*2) Matrix  to a corresponing (3*1) Vector
     * @param rTensor the given symmetric second order stress tensor
     * @return the corresponding stress tensor in vector form
     */
    
    static inline Vector SymmetricTensorToVector(
        const Matrix& rTensor, 
        unsigned int rSize = 0
        )
    {
        KRATOS_TRY;
        
        Vector vector;
        
        if(rSize == 0)
        {
            if(rTensor.size1() == 2)
            {
                rSize = 3;
            }
            else if(rTensor.size1() == 3)
            {
                rSize = 6;
            }
        }
        
        if (rSize == 3)
        {
            vector.resize(3,false);
            vector[0]= rTensor(0,0);
            vector[1]= rTensor(1,1);
            vector[2]= rTensor(0,1);
            
        }
        else if (rSize==4)
        {
            vector.resize(4,false);
            vector[0]= rTensor(0,0);
            vector[1]= rTensor(1,1);
            vector[2]= rTensor(2,2);
            vector[3]= rTensor(0,1);
        }
        else if (rSize==6)
        {
            vector.resize(6);
            vector[0]= rTensor(0,0);
            vector[1]= rTensor(1,1);
            vector[2]= rTensor(2,2);
            vector[3]= rTensor(0,1);
            vector[4]= rTensor(1,2);
            vector[5]= rTensor(0,2);
        }
        
        return vector;
        
        KRATOS_CATCH("");
     }

    /**
     * Calculates the eigenvectors and eigenvalues of given symmetric TDimxTDim matrix
     * The eigenvectors and eigenvalues are calculated using the iterative Gauss-Seidel-method
     * @param A: The given symmetric matrix the eigenvectors are to be calculated.
     * @return eigen_vector_matrix: The result matrix (will be overwritten with the eigenvectors)
     * @return eigen_values_matrix: The result diagonal matrix with the eigenvalues
     * @param tolerance: The largest value considered to be zero
     * @param max_iterations: Maximum number of iterations
     */

    template<unsigned int TDim>  
    static inline bool EigenSystem(
            const boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim>& A,
            boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim>& eigen_vector_matrix,
            boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim>& eigen_values_matrix,
            const TDataType tolerance = 1.0e-18,
            const unsigned int max_iterations = 20
            )
    {
        bool is_converged = false;
        eigen_values_matrix = ZeroMatrix(TDim);
        boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> TempMat = A;
        boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> AuxA;

        const boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> Indentity = IdentityMatrix(TDim, TDim);
        boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> V = Indentity;
        boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> Vaux;
        boost::numeric::ublas::bounded_matrix<TDataType, TDim, TDim> Rotation;

        for(unsigned int iterations = 0; iterations < max_iterations; iterations++)
        {
            is_converged = true;

            TDataType a = 0.0;
            unsigned int index1 = 0;
            unsigned int index2 = 1;

            for(unsigned int i = 0; i < TDim; i++)
            {
                for(unsigned int j = (i + 1); j < TDim; j++)
                {
                    if((std::abs(TempMat(i, j)) > a ) && (std::abs(TempMat(i, j)) > tolerance))
                    {
                        a = std::abs(TempMat(i,j));
                        index1 = i;
                        index2 = j;
                        is_converged = false;
                    }
                }
            }

            if(is_converged)
            {
                break;
            }

            // Calculation of Rotation angle
            TDataType gamma = (TempMat(index2, index2)-TempMat(index1, index1)) / (2 * TempMat(index1, index2));
            TDataType u = 1.0;

            if(std::abs(gamma) > tolerance && std::abs(gamma)< (1/tolerance))
            {
                u = gamma / std::abs(gamma) * 1.0 / (std::abs(gamma) + std::sqrt(1.0 + gamma * gamma));
            }
            else
            {
                if  (std::abs(gamma) >= (1.0/tolerance))
                {
                    u = 0.5 / gamma;
                }
            }

            TDataType c = 1.0 / (std::sqrt(1.0 + u * u));
            TDataType s = c * u;
            TDataType teta = s / (1.0 + c);

            // Rotation of the Matrix
            AuxA = TempMat;
            AuxA(index2, index2) = TempMat(index2,index2) + u * TempMat(index1, index2);
            AuxA(index1, index1) = TempMat(index1,index1) - u * TempMat(index1, index2);
            AuxA(index1, index2) = 0.0;
            AuxA(index2, index1) = 0.0;

            for(unsigned int i = 0; i < TDim; i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    AuxA(index2, i) = TempMat(index2, i) + s * (TempMat(index1, i)- teta * TempMat(index2, i));
                    AuxA(i, index2) = TempMat(index2, i) + s * (TempMat(index1, i)- teta * TempMat(index2, i));
                    AuxA(index1, i) = TempMat(index1, i) - s * (TempMat(index2, i) + teta * TempMat(index1, i));
                    AuxA(i, index1) = TempMat(index1, i) - s * (TempMat(index2, i) + teta * TempMat(index1, i));
                }
            }

            TempMat = AuxA;

            // Calculation of the eigeneigen_vector_matrix V
            Rotation = Indentity;
            Rotation(index2, index1) = -s;
            Rotation(index1, index2) =  s;
            Rotation(index1, index1) =  c;
            Rotation(index2, index2) =  c;

            Vaux = ZeroMatrix(TDim, TDim);

            for(unsigned int i = 0; i < TDim; i++)
            {
                for(unsigned int j = 0; j < TDim; j++)
                {
                    for(unsigned int k = 0; k < TDim; k++)
                    {
                        Vaux(i, j) += V(i, k) * Rotation(k, j);
                    }
                }
            }
            V = Vaux;
        }

        if(!(is_converged))
        {
            std::cout<<" WARNING: Spectral decomposition not converged "<<std::endl;
        }

        for(unsigned int i = 0; i < TDim; i++)
        {
            eigen_values_matrix(i, i) = TempMat(i, i);
            for(unsigned int j = 0; j < TDim; j++)
            {
                eigen_vector_matrix(i, j) = V(j, i);
            }
        }
    
        return is_converged;
    }
    
    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

private:
    
    ///@name Private static Member Variables
    ///@{

    ///@}
    ///@name Private member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Private LifeCycle
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    MathUtils(void);

    MathUtils(MathUtils& rSource);

}; /* Class MathUtils */

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}  /* namespace Kratos.*/

#endif /* KRATOS_MATH_UTILS  defined */

