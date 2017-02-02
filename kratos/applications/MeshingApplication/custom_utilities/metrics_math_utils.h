// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_METRICS_MATH_UTILS)
#define KRATOS_METRICS_MATH_UTILS

// Project includes
#include "utilities/math_utils.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"

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

//// This class is used to compute some mathematical operations needed for the metrics computing

template<unsigned int TDim>  
class MetricsMathUtils
{
public:

    ///@name Type Definitions
    ///@{
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
 
    /**
     * Calculates the determinant of a 2x2 or 3x3 matrix 
     * @param A: The matrix to calculate
     * @return DetA: The determinant of the matrix
     */

    static inline double DetMat(const boost::numeric::ublas::bounded_matrix<double, TDim, TDim>& A)
    {
        double det;
        
        if (TDim == 2)
        {
            det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
        }
        else
        {
            det = A(0, 0) * A(1, 1) * A(2, 2)
                + A(1, 0) * A(2, 1) * A(0, 2)
                + A(0, 1) * A(1, 2) * A(2, 0)
                - A(2, 0) * A(1, 1) * A(0, 2)
                - A(2, 1) * A(1, 2) * A(0, 0)
                - A(1, 0) * A(0, 1) * A(2,2);
        }
        
        return det;
    }
    
    /**
     * Calculates the inverse of a 2x2 or 3x3 matrix 
     * @param A: The matrix to invert
     * @return InvA: The inverted matrix
     */

    static inline boost::numeric::ublas::bounded_matrix<double, TDim, TDim> InvMat(
            const boost::numeric::ublas::bounded_matrix<double, TDim, TDim>& A,
            const double tolerance = 1.0e-18
            )
    {
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> InvA;
        
        /* Compute determinant of the matrix */
        const double det = DetMat(A);
        
        if (std::abs(det) < tolerance)
        {
            KRATOS_WATCH(A);
            KRATOS_THROW_ERROR( std::invalid_argument," Determinant of the matrix is zero or almost zero!!!, det = ", det);
        }
        
        if (TDim == 2)
        {            
            /* Compute inverse of the Matrix */
            InvA(0, 0) =   A(1, 1) / det;
            InvA(0, 1) = - A(0, 1) / det;
            InvA(1, 0) = - A(1, 0) / det;
            InvA(1, 1) =   A(0, 0) / det;
        }
        else
        {
            /* Compute inverse of the Matrix */
            InvA(0, 0) =   (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) / det;
            InvA(1, 0) = - (A(1, 0) * A(2, 2) - A(2, 0) * A(1, 2)) / det;
            InvA(2, 0) =   (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)) / det;

            InvA(0, 1) = - (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) / det;
            InvA(1, 1) =   (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) / det;
            InvA(2, 1) = - (A(0, 0) * A(2, 1) - A(0, 1) * A(2, 0)) / det;

            InvA(0, 2) =   (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) / det;
            InvA(1, 2) = - (A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0)) / det;
            InvA(2, 2) =   (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) / det;
        }
        
        return InvA;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Calculates the eigenvectors and eigenvalues of given symmetric TDimxTDim matrix
     * The eigenvectors and eigenvalues are calculated using the iterative Gauss-Seidel-method
     * @param A: The given symmetric matrix the eigenvectors are to be calculated.
     * @return eigen_vector_matrix: The result matrix (will be overwritten with the eigenvectors)
     * @return eigen_values_matrix: The result diagonal matrix with the eigenvalues
     * @param tolerance: The largest value considered to be zero
     * @param max_iterations: Maximum number of iterations
     */

    static inline bool EigenSystem(
            const boost::numeric::ublas::bounded_matrix<double, TDim, TDim>& A,
            boost::numeric::ublas::bounded_matrix<double, TDim, TDim>& eigen_vector_matrix,
            boost::numeric::ublas::bounded_matrix<double, TDim, TDim>& eigen_values_matrix,
            const double tolerance = 1.0e-18,
            const unsigned int max_iterations = 20
            )
    {
        bool is_converged = false;
        eigen_values_matrix = ZeroMatrix(TDim);
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> TempMat = A;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> AuxA;

        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Indentity = IdentityMatrix(TDim, TDim);
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> V = Indentity;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Vaux;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Rotation;

        for(unsigned int iterations = 0; iterations < max_iterations; iterations++)
        {
            is_converged = true;

            double a = 0.0;
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
            double gamma = (TempMat(index2, index2)-TempMat(index1, index1)) / (2 * TempMat(index1, index2));
            double u = 1.0;

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

            double c = 1.0 / (std::sqrt(1.0 + u * u));
            double s = c * u;
            double teta = s / (1.0 + c);

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
    
    /***********************************************************************************/
    /***********************************************************************************/
       
    /**
     * This function converts tensors to vector and viceversa
     */
    
    static inline Vector TensorToVector(const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Tensor)
    {
        Vector TensorCond;
        TensorCond.resize(3 * (TDim - 1), false);
        
        // Common part
        TensorCond[0] = Tensor(0, 0);
        TensorCond[1] = Tensor(0, 1);
            
        if (TDim == 2)
        {
            TensorCond[2] = Tensor(1, 1);
        }
        else
        {
            TensorCond[2] = Tensor(0, 2);
            TensorCond[3] = Tensor(1, 1);
            TensorCond[4] = Tensor(1, 2);
            TensorCond[5] = Tensor(2, 2);
        }

        return TensorCond;
    }
    
    static inline boost::numeric::ublas::bounded_matrix<double, TDim, TDim>  VectorToTensor(const Vector TensorCond)
    {
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Tensor;
        
        // Common part
        Tensor(0, 0) = TensorCond[0];
        Tensor(0, 1) = TensorCond[1];
        Tensor(1, 0) = TensorCond[1];
        if (TDim == 2)
        {
            Tensor(1, 1) = TensorCond[2];
        }
        else
        {
            Tensor(0, 2) = TensorCond[2];
            Tensor(2, 0) = TensorCond[2];
            Tensor(1, 1) = TensorCond[3];
            Tensor(1, 2) = TensorCond[4];
            Tensor(2, 1) = TensorCond[4];
            Tensor(2, 2) = TensorCond[5];
        }

        return Tensor;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * It computes the intersection between two metrics 
     * @param Metric1: The first metric
     * @param Metric2: The second metric
     */
        
    static inline Vector IntersectMetrics(
        const Vector Metric1,
        const Vector Metric2
    )
    {
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Metric1Matrix = VectorToTensor(Metric1);
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Metric2Matrix = VectorToTensor(Metric2);
        
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> auxmat;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> emat;
        
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> invMetric1Matrix = InvMat(Metric1Matrix);
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> NMatrix = prod(invMetric1Matrix, Metric2Matrix);
        
        EigenSystem(NMatrix, emat, auxmat, 1e-18, 20);
        
        typedef boost::numeric::ublas::bounded_matrix<double, TDim, TDim> temp_type;
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> lambdamat =  prod(trans(emat), prod<temp_type>(Metric1Matrix, emat));
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> mumat =  prod(trans(emat), prod<temp_type>(Metric2Matrix, emat));
        
        for (unsigned int i = 0; i < TDim; i++)
        {
            auxmat(i, i) = MathUtils<double>::Max(lambdamat(i, i), mumat(i, i));
        }
        
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> invemat = InvMat(emat);
        
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> IntersectionMatrix =  prod(trans(invemat), prod<temp_type>(auxmat, invemat));
        
        const Vector Intersection = TensorToVector(IntersectionMatrix);
        
        return Intersection;
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
    
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
};// class MetricsMathUtils
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

}// namespace Kratos.
#endif /* KRATOS_METRICS_MATH_UTILS defined */
