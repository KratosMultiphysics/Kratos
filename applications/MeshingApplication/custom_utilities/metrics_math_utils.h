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
            const double tolerance = 1.0e-16,
            const unsigned int max_iterations = 10
            )
    {
        bool is_converged = true;
        if (TDim == 2)
        {        
            // Taking the coefficients of the matrix
            const double a11 = A(0, 0);
            const double a12 = A(0, 1) + tolerance; // NOTE: To avoid problems with diagonal or almost diagonal matrices
            const double a22 = A(1, 1);
            
            // The eigen values matrix 
            eigen_values_matrix(0, 0) = 0.5 * (a11 + a22 - std::sqrt(a11 * a11 + 4.0 * a12 * a12 - 2.0 * a11 * a22 + a22 * a22));
            eigen_values_matrix(0, 1) = 0.0; 
            eigen_values_matrix(1, 0) = 0.0; 
            eigen_values_matrix(1, 1) = 0.5 * (a11 + a22 + std::sqrt(a11 * a11 + 4.0 * a12 * a12 - 2.0 * a11 * a22 + a22 * a22));
            
            // The eigen vector matrix
            array_1d<double, 2> aux_array;
            
            // First eigen vector
            aux_array[0] = -(-a11 + a22 + std::sqrt(a11 * a11 + 4.0 * a12 * a12 - 2.0 * a11 * a22 + a22 * a22))/(2.0 * a12);
            aux_array[1] = 1.0;
            aux_array /= norm_2(aux_array);
            eigen_vector_matrix(0, 0) = aux_array[0];
            eigen_vector_matrix(0, 1) = aux_array[1];
            
            // Second eigen vector
            aux_array[0] = -(-a11 + a22 - std::sqrt(a11 * a11 + 4.0 * a12 * a12 - 2.0 * a11 * a22 + a22 * a22))/(2.0 * a12);
            aux_array[1] = 1.0;
            aux_array /= norm_2(aux_array);
            eigen_vector_matrix(1, 0) = aux_array[0];
            eigen_vector_matrix(1, 1) = aux_array[1];
        }
        else
        {
            is_converged = false;
            boost::numeric::ublas::bounded_matrix<double, 3, 3> auxA;
            const boost::numeric::ublas::bounded_matrix<double, 3, 3> Identity = IdentityMatrix(3, 3);
            boost::numeric::ublas::bounded_matrix<double, 3, 3> V = Identity;
            boost::numeric::ublas::bounded_matrix<double, 3, 3> Rotation;

            for(unsigned int iterations = 0; iterations < max_iterations; iterations++)
            {
                is_converged = true;

                double a = 0.0;
                unsigned int index1 = 0;
                unsigned int index2 = 0;

                for(unsigned int i = 0; i < 3; i++)
                {
                    for(unsigned int j = (i + 1); j < 3; j++)
                    {
                        if((std::abs(A(i, j)) > a ) && (std::abs(A(i, j)) > tolerance))
                        {
                            a = std::abs(A(i,j));
                            index1 = i;
                            index2 = j;
                            is_converged = false;
                        }
                    }
                }

                if(is_converged == true)
                {
                    break;
                }

                // Calculation of Rotation angle
                const double gamma = (A(index2, index2)-A(index1, index1)) / (2.0 * A(index1, index2));
                double u = 1.0;

                if(std::abs(gamma) > tolerance && std::abs(gamma)< (1.0/tolerance))
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

                const double c = 1.0 / (sqrt(1.0 + u * u));
                const double s = c * u;
                const double teta = s / (1.0 + c);

                // Rotation of the Matrix
                auxA = A;
                auxA(index2, index2) = A(index2,index2) + u * A(index1, index2);
                auxA(index1, index1) = A(index1,index1) - u * A(index1, index2);
                auxA(index1, index2) = 0.0;
                auxA(index2, index1) = 0.0;

                for(unsigned int i = 0; i < 3; i++)
                {
                    if((i!= index1) && (i!= index2))
                    {
                        auxA(index2, i) = A(index2, i) + s * (A(index1, i) - teta * A(index2, i));
                        auxA(i, index2) = A(index2, i) + s * (A(index1, i) - teta * A(index2, i));
                        auxA(index1, i) = A(index1, i) - s * (A(index2, i) + teta * A(index1, i));
                        auxA(i, index1) = A(index1, i) - s * (A(index2, i) + teta * A(index1, i));
                    }
                }

                // Calculation of the eigenvectors V
                Rotation = Identity;
                Rotation(index2, index1) = -s;
                Rotation(index1, index2) =  s;
                Rotation(index1, index1) =  c;
                Rotation(index2, index2) =  c;

                eigen_vector_matrix = ZeroMatrix(3, 3);
                for(unsigned int i = 0; i < 3; i++)
                {
                    for(unsigned int j = 0; j < 3; j++)
                    {
                        for(unsigned int k = 0; k < 3; k++)
                        {
                            eigen_vector_matrix(j, i) += V(i, k) * Rotation(k, j);
                        }
                    }
                }
            }

            if(!(is_converged))
            {
                std::cout<<" WARNING: Spectral decomposition not converged "<<std::endl;
            }

            for(unsigned int i = 0; i < 3; i++)
            {
                for(unsigned int j = 0; j < 3; j++)
                {
                    if (i == j)
                    {
                        eigen_values_matrix(i, i) = auxA(i, i);
                    }
                    else
                    {
                        eigen_values_matrix(i, j) = 0.0;
                    }
                }
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
        Vector Metric1,
        Vector Metric2
    )
    {
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Metric1Matrix = VectorToTensor(Metric1);
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> Metric2Matrix = VectorToTensor(Metric2);
        
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> auxmat;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> emat;
        
        double det;
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> invMetric1Matrix;
        MathUtils<double>::InvertMatrix( Metric1Matrix, invMetric1Matrix, det);
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> NMatrix = prod(invMetric1Matrix, Metric2Matrix);
        
        EigenSystem(NMatrix, emat, auxmat, 1e-18, 10);
        
        typedef boost::numeric::ublas::bounded_matrix<double, TDim, TDim> temp_type;
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> lambdamat =  prod(trans(emat), prod<temp_type>(Metric1Matrix,emat));
        const boost::numeric::ublas::bounded_matrix<double, TDim, TDim> mumat =  prod(trans(emat), prod<temp_type>(Metric2Matrix,emat));
        
        for (unsigned int i = 0; i < TDim; i++)
        {
            auxmat(i, i) = MathUtils<double>::Max(lambdamat(i, i), mumat(i, i));
        }
        
        boost::numeric::ublas::bounded_matrix<double, TDim, TDim> invemat;
        MathUtils<double>::InvertMatrix(emat, invemat, det);
        
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
