//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_RANDOM_UTILITY_INITIALIZER_H_INCLUDED )
#define  KRATOS_RANDOM_UTILITY_INITIALIZER_H_INCLUDED

// System includes

#include <random>
#include <boost/range/algorithm.hpp>

// External includes

// Project includes
#include "spaces/ublas_space.h"
#include "includes/ublas_interface.h"

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
    
/// Utility to initialize a random vector
/**
 * Defines several utility functions
 */
template<class TDataType>
class RandomInitializeUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef UblasSpace<TDataType, CompressedMatrix, Vector> SparseSpaceType;
    
    typedef UblasSpace<TDataType, Matrix, Vector> LocalSpaceType;
    
    typedef typename SparseSpaceType::MatrixType SparseMatrixType;

    typedef typename SparseSpaceType::VectorType VectorType;

    typedef typename LocalSpaceType::MatrixType DenseMatrixType;

    typedef typename LocalSpaceType::VectorType DenseVectorType;
    
    typedef std::size_t SizeType;

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
     * This method initializes a vector using a normal normal distribution
     * @param R: The vector to fill with random values
     * @param MeanValue: The mean value used in the normal distribution
     * @param VarianceValue: The variance value used in the normal distribution
     */
    static inline void NormalDestributionRandom(
        DenseVectorType& R,
        const TDataType& MeanValue,
        const TDataType& VarianceValue
        )
    {
        // We create a random vector as seed of the method
        std::random_device this_random_device;
        std::mt19937 generator(this_random_device());
        
        std::normal_distribution<> normal_distribution(MeanValue, VarianceValue);
        
        for (SizeType i = 0; i < R.size(); ++i)
        {
            R[i] = normal_distribution(generator);
        }
    }
    
    
    /**
     * This method initializes a vector using a normal distribution. The mean and the variance is taken from the norm of the matrix
     * @param K: The stiffness matrix
     * @param R: The vector to initialize
     * @param Inverse: If consider the inverse pf the matrix norm or not
     */
    static inline void RandomInitialize(
        const SparseMatrixType& K,
        DenseVectorType& R,
        const bool Inverse = false 
        )
    {
        const TDataType normK = SparseSpaceType::TwoNorm(K);
        const TDataType aux_value = (Inverse == false) ? normK : 1.0/normK;
        NormalDestributionRandom(R, aux_value, 0.25 * aux_value);
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

    RandomInitializeUtility(void);

    RandomInitializeUtility(RandomInitializeUtility& rSource);

}; /* Class RandomInitializeUtility */

}  // namespace Kratos.

#endif // KRATOS_RANDOM_UTILITY_INITIALIZER_H_INCLUDED defined
