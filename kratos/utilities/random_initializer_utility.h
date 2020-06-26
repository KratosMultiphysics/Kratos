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

// External includes

// Project includes
#include "spaces/ublas_space.h"

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
    
/**
 * @class RandomInitializeUtility
 * @ingroup KratosCore
 * @brief Utility to initialize a random vector
 * @details Defines several utility functions related with the initialization of random matrixes. The class can be adapted for several types of floating numbers via template
 * @author Vicente Mataix Ferrandiz
 */
template<class TDataType>
class RandomInitializeUtility
{
public:

    ///@name Type Definitions
    ///@{

    typedef UblasSpace<TDataType, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    
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
     * @brief This method initializes a vector using a normal normal distribution
     * @param R The vector to fill with random values
     * @param MeanValue The mean value used in the normal distribution
     * @param VarianceValue The variance value used in the normal distribution
     */
    static inline void NormalDestributionRandom(
        VectorType& R,
        const TDataType& MeanValue,
        const TDataType& VarianceValue
        )
    {
        // We create a random vector as seed of the method
        unsigned int seed = 1; // Constant seed
        std::default_random_engine generator(seed);
        
        std::normal_distribution<TDataType> normal_distribution(MeanValue, VarianceValue);
        
        for (SizeType i = 0; i < R.size(); ++i)
            R[i] = normal_distribution(generator);
    }
    
    
    /**
     * @brief This method initializes a vector using a normal distribution. The mean and the variance is taken from the norm of the matrix
     * @param K The stiffness matrix
     * @param R The vector to initialize
     * @param Inverse If consider the inverse pf the matrix norm or not
     */
    static inline void RandomInitialize(
        const SparseMatrixType& K,
        VectorType& R,
        const bool Inverse = false 
        )
    {
        const TDataType threshold = std::numeric_limits<TDataType>::epsilon();
        const TDataType normK = SparseSpaceType::TwoNorm(K);
        const TDataType aux_value = (Inverse == false) ? normK : (normK > threshold) ? 1.0/normK : 1.0;
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
