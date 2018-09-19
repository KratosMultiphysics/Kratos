// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_METRICS_MATH_UTILS)
#define KRATOS_METRICS_MATH_UTILS

// Project includes
#include "utilities/math_utils.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
    typedef std::size_t SizeType;

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
 * @class MetricsMathUtils
 * @ingroup MeshingApplication
 * @brief This class is used to compute some mathematical operations needed for the metrics computing
 * @author Vicente Mataix Ferrandiz
 */
template<SizeType TDim>
class MetricsMathUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// The type of array considered for the tensor
    typedef typename std::conditional<TDim == 2, array_1d<double, 3>, array_1d<double, 6>>::type TensorArrayType;

    /// The definition of the matrix type
    typedef BoundedMatrix<double, TDim, TDim> MatrixType;

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
     * @brief It computes the intersection between two metrics
     * @param Metric1 The first metric
     * @param Metric2 The second metric
     * @return The intersected metric
     */
    static inline TensorArrayType IntersectMetrics(
        const TensorArrayType& Metric1,
        const TensorArrayType& Metric2
        )
    {
        const MatrixType& metric1_matrix = MathUtils<double>::VectorToSymmetricTensor<TensorArrayType, MatrixType>(Metric1);
        const MatrixType& metric2_matrix = MathUtils<double>::VectorToSymmetricTensor<TensorArrayType, MatrixType>(Metric2);

        MatrixType auxmat, emat;

        double auxdet;
        const MatrixType& inv_metric1_matrix = MathUtils<double>::InvertMatrix<TDim>(metric1_matrix, auxdet);
        const MatrixType n_matrix = prod(inv_metric1_matrix, metric2_matrix);

        MathUtils<double>::EigenSystem<TDim>(n_matrix, emat, auxmat, 1e-18, 20);

        typedef MatrixType temp_type;
        const MatrixType lambdamat =  prod(trans(emat), prod<temp_type>(metric1_matrix, emat));
        const MatrixType mumat =  prod(trans(emat), prod<temp_type>(metric2_matrix, emat));

        for (std::size_t i = 0; i < TDim; ++i)
            auxmat(i, i) = MathUtils<double>::Max(lambdamat(i, i), mumat(i, i));

        const MatrixType& invemat = MathUtils<double>::InvertMatrix<TDim>(emat, auxdet);

        MatrixType IntersectionMatrix =  prod(trans(invemat), prod<temp_type>(auxmat, invemat));

        const TensorArrayType& Intersection = MathUtils<double>::StressTensorToVector<MatrixType, TensorArrayType>(IntersectionMatrix);

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
