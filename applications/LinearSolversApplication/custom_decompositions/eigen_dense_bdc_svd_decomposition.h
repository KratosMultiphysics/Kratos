//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_EIGEN_DENSE_BDCSVD_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_BDCSVD_H_INCLUDED

// External includes
#include <Eigen/SVD>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "utilities/dense_svd_decomposition.h"

namespace Kratos {

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

template<class TDenseSpace>
class EigenDenseBDCSVD : public DenseSingularValueDecomposition<TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(EigenDenseBDCSVD);

    typedef typename TDenseSpace::DataType DataType;
    typedef typename TDenseSpace::VectorType VectorType;
    typedef typename TDenseSpace::MatrixType MatrixType;

    using EigenVector = Eigen::Matrix<DataType, Eigen::Dynamic, 1>;
    using EigenMatrix = Eigen::Matrix<DataType, Eigen::Dynamic, Eigen::Dynamic>;
    using DecompositionOptions = Eigen::DecompositionOptions;

    ///@}
    ///@name Life Cycle
    ///@{

    EigenDenseBDCSVD() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static std::string Name()
    {
        return "dense_bidiagonal_divide_and_conquer_singular_value_decomposition";
    }

    void Compute(
        MatrixType& rInputMatrix,
        Parameters Settings)
    {
        // Check user defined options
        Parameters default_settings = Parameters(R"(
        {
            "compute_U" : true,
            "compute_V" : true,
            "compute_thin_U" : true,
            "compute_thin_V" : true
        })");
        Settings.ValidateAndAssignDefaults(default_settings);

        // Set the Eigen decomposition options
        int decomposition_options = 0;
        if (Settings["compute_U"].GetBool()) {
            decomposition_options |= Settings["compute_thin_U"].GetBool() ? Eigen::ComputeThinU : Eigen::ComputeFullU;
        }
        if (Settings["compute_V"].GetBool()) {
            decomposition_options |= Settings["compute_thin_V"].GetBool() ? Eigen::ComputeThinV : Eigen::ComputeFullV;
        }

        // Compute the Bidiagonal Divide and Conquer Singular Value Decomposition (BDCSVD)
        Eigen::Map<EigenMatrix> eigen_input_matrix_map(rInputMatrix.data().begin(), rInputMatrix.size1(), rInputMatrix.size2());
        mBDCSVD.compute(eigen_input_matrix_map, decomposition_options);
    }

    void Compute(
        MatrixType& rInputMatrix,
        VectorType& rVectorS,
        MatrixType& rMatrixU,
        MatrixType& rMatrixV,
        Parameters Settings)
    {
        Compute(rInputMatrix, Settings);
        SingularValues(rVectorS);
        MatrixU(rMatrixU);
        MatrixV(rMatrixV);
    }

    void MatrixU(MatrixType& rMatrixU)
    {
        KRATOS_ERROR_IF_NOT(mBDCSVD.computeU()) << "Matrix U has not been computed. Switch \'compute_U\' to \'true\' in the \'Compute\' input settings." << std::endl;

        const auto& r_matrix_U = mBDCSVD.matrixU();
        const std::size_t m = r_matrix_U.rows();
        const std::size_t n = r_matrix_U.cols();
        if (rMatrixU.size1() != m || rMatrixU.size2() != n) {
            rMatrixU.resize(m,n);
        }

        Eigen::Map<EigenMatrix> matrix_u_map(rMatrixU.data().begin(), rMatrixU.size1(), rMatrixU.size2());
        matrix_u_map = r_matrix_U;
    }

    void MatrixV(MatrixType& rMatrixV)
    {
        KRATOS_ERROR_IF_NOT(mBDCSVD.computeV()) << "Matrix V has not been computed. Switch \'compute_V\' to \'true\' in the \'Compute\' input settings." << std::endl;

        const auto& r_matrix_V = mBDCSVD.matrixV();
        const std::size_t m = r_matrix_V.rows();
        const std::size_t n = r_matrix_V.cols();
        if (rMatrixV.size1() != m || rMatrixV.size2() != n) {
            rMatrixV.resize(m,n);
        }

        Eigen::Map<EigenMatrix> matrix_v_map(rMatrixV.data().begin(), rMatrixV.size1(), rMatrixV.size2());
        matrix_v_map = r_matrix_V;
    }

    void SingularValues(VectorType& rVectorS)
    {
        const auto& r_vector_S = mBDCSVD.singularValues();
        const std::size_t m = r_vector_S.rows();
        if (rVectorS.size() != m) {
            rVectorS.resize(m);
        }

        Eigen::Map<EigenVector> vector_s_map(rVectorS.data().begin(), rVectorS.size());
        vector_s_map = r_vector_S;
    }

    std::size_t NonZeroSingularValues()
    {
        return mBDCSVD.nonzeroSingularValues();
    }

    std::size_t Rank()
    {
        return mBDCSVD.rank();
    }

    void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "EigenDecomposition <" << Name() << "> finished.";
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


    ///@}
private:
    ///@name Private static Member Variables
    ///@{


    ///@}
    ///@name Private member Variables
    ///@{

    Eigen::BDCSVD<EigenMatrix> mBDCSVD;

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


    ///@}
};

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
} // namespace Kratos

#endif // defined(KRATOS_EIGEN_DENSE_BDCSVD_H_INCLUDED)
