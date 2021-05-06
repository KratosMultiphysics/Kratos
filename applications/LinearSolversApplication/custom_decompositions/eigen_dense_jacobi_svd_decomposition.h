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

#if !defined(KRATOS_EIGEN_DENSE_JACOBI_SVD_H_INCLUDED)
#define KRATOS_EIGEN_DENSE_JACOBI_SVD_H_INCLUDED

// External includes
#include <Eigen/SVD>

// Project includes
#include "includes/define.h"
#include "linear_solvers_define.h"
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
class EigenDenseJacobiSVD : public DenseSingularValueDecomposition<TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    /// Definition of the shared pointer of the class
    KRATOS_CLASS_POINTER_DEFINITION(EigenDenseJacobiSVD);

    typedef typename TDenseSpace::DataType DataType;
    typedef typename TDenseSpace::VectorType VectorType;
    typedef typename TDenseSpace::MatrixType MatrixType;

    using EigenVector = Kratos::EigenDynamicVector<DataType>;
    using EigenMatrix = Kratos::EigenDynamicMatrix<DataType>;
    using DecompositionOptions = Eigen::DecompositionOptions;

    ///@}
    ///@name Life Cycle
    ///@{

    EigenDenseJacobiSVD() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static std::string Name()
    {
        return "dense_jacobi_singular_value_decomposition";
    }

    void Compute(
        MatrixType& rInputMatrix,
        Parameters Settings) override
    {
        // Check user defined options
        Parameters default_settings = Parameters(R"(
        {
            "compute_u" : true,
            "compute_v" : true,
            "compute_thin_u" : true,
            "compute_thin_v" : true
        })");
        Settings.ValidateAndAssignDefaults(default_settings);

        // Set the Eigen decomposition options
        int decomposition_options = 0;
        if (Settings["compute_u"].GetBool()) {
            decomposition_options |= Settings["compute_thin_u"].GetBool() ? Eigen::ComputeThinU : Eigen::ComputeFullU;
        }
        if (Settings["compute_v"].GetBool()) {
            decomposition_options |= Settings["compute_thin_v"].GetBool() ? Eigen::ComputeThinV : Eigen::ComputeFullV;
        }

        // Compute the Bidiagonal Divide and Conquer Singular Value Decomposition (BDCSVD)
        Eigen::Map<EigenMatrix> eigen_input_matrix_map(rInputMatrix.data().begin(), rInputMatrix.size1(), rInputMatrix.size2());
        mJacobiSVD.compute(eigen_input_matrix_map, decomposition_options);
    }

    void Compute(
        MatrixType& rInputMatrix,
        VectorType& rVectorS,
        MatrixType& rMatrixU,
        MatrixType& rMatrixV,
        Parameters Settings) override
    {
        Compute(rInputMatrix, Settings);
        SingularValues(rVectorS);
        MatrixU(rMatrixU);
        MatrixV(rMatrixV);
    }

    void MatrixU(MatrixType& rMatrixU) override
    {
        KRATOS_ERROR_IF_NOT(mJacobiSVD.computeU()) << "Matrix U has not been computed. Switch \'compute_u\' to \'true\' in the \'Compute\' input settings." << std::endl;

        const auto& r_matrix_U = mJacobiSVD.matrixU();
        const std::size_t m = r_matrix_U.rows();
        const std::size_t n = r_matrix_U.cols();
        if (rMatrixU.size1() != m || rMatrixU.size2() != n) {
            rMatrixU.resize(m,n);
        }

        Eigen::Map<EigenMatrix> matrix_u_map(rMatrixU.data().begin(), rMatrixU.size1(), rMatrixU.size2());
        matrix_u_map = r_matrix_U;
    }

    void MatrixV(MatrixType& rMatrixV) override
    {
        KRATOS_ERROR_IF_NOT(mJacobiSVD.computeV()) << "Matrix V has not been computed. Switch \'compute_v\' to \'true\' in the \'Compute\' input settings." << std::endl;

        const auto& r_matrix_V = mJacobiSVD.matrixV();
        const std::size_t m = r_matrix_V.rows();
        const std::size_t n = r_matrix_V.cols();
        if (rMatrixV.size1() != m || rMatrixV.size2() != n) {
            rMatrixV.resize(m,n);
        }

        Eigen::Map<EigenMatrix> matrix_v_map(rMatrixV.data().begin(), rMatrixV.size1(), rMatrixV.size2());
        matrix_v_map = r_matrix_V;
    }

    void SingularValues(VectorType& rVectorS) override
    {
        const auto& r_vector_S = mJacobiSVD.singularValues();
        const std::size_t m = r_vector_S.rows();
        if (rVectorS.size() != m) {
            rVectorS.resize(m);
        }

        Eigen::Map<EigenVector> vector_s_map(rVectorS.data().begin(), rVectorS.size());
        vector_s_map = r_vector_S;
    }

    std::size_t NonZeroSingularValues() override
    {
        return mJacobiSVD.nonzeroSingularValues();
    }

    std::size_t Rank() override
    {
        return mJacobiSVD.rank();
    }

    void PrintInfo(std::ostream &rOStream) const override
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

    Eigen::JacobiSVD<EigenMatrix> mJacobiSVD;

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

#endif // defined(KRATOS_EIGEN_DENSE_JACOBI_SVD_H_INCLUDED)
