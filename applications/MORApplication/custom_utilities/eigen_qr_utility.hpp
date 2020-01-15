#ifndef EIGEN_QR_UTILITY_HPP
#define EIGEN_QR_UTILITY_HPP

// System includes
#include <vector>

// External includes
#include <Eigen/QR>
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"
// #include "custom_utilities/ublas_wrapper.h"
    
//     namespace Kratos {

// template <typename TScalar = double>
// class EigenSparseQRSolver
// {
// public:
//     using Scalar = TScalar;
//     using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor, int>;
//     using Vector = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

// private:
//     Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> m_solver;

// public:
//     static std::string Name()
//     {
//         return "eigen_sparse_qr";
//     }

/*void Eigentest()
{
    std::cout<<"Eigen Utility called"<<std::endl;
}*/

namespace Kratos
{
template <typename TSparseSpaceType>
class EigenQrUtility
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(EigenQrUtility);


    typedef typename TSparseSpaceType::DataType ScalarType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    // typename EigenSparseMatrixType = Eigen::SparseMatrix<ScalarType, Eigen::RowMajor, int>;

// EigenQrUtility(unsigned int system_size, const std::size_t n_sampling_points)
    EigenQrUtility()
    {
        // std::cout<<"Eigen constructor called\n";

        // std::cout<<"row size "<<system_size<<std::endl;
        // std::cout<<"column size "<<n_sampling_points<<std::endl;
        // Eigen::MatrixXd B(system_size, n_sampling_points);
        //  std::cout<<B.size()<<std::endl;

        //Eigen::Map<>
    }

    ~EigenQrUtility(){}

    void Initialize()
    {

    }

    void MatrixQ(SparseMatrixType& rA)
    {
        //sparse
        // std::vector<int> index_1 = std::vector<int>(rA.index1_data().begin(), rA.index1_data().end());
        // std::vector<int> index_2 = std::vector<int>(rA.index2_data().begin(), rA.index2_data().end());
        // Eigen::Map<const Eigen::SparseMatrix<ScalarType, Eigen::RowMajor, int>> A = Eigen::Map<const Eigen::SparseMatrix<ScalarType, Eigen::RowMajor, int>>
        //     ( rA.size1(), rA.size2(), rA.nnz(), index_1.data(), index_2.data(), rA.value_data().begin() );

        // BuiltinTimer aa;
        // Eigen::SparseMatrix<ScalarType> Q;
        // Q = Eigen::SparseQR<const Eigen::SparseMatrix<ScalarType, Eigen::RowMajor, int>, Eigen::COLAMDOrdering<int>>(A).matrixQ();
        // KRATOS_WATCH(aa.ElapsedSeconds())
        // KRATOS_WATCH(Q.outerSize())
        // KRATOS_WATCH(Q.innerSize())
        // SparseMatrixType tmp(Q.outerSize(), Q.innerSize(), Q.nonZeros());

        //dense
        Eigen::Map<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> A(rA.data().begin(), rA.size1(), rA.size2());
        Eigen::HouseholderQR<Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>> qr(A);
        Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic> Q;
        Eigen::MatrixXcd thinQ;
        Eigen::MatrixXd I(Eigen::MatrixXd::Identity(rA.size1(),rA.size2()));
        BuiltinTimer aa;
        Q = qr.householderQ();
        KRATOS_WATCH(aa.ElapsedSeconds())

        BuiltinTimer bb;
        thinQ = qr.householderQ() * I;
        KRATOS_WATCH(bb.ElapsedSeconds())
        KRATOS_WATCH(thinQ(1,2))

        // rA.data().begin() = thinQ.data();
        A = thinQ;


    }

private:
    // Eigen::SparseQR<const Eigen::SparseMatrix<ScalarType, Eigen::RowMajor, int>, Eigen::COLAMDOrdering<int>> mQR;


protected:

};

} // end class 
 // end namespace
#endif