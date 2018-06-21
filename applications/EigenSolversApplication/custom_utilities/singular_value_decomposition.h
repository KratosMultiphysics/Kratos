/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Long Chen
//
*/

#ifndef SINGULAR_VALUE_DECOMPOSITION_H
#define SINGULAR_VALUE_DECOMPOSITION_H

// External includes
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#if defined EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#endif
#include <Eigen/Sparse>

// Project includes
#include "spaces/ublas_space.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
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

    /// Short class definition.
    /** Detail class definition.
    */

    

    class SingularValueDecomposition
    {
    public:
        ///@name Type Definitions
        ///@{

        typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> dMatrix;
        typedef UblasSpace<double, Matrix, Vector> TDenseSpaceType;
        typedef TDenseSpaceType::MatrixType DenseMatrixType;
        typedef UblasSpace<double, CompressedMatrix, Vector> TSparseSpaceType;
        typedef TSparseSpaceType::VectorType VectorType;

        //typedef UblasSpace<double, CompressedMatrix, Vector> TSparseSpaceType;
        //typedef TSparseSpaceType::VectorType VectorType;

        //typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

        /// Pointer definition of MapperVertexMorphing
        KRATOS_CLASS_POINTER_DEFINITION(SingularValueDecomposition);

        /// Default constructor.
        SingularValueDecomposition()
        {
        }

        /// Destructor.
        virtual ~SingularValueDecomposition()
        {
        }

        void Solve(
            DenseMatrixType& rM,
            //DenseMatrixType& rM2,
            DenseMatrixType& rU,
            DenseMatrixType& rV,
            DenseMatrixType& rS)
        {
            std::cout << " SingularValueDecomposition::Solve starts!! "<< std::endl;
            std::cout << " SingularValueDecomposition::Solve starts!! "<< std::endl;

            using scalar_t = double;
            using vector_t = Eigen::VectorXd;
            using matrix_t = Eigen::MatrixXd;

            //rM = ZeroMatrix(2,4);
            //rM(0, 0) =  3;
            //rM(0, 1) =  4;
            //rM(0, 2) =  1;
            //rM(0, 3) =  9;

            //rM(1, 0) = -5;
            //rM(1, 1) =  -7;
            //rM(1, 2) =  -6;
            //rM(1, 3) =  4;

            VectorType M1 = row(rM, 0);
            VectorType M2 = row(rM, 1);

            std::cout << " M1 Norm:: "<< norm_2(M1) <<std::endl;

            std::cout << " M2 Norm:: "<< norm_2(M2) <<std::endl;
            
            double norm_M1 = norm_2(M1);
            double norm_M2 = norm_2(M2);
            if (norm_M1 != 0.0)
                M1 = M1 / norm_M1; // norm_2(M1);
            if (norm_M2 != 0.0)
                M2 = M2 / norm_M2; // norm_2(M2);
            
            row(rM, 0) = M1;
            row(rM, 1) = M2;


            
            Eigen::Map<matrix_t> svdMatrix(rM.data().begin(), rM.size1(), rM.size2());
            std::cout << " SingularValueDecomposition::1 "<< std::endl;
            Eigen::BDCSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > svd(svdMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            //Eigen::BDCSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > svd(svdMatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            std::cout << " SingularValueDecomposition::2 "<< std::endl;


            dMatrix V = svd.matrixV();
            rV.resize(V.rows(),V.cols());
            Eigen::Map<matrix_t> tmpV(rV.data().begin(), rV.size1(), rV.size2());
            tmpV = V;
            std::cout << "rV rows: " << std::endl;
            std::cout << V.rows() << std::endl;
            std::cout << "rV columns: " << std::endl;
            std::cout << V.cols() << std::endl;
            //std::cout << rV << std::endl;

            
            dMatrix U = svd.matrixU();
            rU.resize(U.rows(),U.cols());
            Eigen::Map<matrix_t> tmpU(rU.data().begin(), rU.size1(), rU.size2());
            tmpU = U;
            std::cout << "rU: " << std::endl;
            std::cout << rU << std::endl;    


            
            dMatrix S = svd.singularValues();
            rS.resize(S.rows(), S.cols());
            Eigen::Map<matrix_t> tmpS(rS.data().begin(), rS.size1(), rS.size2());
            //Eigen::Map<matrix_t> tmpU(rU.data().begin(), rU.size1(), rU.size2());
            tmpS = S;
            std::cout << "rS: " << std::endl;
            std::cout << rS << std::endl;   

            

        }
    };
} // namespace Kratos

#endif // defined(SINGULAR_VALUE_DECOMPOSITION_H)