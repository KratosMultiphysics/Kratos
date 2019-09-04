/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Authors: Matthias Ebert
//           based on code of
//             Thomas Oberbichler
//             Armin Geiser
*/

#if !defined(KRATOS_GEN_EIGENSYSTEM_SOLVER_H_INCLUDED)
#define KRATOS_GEN_EIGENSYSTEM_SOLVER_H_INCLUDED

// External includes
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/openmp_utils.h"
#include "custom_utilities/ublas_wrapper.h"

namespace Kratos
{

template<
    //class TSolverType,
    class TSparseSpaceType, // = typename TSolverType::TGlobalSpace,
    class TDenseSpaceType, // = typename TSolverType::TLocalSpace,
    class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
    class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType>>
class GenEigensystemSolver
    : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
    Parameters mParam;

  public:
    KRATOS_CLASS_POINTER_DEFINITION(GenEigensystemSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    GenEigensystemSolver(
        Parameters param
    ) : mParam(param)
    {
        Parameters default_params(R"(
        {
            "solver_type": "eigen_eigensystem",
            "number_of_eigenvalues": 1,
            "normalize_eigenvectors": false,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "echo_level": 1
        })");

        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance(mParam["tolerance"].GetDouble());
        BaseType::SetMaxIterationsNumber(mParam["max_iteration"].GetInt());
    }

    ~GenEigensystemSolver() override {}

    /**
     * Solve the generalized eigenvalue problem.
     *
     * K is a (non-)symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
        SparseMatrixType& rK,
        SparseMatrixType& rM,
        VectorType& rEigenvalues,
        //VectorType& rEigenvalues_im,
        //DenseMatrixType& rEigenvalues,
        DenseMatrixType& rEigenvectors) override
    {
        using scalar_t = double;
        //using vector_t = Eigen::VectorXcd;  // complex
        //using vector_t = Eigen::VectorXd;
        using matrix_t = Eigen::MatrixXd;

        // --- get settings

        const int nroot = mParam["number_of_eigenvalues"].GetInt();
        //const int max_iteration = BaseType::GetMaxIterationsNumber();
        //const double tolerance = BaseType::GetTolerance();
        const int echo_level = mParam["echo_level"].GetInt();

KRATOS_INFO("---this is the new one----");
        // --- wrap ublas matrices

        UblasWrapper<scalar_t> a_wrapper(rK);
        UblasWrapper<scalar_t> b_wrapper(rM);

        const auto& a = a_wrapper.matrix();
        const auto& b = b_wrapper.matrix();


        // --- timer

        double start_time = OpenMPUtils::GetCurrentTime();

        KRATOS_INFO_IF("GenEigensystemSolver:", echo_level > 0) << "Start"  << std::endl;

        // --- calculation
        //int nn = a.rows();  //TODO: rename nn      and define nroot*2 for n_complex_root or so

        Eigen::GeneralizedEigenSolver<matrix_t> eig;
        eig.compute(a, b);

            if(eig.info() != Eigen::Success) {
                KRATOS_WARNING("GenEigensystemSolver:") << "Eigen solution was not successful!" << std::endl;
            }

        // debug print, to delete
        KRATOS_WATCH(nroot);
        KRATOS_WATCH(eig.eigenvalues());
        KRATOS_WATCH(eig.eigenvectors());


        if (static_cast<int>(rEigenvalues.size()) != nroot*2) {
            rEigenvalues.resize(nroot*2);
        }
        if (static_cast<int>(rEigenvectors.size1()) != nroot || static_cast<int>(rEigenvectors.size2()) != nroot*2) {
            rEigenvectors.resize(nroot, nroot*2);
        }

 
        // complex format in double matrix
        // // if (static_cast<int>(rEigenvalues.size1()) != nroot || static_cast<int>(rEigenvalues.size2()) != 2) {
        // //     rEigenvalues.resize(nroot,2);
        // // }
        // // if (static_cast<int>(rEigenvectors.size1()) != nroot*2 || static_cast<int>(rEigenvectors.size2()) != nn) {
        // //     rEigenvectors.resize(nroot*2, nn);
        // // }



        // // // // // KRATOS_WATCH(rEigenvalues);

        // // // // // //rEigenvalues(1) = 2.0;


        // // // // // KRATOS_WATCH(rEigenvalues);

        // // // // // //Eigen::Map<matrix_t> rEigenvalues (eig.eigenvalues(),nroot,2, Eigen::Stride<0,0>);

        // // // // // KRATOS_WATCH(rEigenvalues);


        // // // // // Eigen::VectorXcd Vt = eig.eigenvalues();

        // // // // // KRATOS_WATCH(Vt);

        // // // // // KRATOS_WATCH(Vt(0));
        // // // // // KRATOS_WATCH(Vt(1));

        // // // // // KRATOS_WATCH(Vt(0,0));
        // // // // // KRATOS_WATCH(Vt(0,1));

        // // // // // std::complex<double> t_c = Vt(0);

        // // // // // KRATOS_WATCH(t_c);

        // // // // // KRATOS_WATCH(t_c.real());
        // // // // // KRATOS_WATCH(t_c.imag());


        // // // // // double t_a[2] = {0,0};
        // // // // // t_a[0] = t_c.real();
        // // // // // t_a[0] = t_c.imag();

        // // // // // KRATOS_WATCH(t_a);

        // // // // // rEigenvalues(0) = t_c.real();
        // // // // // rEigenvalues(1) = t_c.imag();
        // // // // // KRATOS_WATCH(rEigenvalues);

        // // // // // rEigenvalues(2) = Vt(0).real();
        // // // // // rEigenvalues(3) = Vt(0).imag();
        // // // // // KRATOS_WATCH(rEigenvalues);

        // // // // // //rEigenvalues_im(0) = Vt(0).imag();
        // // // // // //rEigenvalues_im(1) = Vt(1).imag();


        Eigen::VectorXcd eigvals_tmp = eig.eigenvalues();
        for(int i=0; i<nroot; i++){
            rEigenvalues(2*i) = eigvals_tmp(i).real();
            rEigenvalues(2*i+1) = eigvals_tmp(i).imag();
        }


        Eigen::MatrixXcd eigvecs_tmp = eig.eigenvectors();
        for(int i=0; i<nroot; i++){
            for(int j=0; j<nroot; j++){
                rEigenvectors(i,2*j) = eigvecs_tmp(i,j).real();
                rEigenvectors(i,2*j+1) = eigvecs_tmp(i,j).imag();
            }
        }
    



 











        //Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size(),2, Eigen::Innerstride<2>);
        //Eigen::Map<vector_t> eigvals (rEigenvalues.data().begin(), rEigenvalues.size());
        //Eigen::Map<matrix_t> eigvecs (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());
        //Eigen::Map<matrix_t> rEigenvectors (rEigenvectors.data().begin(), rEigenvectors.size1(), rEigenvectors.size2());
//eigvals = eig.eigenvalues();

//KRATOS_WATCH(eigvals);

// what happens heree? solver = Eigen::SparseLU
/*
        eigvals = eig.eigenvalues().head(nroot);

        for (int i = 0; i != nroot; ++i) {
            tmp = r.col(i);
            eigvecs.row(i) = solver.solve(tmp).normalized();
        }
*/

// TODO: could be reused
/*

        // --- normalization
        // Given generalized eigenvalue problem (A - eigenvalue * B) * eigenvector = 0,
        // eigenvector is normalized such that eigenvector^T * B * eigenvector = 1
        if(mParam["normalize_eigenvectors"].GetBool())
        {
            for (int i = 0; i != nroot; ++i)
            {
                const double tmp = eigvecs.row(i) * b * eigvecs.row(i).transpose();
                const double factor = 1.0 / std::sqrt(tmp);
                eigvecs.row(i) *=  factor;
                KRATOS_INFO_IF("GenEigensystemSolver:", echo_level > 0) << "Eigenvector " << i+1 << " is normalized - used factor: " << factor << std::endl;
            }
        }
*/

// TODO: print eigenvalues
        // --- output
        if (echo_level > 0) {
            double end_time = OpenMPUtils::GetCurrentTime();
            double duration = end_time - start_time;

            Eigen::IOFormat fmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");

            KRATOS_INFO("GenEigensystemSolver:") << "Completed in " << duration << " seconds" << std::endl;
          //            << "                   Eigenvalues = " << eigvals.transpose().format(fmt) << std::endl;
        }
    }

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "GenEigensystemSolver";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

}; // class GenEigensystemSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(
    std::istream& rIStream,
    GenEigensystemSolver<TSparseSpaceType,
    TDenseSpaceType,
    TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(
    std::ostream& rOStream,
    const GenEigensystemSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_GEN_EIGENSYSTEM_SOLVER_H_INCLUDED)
