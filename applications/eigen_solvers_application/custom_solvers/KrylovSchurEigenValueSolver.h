/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Armin Geiser
*/

// TODO check naming: is this really a Krylov-Schur algorithm?
#if !defined(KRATOS_KRYLOV_SCHUR_EIGEN_VALUE_SOLVER_H_INCLUDED)
#define KRATOS_KRYLOV_SCHUR_EIGEN_VALUE_SOLVER_H_INCLUDED

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "linear_solvers/iterative_solver.h"

#include <Eigen/Core>
#include <Eigen/Eigenvalues>


namespace Kratos
{

///@name Kratos Classes
///@{

/// Adapter TODO
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class KrylovSchurEigenValueSolver: public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{

  public:
    KRATOS_CLASS_POINTER_DEFINITION(KrylovSchurEigenValueSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    KrylovSchurEigenValueSolver(Parameters::Pointer pParam,
        typename TLinearSolverType::Pointer pLinearSolver) : mpParam(pParam), mpLinearSolver(pLinearSolver)

    {

        Parameters default_params(R"(
        {
            "solver_type": "KrylovSchurEigenValueSolver",
            "number_of_eigenvalues": 1,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "verbosity": 1,
            "linear_solver_settings": {}
        })");

        // don't validate linear_solver_settings here
        mpParam->ValidateAndAssignDefaults(default_params);

        Parameters& settings = *mpParam;
        BaseType::SetTolerance( settings["tolerance"].GetDouble() );
        BaseType::SetMaxIterationsNumber( settings["max_iteration"].GetInt() );

    }

    ~KrylovSchurEigenValueSolver() override {}


    /**
     * Solve the generalized eigenvalue problem using an eigen subspace iteration method
     * The implementation follows the code from
     * K. J. Bathe, Finite Element Procedures second Edition, ISBN-13: 978-0979004957
     * page 954 and following
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& rK,
            SparseMatrixType& rM,
            VectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors) override
    {
        std::cout << "Start solving for eigen values."  << std::endl;

        // settings
        Parameters& settings = *mpParam;

        const int verbosity = settings["verbosity"].GetInt();
        int NITEM = settings["max_iteration"].GetInt();
        int _NROOT = settings["number_of_eigenvalues"].GetInt(); // number of eigenvalues requested
        double RTOL = settings["tolerance"].GetDouble();

        int nn;  // size of problem (order of stiffness and mass matrix)
        int nc;  // number of iteration vectors used (automatically computed)
        int nc1; // = nc-1
        int ij=0;
        int nite;
        int is;
        double art;
        double brt;
        double rt;
        double eigvt;
        double dif;
        bool eigen_solver_successful = true;

        nn = rK.size1();
        nc = std::min(2*_NROOT,_NROOT+8);
        if (nc > nn ) nc = nn;  //nc cannot be larger than number of mass degree of freedom
        nc1 = nc - 1 ;


        // ar: working matrix storing projection of rK
        DenseMatrixType ar = ZeroMatrix(nc,nc);

        // br: working matrix storing projection of rM
        DenseMatrixType br = ZeroMatrix(nc,nc);

        // declaration and initialization of vectors
        VectorType w = ZeroVector(nn);                  // working vector

        VectorType d = ZeroVector(nc);                  // working vector

        VectorType tt = ZeroVector(nn);                 // working vector

        VectorType temp_tt(nn);                         // working vector (contains the exited DOF on exit)

        VectorType rtolv = ZeroVector(nc);              // working vector

        VectorType temp(nn);           // working vector, not needed at Bathe

        // vec: eigenvectors for reduced problem
        DenseMatrixType vec = ZeroMatrix(nc,nc);        // working matrix

        // TODO here maybe teh rEigenvAlue... can be used
        // create working arrays for eigenvectors and eigenvalues
        DenseMatrixType r = ZeroMatrix(nn,nc);          // storage for eigenvectors

        VectorType eigv = ZeroVector(nc);               // storage for eigenvalues

        //------------------------------------------------------------------------------------
        // establish starting iteration vectors (column-wise in matrix r):
        //
        // starting iteration vectors should be created such that they excite DOFs with large
        // mass and small stiffness.
        //
        // Alternative A: Approach according to Bathe:
        // The first column of r is filled with diagonal entries of mass matrix, the others
        // (except the last) with unit vectors containing the 1 at the positions of DOFs with
        // highest ratio mass/stiffness.
        //
        // Alternative B: Approach according to oofem:
        // The first column of r is filled with the result of mass matrix times a vector with
        // all entries equal 1. The other vectors are filled with the result of mass matrix
        // times unity vectors containing the 1 at the positions of DOFs with highest ratio
        // mass/stiffness.
        //------------------------------------------------------------------------------------

        //Alternative A: Approach according to Bathe:
        //-------------------------------------------

        // case: mass matrix is a diagonal matrix (lumped)
        if(int(rM.nnz()) == nn)
        {
            int j=0;
            for (int i = 1; i<= nn ; i++)
            {
                r(i-1,0)=rM(i-1,i-1);  // first vector is diagonal of mass matrix
                if(rM(i-1,i-1)>0.0)
                    j++;
                w(i-1) = rM(i-1,i-1)/rK(i-1,i-1);
            }
            if(nc>j)
            {
                KRATOS_ERROR <<"Number of needed iteration vectors (nc) is larger than number of mass DOF." << std::endl;
            }
        }
        else // case: mass matrix is not a diagonal matrix
        {
            for (int i = 1; i<= nn ; i++)
            {
                r(i-1,0)=rM(i-1,i-1);  // first vector is diagonal of mass matrix
                w(i-1) = rM(i-1,i-1)/rK(i-1,i-1);
            }
        }

        // For the computation of large systems, suitable distances between the unity entries
        // are needed in the starting iteration vectors in r. This is taken into account via
        // variable l.
        int nd = nn/nc;
        int l = nn - nd;

        for (int j = 2; j<= nc ; j++)
        {
            rt = 0.0;
            for (int i = 1; i<= l ; i++)
            {
                if (w(i-1) >= rt)
                {
                    rt = w(i-1);
                    ij = i;
                }
            }
            for (int i = l; i<= nn ; i++)
            {
                if (w(i-1) > rt)
                {
                    rt = w(i-1);
                    ij = i;
                }
            }
            tt(j-1) = double(ij); //gets the equation-nrs of excited DOFs
            w(ij-1) = 0.0;

            l -= nd;

            r(ij-1,j-1) = 1.0;
        }
        if (verbosity>1) {
            std::cout << "Degrees of freedom excited by unit starting iteration vectors:" << std::endl;
            for (int j = 2; j <= nc; j++)
                std::cout << tt(j-1) << std::endl << std::endl;
        }
        // r has now the initial values for the eigenvectors

        // TODO add random numbers to last column of r(:,nc) -> Bathe SSP00131

        if (verbosity>1) std::cout << "Initialize linear solver." << std::endl;
        mpLinearSolver->InitializeSolutionStep(rK, temp_tt, tt);

        //------------------------------------------------------------------------------------
        // start of iteration loop
        //------------------------------------------------------------------------------------
        nite = 0;
        do // label 100
        {
            nite ++;

            if (verbosity>1) std::cout << "Subspace Iteration Eigenvalue Solver: Iteration no. " << nite <<std::endl;

            // compute projection ar of matrix casted_A
            for (int j = 1; j<= nc; j++)
            {
                for (int k = 1; k<= nn; k++) // tt gets an iteration eigenvector
                    tt(k-1) = r(k-1,j-1);

                // K*temp_tt = tt
                if (verbosity>1) std::cout << "Backsubstitute using initialized linear solver." << std::endl;
                mpLinearSolver->PerformSolutionStep(rK, temp_tt, tt);
                for(int k=1; k<=nn; k++) tt(k-1) = temp_tt(k-1);

                for (int i = j; i<= nc; i++) //SSP00161
                {
                    art = 0.;
                    for (int k = 1; k<= nn; k++)
                        art += r(k-1,i-1)*tt(k-1);
                    ar(j-1,i-1) = art;
                    ar(i-1,j-1) = art;
                }
                for (int k = 1; k<= nn; k++)
                    r(k-1,j-1) = tt(k-1);   // (r = xbar)
            }

            // compute projection br of matrix casted_B
            for (int j = 1; j<= nc; j++) //SSP00170
            {
                for (int k = 1; k<= nn ; k++)  // temp gets an iteration eigenvector
                    temp(k-1)=r(k-1,j-1);

                // tt = rM*temp
                TSparseSpaceType::Mult(rM, temp, tt);

                for (int i = j; i<= nc; i++)
                {
                    brt = 0.;
                    for (int k = 1; k<= nn; k++)
                        brt += r(k-1,i-1)*tt(k-1);
                    br(j-1,i-1) = brt;
                    br(i-1,j-1) = brt;
                }                         // label 180

                for (int k = 1; k<= nn; k++)
                    r(k-1,j-1)=tt(k-1);

            }                           // label 160

            //solve sub-problem
            Eigen::Map<Eigen::MatrixXd> AR(&ar(0,0),nc,nc);
            Eigen::Map<Eigen::MatrixXd> BR(&br(0,0),nc,nc);

            Eigen::VectorXd evalues;
            Eigen::MatrixXd evecs;

            if (verbosity>1) std::cout << "Start eigen value solution of a nxn matrix with n=" << nc << std::endl;
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(AR, BR);

            if (ges.info() != Eigen::Success)
                eigen_solver_successful = false;

            evalues = ges.eigenvalues();
            evecs = ges.eigenvectors();

            eigv.resize(evalues.rows());
            vec.resize(evecs.rows(), evecs.cols());

            for (int i=0; i<evalues.rows(); i++){
                eigv(i) = evalues(i);
                for (int j=0; j<evecs.cols(); j++){
                    vec(i,j) = evecs(i,j);
                }
            }

            if(!eigen_solver_successful)
            {
                std::cout << "Eigen solution was not successful!" << std::endl;
                break;
            }

            // sorting eigenvalues according to their values
            do
            {
                is = 0; // label 350
                for (int i = 1; i<= nc1; i++)
                {
                    if ((eigv(i)) < (eigv(i-1)))
                    {
                        is++;
                        eigvt = eigv(i);
                        eigv(i) = eigv(i-1) ;
                        eigv(i-1) = eigvt ;
                        for (int k = 1; k<= nc; k++)
                        {
                            rt = vec(k-1,i) ;
                            vec(k-1,i) = vec(k-1,i-1) ;
                            vec(k-1,i-1)   = rt ;
                        }
                    }
                }                         // label 360
            } while (is != 0);

            // compute eigenvectors
            for (int i = 1; i<= nn; i++)  // label 375
            {
                for (int j = 1; j<= nc; j++)
                    tt(j-1) = r(i-1,j-1);
                for (int k = 1; k<= nc; k++)
                {
                    rt = 0.0;
                    for (int j = 1; j<= nc; j++)
                        rt += tt(j-1)*vec(j-1,k-1) ;
                    r(i-1,k-1) = rt ;
                }
            }                           // label 420   (r = z)

            // convergency check alternative A:
            for (int i = 1; i<= nc; i++)
            {
                dif = (eigv(i-1) - d(i-1));
                rtolv(i-1) = fabs(dif / eigv(i-1));
            }

            // convergency check alternative B: according to BATHE (more restrictive)
            // for (int i = 1; i<=nc; i++)
            // {
            //     double vdot = 0.0;
            //     for (int j=1; j<=nc; j++)
            //         vdot += vec(j-1,i-1);
            //     double eigv2 = eigv(i-1) * eigv(i-1);
            //     dif = vdot-eigv2;
            //     rtolv(i-1) = sqrt(std::max(dif, 1e-24*eigv2)/vdot);
            // }

            for (int i = 1; i<= _NROOT; i++)
            {
                if (rtolv(i-1) > RTOL) goto label400 ;
            }

            std::cout << "Convergence reached after " << nite << " iterations within a relative tolerance: " << RTOL << std::endl;

            break;

            label400:      // "not converged so far"
            if (nite >= NITEM)
            {
                std::cout << "Convergence not reached in " << NITEM << " iterations." << std::endl;
                break;
            }

            for (int i = 1; i<= nc ; i++)
                d(i-1) = eigv(i-1); // label 410 and 440

            //GramSchmidt Orthogonalization
            // GramSchmidt(r);
        } while (1);

        mpLinearSolver->FinalizeSolutionStep(rK, temp_tt, tt);

        if(eigen_solver_successful)
        {
            // compute eigenvectors
            for (int j = 1; j<= nc; j++)
            {
                for (int k = 1; k<= nn; k++)
                    tt(k) = r(k-1,j-1);

                // K*temp_tt = tt
                mpLinearSolver->PerformSolutionStep(rK, temp_tt, tt);
                for(int k=1; k<=nn; k++)
                    tt(k-1) = temp_tt(k-1);

                for (int k = 1; k<= nn; k++)
                    r(k-1,j-1) = tt(k-1);   // (r = xbar)
            }

            // copy results to function parameters
            rEigenvalues.resize(_NROOT);
            rEigenvectors.resize(_NROOT, nn);
            for (int i = 0; i< _NROOT; i++)
            {
                rEigenvalues(i) = eigv(i);
                for (int j=0; j< nn; j++)
                    rEigenvectors(i,j) = r(j,i);
            }

            KRATOS_WATCH(rEigenvalues);

        // TODO
        // 1. make vectors M-normalized
        // 2. orient vectors
        // 3. speed up
        // 4. sturm sequence check
        }
        else{
            KRATOS_ERROR << "Solution failed!" <<std::endl;
        }

        return;
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "KrylovSchurEigenValueSolver.";
    }

    /**
     * Print object's data.
     */
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///@}

  private:
    ///@name Member Variables
    ///@{

    Parameters::Pointer mpParam;

    typename TLinearSolverType::Pointer mpLinearSolver;

    ///@}

}; // class KrylovSchurEigenValueSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        KrylovSchurEigenValueSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const KrylovSchurEigenValueSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_KRYLOV_SCHUR_EIGEN_VALUE_SOLVER_H_INCLUDED)