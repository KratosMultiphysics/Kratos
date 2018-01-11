//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Armin Geiser
//
//

#if !defined(KRATOS_SUBSPACE_ITERATION_EIGEN_VALUE_SOLVER_H_INCLUDED)
#define KRATOS_SUBSPACE_ITERATION_EIGEN_VALUE_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/iterative_solver.h"
#include "utilities/openmp_utils.h" // for timer

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

/// This class uses the subspace iteration method to obtain the n lowest eigenvalues of a system
/// Gives slightly different results the e.g. scipys geigs function
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SubspaceIterationEigenvalueSolver: public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{

  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SubspaceIterationEigenvalueSolver
    KRATOS_CLASS_POINTER_DEFINITION(SubspaceIterationEigenvalueSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    SubspaceIterationEigenvalueSolver(Parameters param,
        typename TLinearSolverType::Pointer pLinearSolver,
        typename TLinearSolverType::Pointer pEigenValueSolver) : mParam(param),
                                                                 mpLinearSolver(pLinearSolver),
                                                                 mpEigenValueSolver(pEigenValueSolver)

    {

        Parameters default_params(R"(
        {
            "solver_type": "SubspaceIterationEigenvalueSolver",
            "number_of_eigenvalues": 1,
            "max_iteration": 1000,
            "tolerance": 1e-6,
            "mass_normalization": true,
            "orient_eigen_vectors": true,
            "echo_level": 1,
            "linear_solver_settings": {},
            "eigen_sub_solver_settings" :{}
        })");

        // don't validate linear_solver_settings and eigen_sub_solver_settings here
        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance( mParam["tolerance"].GetDouble() );
        BaseType::SetMaxIterationsNumber( mParam["max_iteration"].GetInt() );

        KRATOS_ERROR_IF( mParam["linear_solver_settings"]["solver_type"].GetString() == "eigen_pardiso_llt") <<
            "eigen_pardiso_llt cannot handle negative entries on the main diagonal" << std::endl;

        std::cout << std::endl;
        std::cout << "WARNING: SubspaceIterationEigenvalueSolver showed slightly different results than e.g. Arnoldi method." << std::endl << std::endl;

        std::cout << "WARNING: Make sure the linear solver can handle negative entries on the main diagonal!" << std::endl << std::endl;

    }

    /// Destructor.
    ~SubspaceIterationEigenvalueSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Solve the generalized eigenvalue problem using an eigen subspace iteration method
     * The implementation follows the code from
     * K. J. Bathe, Finite Element Procedures second Edition, ISBN-13: 978-0979004957
     * page 954 and following
     * The naming of the variables is chose according to the reference.
     *
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& rK,
            SparseMatrixType& rM,
            VectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors) override
    {
        double start_time = OpenMPUtils::GetCurrentTime();

        const int echo_level = mParam["echo_level"].GetInt();
        const int nroot = mParam["number_of_eigenvalues"].GetInt(); // number of eigenvalues requested
        const int nitem = BaseType::GetMaxIterationsNumber();
        const double rtol = BaseType::GetTolerance();

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

        if (echo_level > 0) std::cout << "Start subspace iteration to solve for eigen values."  << std::endl;

        nn = rK.size1();
        nc = std::min(2*nroot,nroot+8);
        if (nc > nn)
            nc = nn;  //nc cannot be larger than number of mass degree of freedom
        nc1 = nc - 1 ;

        // ar: working matrix storing projection of rK
        //DenseMatrixType ar = ZeroMatrix(nc,nc);
        SparseMatrixType ar(nc,nc);

        // br: working matrix storing projection of rM
        //DenseMatrixType br = ZeroMatrix(nc,nc);
        SparseMatrixType br(nc,nc);

        // declaration and initialization of vectors
        VectorType w = ZeroVector(nn);                  // working vector

        VectorType tt = ZeroVector(nn);                 // working vector

        VectorType temp_tt(nn);                         // working vector (contains the exited DOF on exit)

        VectorType rtolv = ZeroVector(nc);              // working vector

        VectorType temp(nn);           // working vector, not needed at Bathe

        // vec: eigenvectors for reduced problem
        DenseMatrixType vec = ZeroMatrix(nc,nc);        // working matrix

        // TODO here maybe the rEigenvalue... can be used
        // create working arrays for eigenvectors and eigenvalues
        DenseMatrixType r = ZeroMatrix(nn,nc);          // storage for eigenvectors

        VectorType eigv = ZeroVector(nc);               // storage for eigenvalues

        VectorType prev_eigv = ZeroVector(nc);                  // working vector

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
            for (int i = 0; i< nn ; i++)
            {
                r(i,0)=rM(i,i);  // first vector is diagonal of mass matrix
                if(rM(i,i)>0.0)
                    j++;
                w(i) = rM(i,i)/rK(i,i);
            }
            if(nc>j)
            {
                KRATOS_ERROR <<"Number of needed iteration vectors (nc) is larger than number of mass DOF." << std::endl;
            }
        }
        else // case: mass matrix is not a diagonal matrix
        {
            for (int i = 0; i< nn ; i++)
            {
                r(i,0)=rM(i,i);  // first vector is diagonal of mass matrix
                w(i) = rM(i,i)/rK(i,i);
            }
        }

        // For the computation of large systems, suitable distances between the unity entries
        // are needed in the starting iteration vectors in r. This is taken into account via
        // variable l.
        int nd = nn/nc;
        int l = nn - nd;

        for (int j = 1; j< nc ; j++)
        {
            rt = 0.0;
            for (int i = 0; i< l ; i++)
            {
                if (w(i) >= rt)
                {
                    rt = w(i);
                    ij = i;
                }
            }
            for (int i = l-1; i< nn ; i++)
            {
                if (w(i) > rt)
                {
                    rt = w(i);
                    ij = i;
                }
            }
            tt(j) = double(ij); //gets the equation-nrs of excited DOFs
            w(ij) = 0.0;

            l -= nd;

            r(ij,j) = 1.0;
        }
        if (echo_level > 1) {
            std::cout << "Degrees of freedom excited by unit starting iteration vectors:" << std::endl;
            for (int j = 1; j < nc; j++)
                std::cout << tt(j) << std::endl;
        }
        // r has now the initial values for the eigenvectors

        // TODO add random numbers to last column of r(:,nc) -> Bathe SSP00131

        if (echo_level>1) std::cout << "Initialize linear solver." << std::endl;
        mpLinearSolver->InitializeSolutionStep(rK, temp_tt, tt);

        //------------------------------------------------------------------------------------
        // start of iteration loop
        //------------------------------------------------------------------------------------
        nite = 0;
        do // label 100
        {
            nite ++;

            if (echo_level > 1) std::cout << "Subspace Iteration Eigenvalue Solver: Iteration no. " << nite <<std::endl;

            // compute projection ar of matrix rK
            for (int j = 0; j< nc; j++)
            {
                for (int k = 0; k< nn; k++) // tt gets an iteration eigenvector
                    tt(k) = r(k,j);

                // K*temp_tt = tt
                if (echo_level > 1) std::cout << "Backsubstitute using initialized linear solver." << std::endl;
                mpLinearSolver->PerformSolutionStep(rK, temp_tt, tt);
                for(int k=0; k<nn; k++)
                    tt(k) = temp_tt(k);

                for (int i = j; i< nc; i++) //SSP00161
                {
                    art = 0.;
                    for (int k = 0; k< nn; k++)
                        art += r(k,i)*tt(k);
                    ar(j,i) = art;
                    ar(i,j) = art;
                }
                for (int k = 0; k< nn; k++)
                    r(k,j) = tt(k);   // (r = xbar)
            }

            // compute projection br of matrix
            for (int j = 0; j < nc; j++) //SSP00170
            {
                for (int k = 0; k < nn ; k++)  // temp gets an iteration eigenvector
                    temp(k)=r(k,j);

                // tt = rM*temp
                TSparseSpaceType::Mult(rM, temp, tt);

                for (int i = j; i < nc; i++)
                {
                    brt = 0.;
                    for (int k=0; k<nn; k++)
                        brt += r(k,i)*tt(k);
                    br(j,i) = brt;
                    br(i,j) = brt;
                }                         // label 180

                for (int k=0; k<nn; k++)
                    r(k,j)=tt(k);

            }                           // label 160

            mpEigenValueSolver->Solve(ar,br,eigv,vec);

            if(!eigen_solver_successful)
            {
                //TODO how can this be checked?
                std::cout << "Eigen solution was not successful!" << std::endl;
                break;
            }

            // sorting eigenvalues according to their values
            do
            {
                is = 0; // label 350
                for (int i = 0; i< nc1; i++)
                {
                    if ((eigv(i+1)) < (eigv(i)))
                    {
                        is++;
                        eigvt = eigv(i+1);
                        eigv(i+1) = eigv(i) ;
                        eigv(i) = eigvt ;
                        for (int k = 0; k < nc; k++)
                        {
                            rt = vec(k,i+1) ;
                            vec(k,i+1) = vec(k,i) ;
                            vec(k,i)   = rt ;
                        }
                    }
                }                         // label 360
            } while (is != 0);

            // compute eigenvectors
            for (int i = 0; i< nn; i++)  // label 375
            {
                for (int j = 0; j< nc; j++)
                    tt(j) = r(i,j);
                for (int k = 0; k< nc; k++)
                {
                    rt = 0.0;
                    for (int j = 0; j< nc; j++)
                        rt += tt(j)*vec(j,k) ;
                    r(i,k) = rt ;
                }
            }                           // label 420   (r = z)

            // convergency check alternative A: change in eigenvalues
            for (int i = 0; i< nc; i++)
            {
                dif = (eigv(i) - prev_eigv(i));
                rtolv(i) = std::fabs(dif / eigv(i));
            }

            // convergency check alternative B: according to BATHE (more restrictive)
            // for (int i = 0; i<nc; i++)
            // {
            //     double vdot = 0.0;
            //     for (int j=0; j<nc; j++)
            //         vdot += vec(j,i) * vec(j,i);
            //     double eigv2 = eigv(i) * eigv(i);
            //     dif = vdot-eigv2;
            //     rtolv(i) = std::sqrt(std::max(dif, 1e-24*eigv2)/vdot);
            // }

            bool converged = true;
            for (int i = 0; i< nroot; i++)
            {
                if (rtolv(i) > rtol) converged = false;
            }
            if (converged)
            {
                if (echo_level > 0) std::cout << "Convergence reached after " << nite << " iterations within a relative tolerance: " << rtol << std::endl;
                break;
            }
            else{
                if (nite >= nitem)
                {
                    if (echo_level > 0) std::cout << "Convergence not reached in " << nitem << " iterations." << std::endl;
                    break;
                }
            }

            for (int i = 0; i< nc ; i++)
                prev_eigv(i) = eigv(i); // label 410 and 440

            //TODO evaluate if necessary: GramSchmidt Orthogonalization
        } while (true);

        mpLinearSolver->FinalizeSolutionStep(rK, temp_tt, tt);

        if(eigen_solver_successful)
        {
            // compute eigenvectors
            for (int j = 0; j< nc; j++)
            {
                for (int k = 0; k< nn; k++)
                    tt(k) = r(k,j);

                // K*temp_tt = tt
                mpLinearSolver->PerformSolutionStep(rK, temp_tt, tt);

                for (int k = 0; k< nn; k++)
                    r(k,j) = temp_tt(k);
            }

            // copy results to function parameters
            if (static_cast<int>(rEigenvalues.size()) != nroot)
                rEigenvalues.resize(nroot);
            if (static_cast<int>(rEigenvectors.size1()) != nroot || static_cast<int>(rEigenvectors.size2()) != nn)
                rEigenvectors.resize(nroot, nn);

            for (int i = 0; i< nroot; i++)
            {
                rEigenvalues(i) = eigv(i);
                for (int j=0; j< nn; j++)
                    rEigenvectors(i,j) = r(j,i);
            }

            if (echo_level > 0) KRATOS_WATCH(rEigenvalues);

            // make vectors M-normalized
            if (mParam["mass_normalization"].GetBool())
            {
                MassNormalizeEigenVectors(rM, rEigenvectors);
            }

            // orient vectors
            if (mParam["orient_eigen_vectors"].GetBool()){
                OrientEigenVectors(rEigenvectors);
            }
            // TODO
            // 1. speed up
            // 2. sturm sequence check
        }
        else{
            KRATOS_ERROR << "Solution failed!" <<std::endl;
        }

        if (echo_level > 0) {
            std::cout  << "SubspaceIterationEigenvalueSolver completed in " <<
            OpenMPUtils::GetCurrentTime() - start_time << " seconds" << std::endl;
        }

        return;
    }

    /**
     * This method returns directly the first eigen value obtained
     * @param rK: The stiffness matrix
     * @param rM: The mass matrix
     * @return The first eigenvalue
     */
    double GetEigenValue(
        SparseMatrixType& rK,
        SparseMatrixType& rM
        )
    {
        VectorType eigen_values;
        DenseMatrixType eigen_vectors;

        Solve(rK, rM, eigen_values, eigen_vectors);

        return eigen_values[0];
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

    /**
     * Print information about this object.
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "SubspaceIterationEigenvalueSolver.";
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

    Parameters mParam;
    typename TLinearSolverType::Pointer mpLinearSolver;
    typename TLinearSolverType::Pointer mpEigenValueSolver;

    ///@}
    ///@name Private Operations
    ///@{
    void MassNormalizeEigenVectors(SparseMatrixType& rM, DenseMatrixType& rEigenvectors){
        VectorType tmp(rM.size1());
        VectorType vec_i(rM.size1());
        int size1 = int(rEigenvectors.size1());
        int size2 = int(rEigenvectors.size2());
        for (int i=0; i<size1; ++i)
        {
            for (int j=0; j<size2; ++j)
            {
                vec_i(j) = rEigenvectors(i,j);
            }
            // tmp = rM*vec_i
            TSparseSpaceType::Mult(rM, vec_i, tmp);
            double phi_t_M_phi = TSparseSpaceType::Dot(tmp, vec_i);
            double factor = 1.0 / std::sqrt(std::fabs(phi_t_M_phi));

            for (int j=0; j<size2; ++j)
            {
                rEigenvectors(i,j) *= factor;
            }
        }
    }

    void OrientEigenVectors(DenseMatrixType& rEigenvectors){
        int size1 = int(rEigenvectors.size1());
        int size2 = int(rEigenvectors.size2());
        if (size2 == 0) return;
        for (int i=0; i<size1; ++i)
        {
            if (rEigenvectors(i,0) >= 0.0)
                continue;
            for (int j=0; j<size2; ++j)
            {
                rEigenvectors(i,j) *= -1.0;
            }
        }
    }
    ///@}

}; // class SubspaceIterationEigenvalueSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType, class TReordererType>
inline std::istream& operator >>(std::istream& rIStream,
        SubspaceIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType,  TLinearSolverType, TReordererType>& rThis) {
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TLinearSolverType, class TReordererType>
inline std::ostream& operator <<(std::ostream& rOStream,
        const SubspaceIterationEigenvalueSolver<TSparseSpaceType, TDenseSpaceType, TLinearSolverType, TReordererType>& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos

#endif // defined(KRATOS_SUBSPACE_ITERATION_EIGEN_VALUE_SOLVER_H_INCLUDED)