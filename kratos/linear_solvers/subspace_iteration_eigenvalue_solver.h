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
            "verbosity": 1,
            "linear_solver_settings": {},
            "eigen_sub_solver_settings" :{}
        })");

        // don't validate linear_solver_settings and eigen_sub_solver_settings here
        mParam.ValidateAndAssignDefaults(default_params);

        BaseType::SetTolerance( mParam["tolerance"].GetDouble() );
        BaseType::SetMaxIterationsNumber( mParam["max_iteration"].GetInt() );

        KRATOS_ERROR_IF( mParam["linear_solver_settings"]["solver_type"].GetString() == "eigen_pardiso_llt") <<
            "eigen_pardiso_llt cannot handle negative entries on the main diagonal" << std::endl;

        std::cout << "\nWARNING: Make sure the linear solver can handle negative entries on the main diagonal!\n" << std::endl;

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
     * K is a symmetric matrix. M is a symmetric positive-definite matrix.
     */
    void Solve(
            SparseMatrixType& rK,
            SparseMatrixType& rM,
            VectorType& rEigenvalues,
            DenseMatrixType& rEigenvectors) override
    {
        std::cout << "Start solving for eigen values."  << std::endl;

        const int verbosity = mParam["verbosity"].GetInt();
        int nroot = mParam["number_of_eigenvalues"].GetInt(); // number of eigenvalues requested
        int nitem = BaseType::GetMaxIterationsNumber();
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
                std::cout << tt(j-1) << std::endl;
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

            // convergency check alternative A: change in eigenvalues
            for (int i = 1; i<= nc; i++)
            {
                dif = (eigv(i-1) - prev_eigv(i-1));
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

            for (int i = 1; i<= nroot; i++)
            {
                if (rtolv(i-1) > rtol) goto label400 ;
            }

            std::cout << "Convergence reached after " << nite << " iterations within a relative tolerance: " << rtol << std::endl;

            break;

            label400:      // "not converged so far"
            if (nite >= nitem)
            {
                std::cout << "Convergence not reached in " << nitem << " iterations." << std::endl;
                break;
            }

            for (int i = 1; i<= nc ; i++)
                prev_eigv(i-1) = eigv(i-1); // label 410 and 440

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
            rEigenvalues.resize(nroot);
            rEigenvectors.resize(nroot, nn);
            for (int i = 0; i< nroot; i++)
            {
                rEigenvalues(i) = eigv(i);
                for (int j=0; j< nn; j++)
                    rEigenvectors(i,j) = r(j,i);
            }

            KRATOS_WATCH(rEigenvalues);

        // 1. make vectors M-normalized
        if (mParam["mass_normalization"].GetBool())
        {
            MassNormalizeEigenVectors(rM, rEigenvectors);
        }

        // 2. orient vectors
        if (mParam["orient_eigen_vectors"]){
            OrientEigenVectors(rEigenvectors);
        }
        // TODO
        // 3. speed up
        // 4. sturm sequence check
        }
        else{
            KRATOS_ERROR << "Solution failed!" <<std::endl;
        }

        return;
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
    MassNormalizeEigenVectors(SparseMatrixType& rM, DenseMatrixType& rEigenvectors){
        // GlobalVectorBasis* aux = this->create_GlobalVector(nn);
        // cfloat phi_t_M_phi, multiplier;
        // for (cint i=0; i< _NROOT; ++i)
        // {
        //     aux->be_Product_Of(_mtx_B, _eigen_vecs[i]);
        //     phi_t_M_phi = aux->dot_Product_With(_eigen_vecs[i]);
        //     multiplier =  1.0 / sqrt(fabs(phi_t_M_phi));
        //     _eigen_vecs[i]->scale_With(multiplier);
        // }
    }

    OrientEigenVectors(DenseMatrixType& rEigenvectors){
        if (rEigenvalues.size2() == 0) return;
        for (int i=0; i<rEigenvectors.size1())
        {
            if (rEigenvalues(i,0) > 0.0)
                continue;
            for (int j=0; j<rEigenvectors.size2())
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