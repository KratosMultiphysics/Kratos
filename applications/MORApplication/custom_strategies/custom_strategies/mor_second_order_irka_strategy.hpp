//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Quirin Aumann
//                   Matthias Ebert
//

#if !defined(MOR_SECOND_ORDER_IRKA_STRATEGY)
#define MOR_SECOND_ORDER_IRKA_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"
// #include "utilities/qr_utility.h"
// #include "utilities/openmp_utils.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"

//default builder and solver
// #include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

// includes for linear solver factory
#include "factories/linear_solver_factory.h"
#include "includes/kratos_parameters.h"

// include the eigensolver
// #include "custom_solvers/gen_eigensystem_solver.h"
#include "custom_utilities/generalized_eigenvalue_utility.h"
#include "custom_utilities/orthogonalization_utility.hpp"
#include "custom_utilities/complex_sort_utility.hpp"
// #include "../EigenSolversApplication/custom_solvers/eigen_direct_solver.h"

#include "omp.h"

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
 * @class MorSecondOrderIRKAStrategy
 * @ingroup KratosMORApplication
 * @brief This is the Iterative Rational Krylov Algorithm
 * @details This strategy builds a reduction basis according to
 *      Wyatt, Sarah Alice. Issues in interpolatory model reduction:
 *      Inexact solves, second-order systems and DAEs. Diss. Virginia Tech, 2012.
 * @author Matthias Ebert, based on code of Quirin Aumann
 */
template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver,
          class TReducedSparseSpace,
          class TReducedDenseSpace,
          bool TUseModalDamping
          >
class MorSecondOrderIRKAStrategy
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace >
{
  
    using complex = std::complex<double>;

  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderIRKAStrategy);

    typedef TUblasSparseSpace<complex> ComplexSparseSpaceType;

    typedef MorOfflineSecondOrderStrategy<TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef TLinearSolver TLinearSolverType;

    typedef typename TDenseSpace::VectorType TDenseVectorType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename ComplexSparseSpaceType::MatrixType ComplexSparseMatrixType;

    typedef typename ComplexSparseSpaceType::MatrixPointerType ComplexSparseMatrixPointerType;

    typedef typename ComplexSparseSpaceType::VectorType ComplexSparseVectorType;

    typedef typename ComplexSparseSpaceType::VectorPointerType ComplexSparseVectorPointerType;

    typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedSparseSpace::VectorType TReducedSparseVectorType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;

    typedef typename TReducedDenseSpace::VectorType TReducedDenseVectorType;

    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderIRKAStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        vector< std::complex<double> > SamplingPoints,
        size_t MaxIter,
        double Tolerance,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true, MoveMeshFlag), 
            mMaxIter(MaxIter), mTolerance(Tolerance)
    {
        KRATOS_TRY;

        // create complex conjugates of the sampling points and sort them
        const size_t n_sampling_points = SamplingPoints.size();
        mSamplingPoints = ComplexZeroVector(2*n_sampling_points);
        for( size_t i=0; i<n_sampling_points; ++i )
        {
            mSamplingPoints(i) = SamplingPoints(i);
            mSamplingPoints(n_sampling_points+i) = std::conj(SamplingPoints(i));
        }
        ComplexSortUtility::PairComplexConjugates(mSamplingPoints);

        KRATOS_ERROR_IF( mMaxIter < 1 ) << "Invalid number of maximal iterations provided\n";
        KRATOS_ERROR_IF( (mTolerance >= 1.) || (mTolerance < 0.) ) << "Invalid tolerance provided\n";

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderIRKAStrategy() override
    {

        if (this->mpM != nullptr)
            TSparseSpace::Clear(this->mpM);
        if (this->mpA != nullptr)
            TSparseSpace::Clear(this->mpA);
        if (this->mpS != nullptr)
            TSparseSpace::Clear(this->mpS);
        // if (this->mp != nullptr)
        //     TSparseSpace::Clear(mpb);


        this->Clear();
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();

            const size_t reduced_system_size = mSamplingPoints.size();
            const size_t system_size = this->GetSystemSize();

            TReducedDenseSpace::Resize(this->GetKr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetDr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetMr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetBasis(), system_size, reduced_system_size);

            this->mSolutionStepIsInitialized = true;
        }


        KRATOS_CATCH("");
    }

    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        return 0;

        KRATOS_CATCH("")

    }
    

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        std::cout << "hello! this is where the second order IRKA MOR magic happens" << std::endl;
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename BaseType::TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        const int rank = BaseType::GetModelPart().GetCommunicator().MyPID();

        // get the full system size matrices                
        TSystemMatrixType& r_K = this->GetSystemMatrix();
        TSystemMatrixType& r_M = this->GetMassMatrix();
        TSystemMatrixType& r_D = this->GetDampingMatrix();
        TSystemVectorType& r_RHS = this->GetSystemVector();
        TReducedDenseMatrixType& r_basis = this->GetBasis();

        // get the reduced matrices
        TReducedDenseMatrixType& r_Kr = this->GetKr();
        TReducedDenseMatrixType& r_Mr = this->GetMr();
        TReducedDenseMatrixType& r_Dr = this->GetDr();

        // set the system sizes
        const size_t system_size = this->GetSystemSize();
        const size_t n_sampling_points = mSamplingPoints.size();
        const size_t reduced_system_size = n_sampling_points;

        // copy mass matrix and rhs to the complex space
        ComplexSparseMatrixType r_M_tmp = ComplexSparseMatrixType(r_M);
        ComplexSparseVectorType r_RHS_tmp = ComplexSparseVectorType(r_RHS);

        // create a complex stiffness matrix
        TReducedSparseMatrixType r_K_cplx;
        if( TUseModalDamping )
        {
            r_K_cplx = TReducedSparseMatrixType(r_D);
            complex z(0,1);
            r_K_cplx *= reinterpret_cast<typename TReducedSparseSpace::DataType(&)[2]>(z)[0];   //cast complex to double (real part)
            r_K_cplx += r_K;

            noalias(r_Dr) = ZeroMatrix(r_Dr.size1(), r_Dr.size2());
        }

        // create solution vector
        ComplexSparseVectorPointerType tmp_dx = ComplexSparseSpaceType::CreateEmptyVectorPointer();
        ComplexSparseVectorType& r_tmp_dx   = *tmp_dx;
        ComplexSparseSpaceType::Resize(r_tmp_dx, system_size);
        ComplexSparseSpaceType::Set(r_tmp_dx, 0.0);

        // create dynamic stiffness matrix
        auto kdyn = ComplexSparseSpaceType::CreateEmptyMatrixPointer();
        auto& r_kdyn   = *kdyn;
        ComplexSparseSpaceType::Resize(r_kdyn, system_size, system_size); // n x n

        BuiltinTimer irka_overall_time;

        //initial basis
        for( size_t i=0; i<n_sampling_points/2; ++i )
        {
            if( TUseModalDamping )
                noalias(r_kdyn) = r_K_cplx;
            else
            {
                noalias(r_kdyn) = r_D;
                r_kdyn *= mSamplingPoints(2*i);
                r_kdyn += r_K;
            }
            r_kdyn += std::pow( mSamplingPoints(2*i), 2.0 ) * r_M_tmp;
            this->mpLinearSolver->Solve( r_kdyn, r_tmp_dx, r_RHS_tmp );

            KRATOS_DEBUG_ERROR_IF(isinf(norm_2(r_tmp_dx))) << "No solution could be obtained (norm infinity)!";

            column(r_basis, 2*i)   = real(r_tmp_dx);
            column(r_basis, 2*i+1) = imag(r_tmp_dx);
        }

        OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);

        // create loop variables
        size_t iter = 0;
        double error = 1;
        vector<complex> eigenvalues;
        auto samplingPoints_old = mSamplingPoints;
        double projection_time;
        double construction_time;

        while( iter < mMaxIter && error > mTolerance )
        {
            BuiltinTimer irka_iteration_time;
            BuiltinTimer irka_projection_time;

            // project onto reduced space
            if( TUseModalDamping )
            {
                this->template ProjectMatrix<TReducedSparseMatrixType>(r_K_cplx, r_basis, r_Kr);
                this->template ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_Mr);
            }
            else
            {
                this->template ProjectMatrix<TSystemMatrixType>(r_K, r_basis, r_Kr);
                this->template ProjectMatrix<TSystemMatrixType>(r_D, r_basis, r_Dr);
                this->template ProjectMatrix<TSystemMatrixType>(r_M, r_basis, r_Mr);
            }
            projection_time = irka_projection_time.ElapsedSeconds();

            if( TUseModalDamping )
            {
                // compute generalized eigenvalues
                vector<complex> tmp_eigenvalues;
                GeneralizedEigenvalueUtility::Compute<TReducedDenseSpace>(r_Kr, r_Mr, tmp_eigenvalues);

                std::for_each(tmp_eigenvalues.begin(), tmp_eigenvalues.end(),
                    [](complex &a) {a = std::sqrt(a);});

                eigenvalues.resize(2*reduced_system_size);
                noalias(subrange(eigenvalues, 0, reduced_system_size)) = tmp_eigenvalues;
                std::for_each(tmp_eigenvalues.begin(), tmp_eigenvalues.end(), 
                    [](complex &a) {a = std::conj(a);});
                noalias(subrange(eigenvalues, reduced_system_size, reduced_system_size*2)) = tmp_eigenvalues;
            }
            else
            {
                // compute generalized eigenvalues
                GeneralizedEigenvalueUtility::ComputePolynomial<TReducedDenseSpace>(r_Kr, r_Dr, r_Mr, eigenvalues);

                // use mirror images of first r eigenvalues as expansion points
                eigenvalues *= -1.;
            }

            ComplexSortUtility::PairComplexConjugates(eigenvalues);
            noalias(mSamplingPoints) = subrange(eigenvalues, 0, reduced_system_size);

            // calculate error
            error = norm_2(mSamplingPoints - samplingPoints_old);
            error /= norm_2(samplingPoints_old);

            //update basis
            BuiltinTimer basis_construction_time;
            for( size_t i=0; i<n_sampling_points/2; ++i )
            {
                if( TUseModalDamping )
                    noalias(r_kdyn) = r_K_cplx;
                else
                {
                    noalias(r_kdyn) = r_D;
                    r_kdyn *= mSamplingPoints(2*i);
                    r_kdyn += r_K;
                }
                r_kdyn += std::pow( mSamplingPoints(2*i), 2.0 ) * r_M_tmp;
                ComplexSparseSpaceType::SetToZero(r_tmp_dx);
                
                this->mpLinearSolver->Solve( r_kdyn, r_tmp_dx, r_RHS_tmp );

                KRATOS_DEBUG_ERROR_IF(isinf(norm_2(r_tmp_dx))) << "No solution could be obtained (norm infinity)!";

                column(r_basis, 2*i)   = real(r_tmp_dx);
                column(r_basis, 2*i+1) = imag(r_tmp_dx);
            }
            
            // orthogonalize basis
            OrthogonalizationUtility::OrthogonalizeQR<TReducedDenseSpace>(r_basis);

            construction_time = basis_construction_time.ElapsedSeconds();
            iter++;
            samplingPoints_old = mSamplingPoints;

            KRATOS_INFO_IF("IRKA error after iteration " + iter, BaseType::GetEchoLevel() > 0 && rank == 0)
                << ": e=" << error << "\n";
            KRATOS_INFO_IF("IRKA Iteration Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << irka_iteration_time.ElapsedSeconds() << "\n";
            KRATOS_INFO_IF("\tProjection Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << projection_time << "\n";
            KRATOS_INFO_IF("\tBasis update Time", BaseType::GetEchoLevel() > 0 && rank == 0)
                << construction_time << "\n";
        }

        KRATOS_INFO_IF("IRKA Complete Solve Time", BaseType::GetEchoLevel() > 0 && rank == 0)
            << irka_overall_time.ElapsedSeconds() << "\n";

		return true;

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep() override
    {
        // matrices have been projected inside the IRKA loop
        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_output_vector = this->GetOutputVector();
        auto& r_output_vector_r = this->GetOVr();
        auto& r_basis = this->GetBasis();

        const size_t reduced_system_size = r_basis.size2();

        BuiltinTimer system_projection_time;

        r_force_vector_reduced.resize( reduced_system_size, false);
        r_force_vector_reduced = prod( this->GetSystemVector(), r_basis );

        r_output_vector_r.resize( reduced_system_size, false);
        r_output_vector_r = prod( r_output_vector, r_basis );
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "MorSecondOrderIRKAStrategy";
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
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

  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    vector< complex > mSamplingPoints; //vector of currently used sampling points
    size_t mMaxIter; //maximum iterations
    double mTolerance; // tolerance

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
    ///@name Un accessible methods

    ///@{

    /**
     * Copy constructor.
     */

    MorSecondOrderIRKAStrategy(const MorSecondOrderIRKAStrategy &Other){};

    ///@}

}; /* Class MorSecondOrderIRKAStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_SECOND_ORDER_IRKA_STRATEGY  defined */