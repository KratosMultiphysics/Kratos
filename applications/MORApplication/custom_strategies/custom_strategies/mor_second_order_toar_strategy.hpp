//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//

#if !defined(MOR_SECOND_ORDER_TOAR_STRATEGY)
#define MOR_SECOND_ORDER_TOAR_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/builtin_timer.h"
#include "custom_strategies/custom_strategies/mor_offline_second_order_strategy.hpp"
#include "custom_utilities/orthogonalization_utility.hpp"

//default builder and solver
#include "custom_strategies/custom_builder_and_solvers/system_matrix_builder_and_solver.hpp"

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
 * @class MorSecondOrderTOARStrategy
 * @ingroup KratosCore
 * @brief This is the linear MOR matrix output strategy
 * @details This strategy builds the K and M matrices and outputs them
 * @author Aditya Ghantasala
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
          class TReducedSparseSpace,
          class TReducedDenseSpace
          >
class MorSecondOrderTOARStrategy
    // : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    : public MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace >
{
    using complex = std::complex<double>;

  public:
    ///@name Type Definitions
    ///@{
    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(MorSecondOrderTOARStrategy);

    // typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef MorOfflineSecondOrderStrategy< TSparseSpace, TDenseSpace, TLinearSolver, TReducedSparseSpace, TReducedDenseSpace > BaseType;

    // typedef SystemMatrixBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::MatrixType TDenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType TDenseMatrixPointerType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef TLinearSolver TLinearSolverType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename TReducedSparseSpace::MatrixType TReducedSparseMatrixType;

    typedef typename TReducedSparseSpace::VectorType TReducedSparseVectorType;

    typedef typename TReducedDenseSpace::MatrixType TReducedDenseMatrixType;

    typedef typename TReducedDenseSpace::VectorType TReducedDenseVectorType;


    ///@}
    ///@name Life Cycle

    ///@{

    /**
     * Default constructor for the damped case
     * @param rModelPart The model part of the problem
     * @param pScheme The integration schemed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MorSecondOrderTOARStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename BaseType::TBuilderAndSolverType::Pointer pBuilderAndSolver,
        typename TLinearSolverType::Pointer pNewLinearSolver,
        // typename TLinearSolverType::Pointer pNewSecondLinearSolver,
        std::vector<complex> samplingPoints,
        std::vector<int> orders)
        : BaseType(rModelPart, pScheme, pBuilderAndSolver, pNewLinearSolver, true),
            // mpSecondLinearSolver(pNewSecondLinearSolver),
            mSamplingPoints(samplingPoints), mOrders(orders)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(samplingPoints.size() != orders.size()) << "Vectors of sampling points and orders must have the same length." << std::endl;
        this->mTolerance = 1e-12;

        KRATOS_CATCH("");
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~MorSecondOrderTOARStrategy() override
    {
        this->Clear();
    }


    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (this->mSolutionStepIsInitialized == false)
        {
            BaseType::InitializeSolutionStep();

            std::size_t reduced_system_size = 0;
            for( std::size_t i=0; i<mOrders.size(); ++i ) {
                reduced_system_size += mOrders[i];
            }

            const unsigned int system_size = this->GetSystemSize();

            TReducedDenseSpace::Resize(this->GetKr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetDr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetMr(), reduced_system_size, reduced_system_size);
            TReducedDenseSpace::Resize(this->GetBasis(), system_size, reduced_system_size);


            this->mSolutionStepIsInitialized = true;
        }

        KRATOS_CATCH("");

    }


    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;
        std::cout << "hello! this is where the second order TOAR MOR magic happens" << std::endl;

        TSystemMatrixType& r_K = this->GetK();
        TSystemMatrixType& r_Ki = this->GetKi();
        TSystemMatrixType& r_M = this->GetM();
        TSystemMatrixType& r_Mi = this->GetMi();
        TSystemMatrixType& r_C = this->GetC();
        TSystemVectorType& r_RHS = this->GetSystemVector();
        TReducedDenseMatrixType& r_basis = this->GetBasis();

        const std::size_t system_size = r_K.size1();
        const complex z(0,1);

        // std::cout << "1\n";
        // TReducedSparseMatrixType r_A(system_size, system_size);
        TReducedSparseMatrixType r_K_cplx(system_size, system_size);
        noalias(r_K_cplx) = r_K + z * TReducedSparseMatrixType(r_Ki);
        // std::cout << "2\n";
        TReducedSparseMatrixType r_M_cplx(system_size, system_size);
        noalias(r_M_cplx) = r_M + z * TReducedSparseMatrixType(r_Mi);
        // std::cout << "3\n";
        TReducedSparseMatrixType r_C_cplx(system_size, system_size);
        noalias(r_C_cplx) = r_C;
        // std::cout << "4\n";
        TReducedSparseVectorType r_RHS_tmp = TReducedSparseVectorType(r_RHS);
        TReducedSparseVectorType r_dx(system_size);

        // working matrices and vectors
        TReducedSparseMatrixType tmp_Kt(system_size, system_size);
        TReducedSparseMatrixType tmp_Ct(system_size, system_size);
        TReducedSparseVectorType p(system_size);
        TReducedSparseVectorType q(system_size);
        TReducedDenseMatrixType V;
        TReducedDenseMatrixType RP;
        TReducedDenseMatrixType RQ;
        TReducedDenseMatrixType H;
        TReducedDenseVectorType x;
        TReducedDenseVectorType tmp_vec;
        std::size_t rk = 0;
        std::size_t nrk = 0;
        double alpha = 0;
        double alpha1 = 0;
        double beta = 0;
        complex tmp(0.0, 0.0);
        TReducedDenseVectorType h;
        TReducedDenseVectorType sp;
        TReducedDenseVectorType sq;

        complex sampling_point(0.0, 0.0);
        complex sampling_point2(0.0, 0.0);
        int basis_index = 0;
        std::size_t order;

        std::cout << "begin loop\n";
        for( std::size_t i=0; i<mSamplingPoints.size(); ++i ) {
            sampling_point = mSamplingPoints[i];
            sampling_point2 = std::pow(sampling_point, 2);
            order = static_cast<std::size_t>(mOrders[i]);

            // resize working matrices
            if( V.size1() != system_size || V.size2() != order ) {
                V.resize(system_size, order, false);
            }
            if( RP.size1() != order + 1 || RP.size2() != order + 1 ) {
                RP.resize(order + 1, order + 1);
                noalias(RP) = ZeroMatrix(order + 1, order + 1);
            }
            if( RQ.size1() != order + 1 || RQ.size2() != order + 1 ) {
                RQ.resize(order + 1, order + 1);
                noalias(RQ) = ZeroMatrix(order + 1, order + 1);
            }
            if( H.size1() != order + 1 || H.size2() != order ) {
                H.resize(order + 1, order, false);
                noalias(H) = ZeroMatrix(order + 1, order);
            }

            // build matrices
            noalias(tmp_Kt) = sampling_point2 * r_M_cplx + sampling_point * r_C_cplx + r_K_cplx;
            noalias(tmp_Ct) = 2.0 * sampling_point * r_M_cplx + r_C_cplx;
            std::cout << "finished build\n";

            // decompose tmp_Kt
            this->mpLinearSolver->InitializeSolutionStep(tmp_Kt, r_dx, r_RHS_tmp);
            std::cout << "finished decompose\n";

            // compute p (= v0)
            this->mpLinearSolver->PerformSolutionStep(tmp_Kt, p, r_RHS_tmp);
            KRATOS_WATCH(tmp_Kt)
            KRATOS_WATCH(r_RHS_tmp)
            std::cout << "finished solve\n";
            KRATOS_WATCH(p)
            // KRATOS_WATCH(prod(tmp_Kt, p))
            // this->mpLinearSolver->Solve(tmp_Kt, p, r_RHS_tmp);
            // KRATOS_WATCH(p)
            KRATOS_WATCH(norm_2(p))
            p /= norm_2(p);
            KRATOS_WATCH(p)
            std::cout << "finished norm\n";
            noalias(q) = ZeroVector(system_size);

            noalias(column(V, 0)) = p;
            std::cout << "finished column\n";
            RP(0, 0) = 1.0;
            rk = 1;

            // main loop
            // for( std::size_t j=0; j<2; ++j ) {
            for( std::size_t j=0; j<order-1; ++j ) {
                // KRATOS_WATCH(tmp_Ct)
                // KRATOS_WATCH(p)
                // KRATOS_WATCH(r_dx)
                axpy_prod( tmp_Ct, p, r_dx, true );
                KRATOS_WATCH(r_dx)
                std::cout << "finished product 1\n";
                this->mpLinearSolver->PerformSolutionStep(tmp_Kt, p, r_dx);
                std::cout << "finished solve 1\n";

                axpy_prod( r_M_cplx, q, r_dx, true );
                std::cout << "finished product 2\n";
                this->mpLinearSolver->PerformSolutionStep(tmp_Kt, q, r_dx);
                // KRATOS_WATCH(q)
                std::cout << "finished solve 2\n";

                p -= q;
                p *= -1.0;
                std::cout << "finished add\n";
                alpha1 = norm_2(p);
                std::cout << "finished norm\n";

                // KRATOS_WATCH(p)
                // KRATOS_WATCH(q)
                // KRATOS_WATCH(alpha)
                //% level-one orthogonalization: MGS
                //  x = zeros(rk,1);
                //  for jj = 1: rk
                //      x(jj) = V(:,jj)'*p;
                //      p = p -x(jj)*V(:,jj);
                //  end

                // level one orthogonalization: MGS
                x.resize(rk, false);
                for( std::size_t k=0; k<rk; ++k ) {
                    x[k] = inner_prod(conj(column(V,k)), p);
                    // complex y(0,0);
                    // auto v = column(V,k);
                    // for( std::size_t kk=0; kk<p.size(); ++kk ) {
                    //     y += v[kk] * p[kk];
                    // }
                    // KRATOS_WATCH(column(V,k))
                    // KRATOS_WATCH(p)
                    KRATOS_WATCH(x)
                    // KRATOS_WATCH(y)
                    p -= x[k] * column(V,k);
                }
                std::cout << "finished MGS1\n";
                KRATOS_WATCH(p)

                // % re-orth
                // if norm(p) < 0.717*alpha1
                //     for jj = 1: rk
                //         tmp = V(:,jj)'*p;
                //         p = p -tmp*V(:,jj);
                //         x(jj) = tmp + x(jj);
                //     end
                // end

                // reorthogonalization
                if( norm_2(p) < 0.717*alpha1 ) {
                    std::cout << "enter reorthogonalization!" << std::endl;
                    for( std::size_t k=0; k<rk; ++k ) {
                        tmp = inner_prod(conj(column(V,k)), p);
                        p -= tmp * column(V,k);
                        x[k] += tmp;
                    }
                }
                std::cout << "finished reorth\n";
                KRATOS_WATCH(p)
                KRATOS_WATCH(x)

                alpha = norm_2(p);

                    //  if alpha>tol
                    //     nrk = rk +1;
                    //     V(:,nrk) = p/alpha;
                    // else
                    //     nrk = rk;
                    //     disp('BREAK');
                    // end
                // KRATOS_WATCH(V)
                if( alpha > mTolerance ) {
                    nrk = rk + 1;
                    // nrk = rk;
                    KRATOS_WATCH(nrk)
                    KRATOS_WATCH(alpha)
                    noalias(column(V, nrk-1)) = p / alpha;
                } else {
                    nrk = rk;
                    // nrk = rk - 1;
                    KRATOS_WARNING("TOAR strategy") << "BREAK" << std::endl;
                }
                std::cout << "finished new column\n";
                KRATOS_WATCH(V)

                // % level-two orthogonalization: MGS
                // h =  zeros(i,1);
                // sp = x;
                // sq = RP(1:rk,i);
                // alpha1 = norm([sp;sq]);

                // level two orthogonalization: MGS
                // KRATOS_WATCH(RP)
                // KRATOS_WATCH(rk)
                // KRATOS_WATCH(j)
                // KRATOS_WATCH( MatrixVectorRange(RP, range(0,1), range(0,1)) )
                // MatrixVectorRange aa(RP, range(0, rk), range(j, j+1));
                std::cout << "hier?\n";
                KRATOS_WATCH(j)
                KRATOS_WATCH(rk)
                h.resize(j+1, false);
                noalias(h) = ZeroVector(j+1);
                KRATOS_WATCH(h)
                sp.resize(rk, false);
                noalias(sp) = x;
                KRATOS_WATCH(sp)
                KRATOS_WATCH(x)
                sq.resize(rk, false);
                KRATOS_WATCH(sq)
                KRATOS_WATCH(RP)
                // KRATOS_WATCH(ComplexMatrixVectorRange(RP, range(0, 1), range(1, 3)))
                // KRATOS_WATCH(ComplexMatrixVectorRange(RP, range(0, rk), range(j, j+1)))
                // noalias(sq) = ComplexMatrixVectorRange(RP, range(0, rk), range(j, j+1));
                // noalias(sq) = ComplexMatrixVectorRange(RP, range(0, rk), range(j, j+1));
                noalias(sq) = ComplexMatrixVectorSlice(RP, slice(0, 1, rk), slice(j, 0, rk));

                // two-norm of vector [sp; sq]
                alpha1 = 0.0;
                for( std::size_t k=0; k<rk; ++k ) {
                    alpha1 += std::pow(std::real(sp[k]), 2);
                    alpha1 += std::pow(std::imag(sp[k]), 2);
                    alpha1 += std::pow(std::real(sq[k]), 2);
                    alpha1 += std::pow(std::imag(sq[k]), 2);
                }
                // for( std::size_t k=0; k<rk; ++k ) {
                //     alpha1 += std::pow(sq[k], 2);
                // }
                alpha1 = std::sqrt(alpha1);

                KRATOS_WATCH(h)
                KRATOS_WATCH(sp)
                KRATOS_WATCH(sq)
                KRATOS_WATCH(alpha1)

                // for jj = 1:i
                //     h(jj) = RP(1:rk,jj)'*sp+RQ(1:rk,jj)'*sq;
                //     sp = sp - RP(1:rk,jj)*h(jj);
                //     sq = sq - RQ(1:rk,jj)*h(jj);
                // end
                for( std::size_t k=0; k<j+1; ++k ) {
                    ComplexMatrixVectorSlice tmp_range_p(RP, slice(0, 1, rk), slice(k, 0, rk));
                    ComplexMatrixVectorSlice tmp_range_q(RQ, slice(0, 1, rk), slice(k, 0, rk));
                    // ComplexMatrixVectorSlice tmp_range_q(RQ, range(0, rk), range(k, k+1));
                    h(k) = inner_prod(conj(tmp_range_p), sp);
                    h(k) += inner_prod(conj(tmp_range_q), sq);
                    sp -= tmp_range_p * h(k);
                    sq -= tmp_range_q * h(k);
                }
                std::cout << "finished MGS2\n";
                KRATOS_WATCH(h)
                KRATOS_WATCH(sp)
                KRATOS_WATCH(sq)

                // % re-orth
                // if norm([sp;sq]) < 0.717* alpha1
                //     for jj = 1:i
                //         tmp = RP(1:rk,jj)'*sp+RQ(1:rk,jj)'*sq;
                //         sp = sp - RP(1:rk,jj)*tmp;
                //         sq = sq - RQ(1:rk,jj)*tmp;
                //         h(jj) = tmp + h(jj);
                //     end
                // end

                // reorthogonalization 2
                beta = 0.0;
                for( std::size_t k=0; k<rk; ++k ) {
                    beta += std::pow(std::real(sp[k]), 2);
                    beta += std::pow(std::imag(sp[k]), 2);
                    beta += std::pow(std::real(sq[k]), 2);
                    beta += std::pow(std::imag(sq[k]), 2);
                }
                beta = std::sqrt(beta);

                if( beta < 0.717 * alpha1 ) {
                    for( std::size_t k=0; k<j+1; ++k ) {
                        // ComplexMatrixVectorRange tmp_range_p(RP, range(0, rk), range(k, k+1));
                        // ComplexMatrixVectorRange tmp_range_q(RQ, range(0, rk), range(k, k+1));
                        ComplexMatrixVectorSlice tmp_range_p(RP, slice(0, 1, rk), slice(k, 0, rk));
                        ComplexMatrixVectorSlice tmp_range_q(RQ, slice(0, 1, rk), slice(k, 0, rk));
                        tmp = inner_prod(conj(tmp_range_p), sp);
                        tmp += inner_prod(conj(tmp_range_q), sq);
                        sp -= tmp_range_p * tmp;
                        sq -= tmp_range_q * tmp;
                        h(k) += tmp;
                    }
                }
                std::cout << "finished reorth after MGS2\n";
                KRATOS_WATCH(h)
                KRATOS_WATCH(sp)
                KRATOS_WATCH(sq)

                // beta = norm([sp;sq;alpha]);
                // if beta>tol
                //     sp = sp/beta;
                //     sq = sq/beta;
                //     alpha = alpha/beta;
                // else
                //     disp('break');
                // end

                // check for deflation
                beta = std::pow(alpha, 2);
                for( std::size_t k=0; k<rk; ++k ) {
                    beta += std::pow(std::real(sp[k]), 2);
                    beta += std::pow(std::imag(sp[k]), 2);
                    beta += std::pow(std::real(sq[k]), 2);
                    beta += std::pow(std::imag(sq[k]), 2);
                }
                beta = std::sqrt(beta);

                if( mTolerance < beta ) {
                    sp /= beta;
                    sq /= beta;
                    alpha /= beta;
                } else {
                    KRATOS_WARNING("TOAR strategy") << "BREAK (tol > beta)" << std::endl;
                }
                std::cout << "finished beta\n";
                KRATOS_WATCH(alpha)
                KRATOS_WATCH(beta)
                KRATOS_WATCH(sp)
                KRATOS_WATCH(sq)

                // update matrices
                std::cout << "begin update matrices\n";
                KRATOS_WATCH(RP)
                KRATOS_WATCH(RQ)
                // RP(1:rk,i+1) = sp;
                // RQ(1:rk,i+1) = sq;
                // RP(rk+1,i+1)   = alpha;
                // H(1:i+1,i) = [h;beta];
                // p = V(:,1:rk+1)*[sp;alpha];
                // q = V(:,1:rk)*sq;
                // rk = nrk;

                // noalias(ComplexMatrixVectorRange(RP, range(0, rk), range(j+1, j+2))) = sp;
                // noalias(ComplexMatrixVectorRange(RQ, range(0, rk), range(j+1, j+2))) = sq;
                noalias(ComplexMatrixVectorSlice(RP, slice(0, 1, rk), slice(j+1, 0, rk))) = sp;
                noalias(ComplexMatrixVectorSlice(RQ, slice(0, 1, rk), slice(j+1, 0, rk))) = sq;
                RP(rk, j+1) = alpha;
                // noalias(ComplexMatrixVectorRange(H, range(0, j+1), range(j, j+1))) = h;
                noalias(ComplexMatrixVectorSlice(H, slice(0, 1, j+1), slice(j, 0, j+1))) = h;
                H(j+1, j) = alpha;

                // ComplexMatrixVectorRange tmp_range_p(V, range(0, system_size), range(0, rk+1));
                // sp.push_back(alpha);
                // axpy_prod(tmp_range_p, sp, p);
                KRATOS_WATCH(rk)
                noalias(p) = column(V, rk) * alpha;
                noalias(q) = ZeroVector(system_size);
                // KRATOS_WATCH(p)
                for( std::size_t k=0; k<sp.size(); ++k ) {
                    p += column(V, k) * sp[k];
                    q += column(V, k) * sq[k];
                }
                std::cout << "finished update matrices\n";
                KRATOS_WATCH(RP)
                KRATOS_WATCH(RQ)
                KRATOS_WATCH(H)
                KRATOS_WATCH(p)
                KRATOS_WATCH(q)
                rk = nrk;
            }

            // update overall basis
            for( std::size_t j=0; j<order; ++j ) {
                noalias(column(r_basis, basis_index + j)) = column(V, j);
            }
            basis_index += order;
        }

        std::cout << "FINISHED TOAR!!!!!!!!!\n";
        KRATOS_WATCH(V)
        KRATOS_WATCH(r_basis)

        //TODO reorthogonalize!

		return true;

        KRATOS_CATCH("");
    }

    virtual void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        TSystemMatrixType& r_K = this->GetK();
        TSystemMatrixType& r_Ki = this->GetKi();
        TSystemMatrixType& r_M = this->GetM();
        TSystemMatrixType& r_Mi = this->GetMi();
        TSystemMatrixType& r_C = this->GetC();
        TSystemVectorType& r_RHS = this->GetSystemVector();
        TSystemVectorType& r_OV = this->GetOutputVector();
        auto& r_Kr = this->GetKr();
        auto& r_Cr = this->GetDr();
        auto& r_Mr = this->GetMr();
        auto& r_force_vector_reduced = this->GetRHSr();
        auto& r_output_vector_r = this->GetOVr();
        TReducedDenseMatrixType& r_basis = this->GetBasis();

        const std::size_t reduced_system_size = r_basis.size2();
        const std::size_t system_size = r_K.size1();
        const complex z(0,1);

        TReducedSparseMatrixType r_K_cplx(system_size, system_size);
        noalias(r_K_cplx) = r_K + z * TReducedSparseMatrixType(r_Ki);
        TReducedSparseMatrixType r_M_cplx(system_size, system_size);
        noalias(r_M_cplx) = r_M + z * TReducedSparseMatrixType(r_Mi);
        TReducedSparseMatrixType r_C_cplx(system_size, system_size);
        noalias(r_C_cplx) = r_C;
        // TReducedSparseVectorType r_RHS_cplx = TReducedSparseVectorType(r_RHS);
        // TReducedSparseVectorType r_OV_cplx = TReducedSparseVectorType(r_OV);

        this->template ProjectMatrix<TReducedSparseMatrixType>(r_K_cplx, r_basis, r_Kr);
        this->template ProjectMatrix<TReducedSparseMatrixType>(r_C_cplx, r_basis, r_Cr);
        this->template ProjectMatrix<TReducedSparseMatrixType>(r_M_cplx, r_basis, r_Mr);

        r_force_vector_reduced.resize( reduced_system_size, false);
        axpy_prod(r_RHS, conj(r_basis), r_force_vector_reduced, true);
        //TODO dont use conj()!

        r_output_vector_r.resize( reduced_system_size, false);
        r_output_vector_r = prod( r_OV, r_basis );

        KRATOS_WATCH(r_force_vector_reduced)
        KRATOS_WATCH(r_output_vector_r)

        KRATOS_CATCH("");
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
        return "MorSecondOrderTOARStrategy";
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

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

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    bool mUseComplexFlag;

    // typename TLinearSolverType::Pointer mpSecondLinearSolver;

    std::vector<complex> mSamplingPoints;
    std::vector<int> mOrders;

    double mTolerance;

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

    MorSecondOrderTOARStrategy(const MorSecondOrderTOARStrategy &Other){};

    ///@}

}; /* Class MorSecondOrderTOARStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */

#endif /* MOR_SECOND_ORDER_TOAR_STRATEGY  defined */