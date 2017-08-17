//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $StructuralMechanicsApplication $
//   Last modified by:    $Author: quirin.aumann@tum.de   $
//   Date:                $Date:            August 2017   $
//   Revision:            $Revision:                0.0   $

#if !defined(KRATOS_HARMONIC_ANALYSIS_STRATEGY )
#define  KRATOS_HARMONIC_ANALYSIS_STRATEGY

// System includes
#include<iostream>
#include<vector>
#include<iterator>

// External includes
#include<boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/strategies/solving_strategy.h"

// Application includes
#include "structural_mechanics_application_variables.h"

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

/// Strategy for solving generalized eigenvalue problems.
template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver
         >
class HarmonicAnalysisStrategy
    : public SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(HarmonicAnalysisStrategy);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType::Pointer SchemePointerType;

    typedef typename BaseType::TBuilderAndSolverType::Pointer BuilderAndSolverPointerType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    typedef std::complex<double> ComplexType;

    typedef boost::numeric::ublas::vector<ComplexType> ComplexVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    HarmonicAnalysisStrategy(
        ModelPart& model_part,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpBuilderAndSolver = pBuilderAndSolver;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        mInitializeWasPerformed = false;

        SparseVectorType* AuxForceVector = new SparseVectorType;
        mpForceVector = boost::shared_ptr<SparseVectorType>(AuxForceVector);

        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build only once)
        this->SetRebuildLevel(0);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    HarmonicAnalysisStrategy(const HarmonicAnalysisStrategy& Other) = delete;

    /// Destructor.
    virtual ~HarmonicAnalysisStrategy()
    {
        // Clear() controls order of deallocation to avoid invalid memory access
        // in some special cases.
        // warning: BaseType::GetModelPart() may be invalid here.
        this->Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetIsInitialized(bool val)
    {
        mInitializeWasPerformed = val;
    }

    bool GetIsInitialized() const
    {
        return mInitializeWasPerformed;
    }

    SparseVectorType& GetForceVector()
    {
        return *mpForceVector;
    }

    SparseVectorPointerType& pGetForceVector()
    {
        return mpForceVector;
    }

    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    };

    SchemePointerType& pGetScheme()
    {
        return mpScheme;
    };

    void SetBuilderAndSolver(BuilderAndSolverPointerType pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    BuilderAndSolverPointerType& pGetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        this->pGetBuilderAndSolver()->SetReshapeMatrixFlag(flag);
    }

    bool GetReformDofSetAtEachStepFlag() const
    {
        return this->pGetBuilderAndSolver()->GetReshapeMatrixFlag();
    }

    /// Set verbosity level of the solving strategy.
    /**
     * - 0 -> mute... no echo at all
     * - 1 -> print time and basic information
     * - 2 -> print linear solver data
     * - 3 -> print debug information
     */
    void SetEchoLevel(int Level)
    {
        BaseType::SetEchoLevel(Level);
        this->pGetBuilderAndSolver()->SetEchoLevel(Level);
    }

    /// Initialization to be performed once before using the strategy.
    virtual void Initialize()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Initialize() of HarmonicAnalysisStrategy." << std::endl;

        this->Check();

        auto& pScheme = this->pGetScheme();

        if (pScheme->SchemeIsInitialized() == false)
            pScheme->Initialize(rModelPart);

        if (pScheme->ElementsAreInitialized() == false)
            pScheme->InitializeElements(rModelPart);

        if (pScheme->ConditionsAreInitialized() == false)
            pScheme->InitializeConditions(rModelPart);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Initialize() of HarmonicAnalysisStrategy." << std::endl;

        // set up the system
        auto& pBuilderAndSolver = this->pGetBuilderAndSolver();

        // Reset solution dofs
        boost::timer system_construction_time;
        // Set up list of dofs
        boost::timer setup_dofs_time;
        pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;
        }

        // Set global equation ids
        boost::timer setup_system_time;
        pBuilderAndSolver->SetUpSystem(rModelPart);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "setup_system_time : " << setup_system_time.elapsed() << std::endl;
        }

        // initialize the force vector; this does not change during the computation
        auto& pForceVector = this->pGetForceVector();
        auto& rForceVector = *pForceVector;
        const unsigned int system_size = pBuilderAndSolver->GetEquationSystemSize();
        
        boost::timer force_vector_build_time;
        if (rForceVector.size() != system_size)
            rForceVector.resize(system_size, false);
        pBuilderAndSolver->BuildRHS(pScheme,rModelPart,rForceVector);
        
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "force_vector_build_time : " << force_vector_build_time.elapsed() << std::endl;
        }

        // get the damping coefficients
        auto& rProcessInfo = rModelPart.GetProcessInfo();
        // KRATOS_WATCH(rProcessInfo)
        if( rProcessInfo.Has(RAYLEIGH_ALPHA) )
            mRayleighAlpha = rProcessInfo[RAYLEIGH_ALPHA];

        if( rProcessInfo.Has(RAYLEIGH_BETA) )
            mRayleighBeta = rProcessInfo[RAYLEIGH_BETA];
        
        KRATOS_CATCH("")
    }

    double Solve()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();

        // operations to be done once
        if (this->GetIsInitialized() == false)
        {
            Initialize();
            this->SetIsInitialized(true);
        }

        this->InitializeSolutionStep();

        auto excitation_frequency = rModelPart.GetProcessInfo()[TIME];

        // get eigenvalues and eigenvectors
        DenseVectorType eigenvalues = rModelPart.GetProcessInfo()[EIGENVALUE_VECTOR];
        DenseMatrixType eigenvectors = rModelPart.GetProcessInfo()[EIGENVECTOR_MATRIX];

        const unsigned int n_dofs = eigenvectors.size2();
        const unsigned int n_modes = eigenvalues.size();
        
        auto f = this->GetForceVector();

        ComplexVectorType mode_weight;
        mode_weight.resize(n_modes, false);
        mode_weight = ZeroVector( n_modes );

        double modal_damping = 0.0;

        for( std::size_t i = 0; i < n_modes; ++i )
        {
            modal_damping = mRayleighAlpha / (2 * eigenvalues[i]) + mRayleighBeta * eigenvalues[i] / 2;
            // rows are columns and vice-versa
            ComplexType factor( eigenvalues[i] - pow( excitation_frequency, 2.0 ), 2 * modal_damping * std::sqrt(eigenvalues[i]) * excitation_frequency );
            mode_weight[i] = inner_prod( row( eigenvectors, i ), f ) / factor;
        }

        ComplexVectorType modal_displacement;
        modal_displacement.resize(n_dofs, false);
        modal_displacement = ZeroVector( n_dofs );

        for( std::size_t i = 0; i < n_dofs; ++i )
        {
            for( std::size_t j = 0; j < n_modes; ++j)
            {
                modal_displacement[i] = modal_displacement[i] + mode_weight[j] * eigenvectors(j,i);
            }
        }

        this->AssignVariables(modal_displacement);
        this->FinalizeSolutionStep();

        return 0.0;

        KRATOS_CATCH("")
    }

    /// Clear the strategy.
    virtual void Clear()
    {
        KRATOS_TRY

        // if the preconditioner is saved between solves, it should be cleared here
        auto& pBuilderAndSolver = this->pGetBuilderAndSolver();
        pBuilderAndSolver->GetLinearSystemSolver()->Clear();

        if (this->pGetForceVector() != nullptr)
            this->pGetForceVector() = nullptr;


        // re-setting internal flag to ensure that the dof sets are recalculated
        pBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        pBuilderAndSolver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;
        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;

        KRATOS_CATCH("")
    }

    /// Initialization to be performed before every solve.
    virtual void InitializeSolutionStep()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering InitializeSolutionStep() of HarmonicAnalysisStrategy" << std::endl;

        auto& pBuilderAndSolver = this->pGetBuilderAndSolver();
        auto& pScheme = this->pGetScheme();
        auto& pForceVector = this->pGetForceVector();
        auto& rForceVector = *pForceVector;

        // // initialize dummy vectors
        auto pA = SparseSpaceType::CreateEmptyMatrixPointer();
        auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rA = *pA;
        auto& rDx = *pDx;

        SparseSpaceType::Resize(rA,SparseSpaceType::Size(rForceVector),SparseSpaceType::Size(rForceVector));
        SparseSpaceType::SetToZero(rA);
        SparseSpaceType::Resize(rDx,SparseSpaceType::Size(rForceVector));
        SparseSpaceType::Set(rDx,0.0);

        // initial operations ... things that are constant over the solution step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,rForceVector);

        // initial operations ... things that are constant over the solution step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,rForceVector);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting InitializeSolutionStep() of HarmonicAnalysisStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /// Check whether initial input is valid.
    virtual int Check()
    {
        KRATOS_TRY

        auto& rModelPart = BaseType::GetModelPart();
        const auto rank = rModelPart.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Check() of HarmonicAnalysisStrategy" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(rModelPart);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(rModelPart);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Check() of HarmonicAnalysisStrategy" << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

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

    SchemePointerType mpScheme;

    BuilderAndSolverPointerType mpBuilderAndSolver;

    bool mInitializeWasPerformed;

    SparseVectorPointerType mpForceVector;

    double mRayleighAlpha;

    double mRayleighBeta;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Assign the modal displacement to the displacement dofs
    void AssignVariables(ComplexVectorType& rModalDisplacement, int step=0)
    {
        auto& rModelPart = BaseType::GetModelPart();
        for( auto itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++ )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = itNode->GetDofs();
            
            for( auto itDof = std::begin(rNodeDofs); itDof != std::end(rNodeDofs); itDof++ )
            {
                if( !itDof->IsFixed() )
                {
                    itDof->GetSolutionStepValue(step) = std::abs(rModalDisplacement(itDof->EquationId()));
                }
                else
                {
                    itDof->GetSolutionStepValue(step) = 0.0;
                }
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}

}; /* Class HarmonicAnalysisStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos */

#endif /* KRATOS_HARMONIC_ANALYSIS_STRATEGY  defined */

