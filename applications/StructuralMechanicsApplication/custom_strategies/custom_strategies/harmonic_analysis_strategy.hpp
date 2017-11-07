//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
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

    typedef TDenseSpace DenseSpaceType;

    typedef typename TDenseSpace::VectorType DenseVectorType;

    typedef typename TDenseSpace::MatrixType DenseMatrixType;

    typedef typename TDenseSpace::MatrixPointerType DenseMatrixPointerType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorType SparseVectorType;
    
    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;

    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

    typedef std::complex<double> ComplexType;

    typedef boost::numeric::ublas::vector<ComplexType> ComplexVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    HarmonicAnalysisStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        BuilderAndSolverPointerType pBuilderAndSolver,
        bool UseMaterialDampingFlag = false
        )
        : SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart)
    {
        KRATOS_TRY

        mpScheme = pScheme;

        mpBuilderAndSolver = pBuilderAndSolver;

        // ensure initialization of system matrices in InitializeSolutionStep()
        mpBuilderAndSolver->SetDofSetIsInitializedFlag(false);

        mInitializeWasPerformed = false;

        mpForceVector = SparseSpaceType::CreateEmptyVectorPointer();
        mpModalMatrix = DenseSpaceType::CreateEmptyMatrixPointer();

        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;
        mSystemDamping = 0.0;
        this->SetUseMaterialDampingFlag(UseMaterialDampingFlag);

        // default echo level (mute)
        this->SetEchoLevel(0);

        // default rebuild level (build only once)
        this->SetRebuildLevel(0);

        KRATOS_CATCH("")
    }

    /// Deleted copy constructor.
    HarmonicAnalysisStrategy(const HarmonicAnalysisStrategy& Other) = delete;

    /// Destructor.
    ~HarmonicAnalysisStrategy() override
    {
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

    void SetUseMaterialDampingFlag(bool flag)
    {
        mUseMaterialDamping = flag;
    }

    bool GetUseMaterialDampingFlag() const
    {
        return mUseMaterialDamping;
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

        auto& r_model_part = BaseType::GetModelPart();
        auto& r_process_info = r_model_part.GetProcessInfo();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Initialize() of HarmonicAnalysisStrategy." << std::endl;

        this->Check();

        auto& p_scheme = this->pGetScheme();

        if (p_scheme->SchemeIsInitialized() == false)
            p_scheme->Initialize(r_model_part);

        if (p_scheme->ElementsAreInitialized() == false)
            p_scheme->InitializeElements(r_model_part);

        if (p_scheme->ConditionsAreInitialized() == false)
            p_scheme->InitializeConditions(r_model_part);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting Initialize() of HarmonicAnalysisStrategy." << std::endl;

        // set up the system
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();

        // Reset solution dofs
        boost::timer system_construction_time;
        // Set up list of dofs
        boost::timer setup_dofs_time;
        p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;
        }

        // Set global equation ids
        boost::timer setup_system_time;
        p_builder_and_solver->SetUpSystem(r_model_part);
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "setup_system_time : " << setup_system_time.elapsed() << std::endl;
        }

        // initialize the force vector; this does not change during the computation
        auto& r_force_vector = *mpForceVector;
        const unsigned int system_size = p_builder_and_solver->GetEquationSystemSize();
        
        boost::timer force_vector_build_time;
        if (r_force_vector.size() != system_size)
            r_force_vector.resize(system_size, false);
        r_force_vector = ZeroVector( system_size );
        p_builder_and_solver->BuildRHS(p_scheme,r_model_part,r_force_vector);
        
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "force_vector_build_time : " << force_vector_build_time.elapsed() << std::endl;
        }

        // initialize the modal matrix
        auto& r_modal_matrix = *mpModalMatrix;
        const std::size_t n_modes = r_process_info[EIGENVALUE_VECTOR].size();
        if( r_modal_matrix.size1() != system_size || r_modal_matrix.size2() != n_modes )
            r_modal_matrix.resize( system_size, n_modes, false );
        r_modal_matrix = ZeroMatrix( system_size, n_modes );

        boost::timer modal_matrix_build_time;
        for( std::size_t i = 0; i < n_modes; ++i )
        {
            for( auto& node : r_model_part.Nodes() )
            {
                ModelPart::NodeType::DofsContainerType node_dofs = node.GetDofs();
                const std::size_t n_node_dofs = node_dofs.size();
                const Matrix& r_node_eigenvectors = node.GetValue(EIGENVECTOR_MATRIX);

                if( node_dofs.IsSorted() == false )
                {
                    node_dofs.Sort();
                }

                for( std::size_t j = 0; j < n_node_dofs; ++j )
                {
                    const auto it_dof = std::begin(node_dofs) + j;
                    r_modal_matrix(it_dof->EquationId(), i) = r_node_eigenvectors(i, j);
                }
            }
        }

        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "modal_matrix_build_time : " << modal_matrix_build_time.elapsed() << std::endl;
        }

        // get the damping coefficients if they exist
        for( auto& property : r_model_part.PropertiesArray() )
        {
            if( property->Has(SYSTEM_DAMPING_RATIO) )
            {
                mSystemDamping = property->GetValue(SYSTEM_DAMPING_RATIO);
            }
            
            if( property->Has(RAYLEIGH_ALPHA) && property->Has(RAYLEIGH_BETA) )
            {
                mRayleighAlpha = property->GetValue(RAYLEIGH_ALPHA);
                mRayleighBeta = property->GetValue(RAYLEIGH_BETA);
            }
        }

        // compute the effective material damping if required
        if( mUseMaterialDamping )
        {
            // throw an error, if no submodelparts are present
            if( r_model_part.NumberOfSubModelParts() < 1 )
                KRATOS_ERROR << "No submodelparts detected!" << std::endl;
            
            //initialize all required variables
            r_model_part.GetProcessInfo()[BUILD_LEVEL] = 2;
            mMaterialDampingRatios = ZeroVector( n_modes );
            
            //initialize dummy vectors
            auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
            auto pb = SparseSpaceType::CreateEmptyVectorPointer();
            auto& rDx = *pDx;
            auto& rb = *pb;
            SparseSpaceType::Resize(rDx,system_size);
            SparseSpaceType::Set(rDx,0.0);
            SparseSpaceType::Resize(rb,system_size);
            SparseSpaceType::Set(rb,0.0);

            //loop over all modes and initialize the material damping ratio per mode
            boost::timer material_damping_build_time;
        
            for( std::size_t i = 0; i < n_modes; ++i )
            {
                double up = 0.0;
                double down = 0.0;
                auto modal_vector = column( r_modal_matrix, i );
                for( auto& sub_model_part : r_model_part.SubModelParts() )
                {
                    double damping_coefficient = 0.0;
                    for( auto& property : sub_model_part.PropertiesArray() )
                    {
                        if( property->Has(SYSTEM_DAMPING_RATIO) )
                        {
                            damping_coefficient = property->GetValue(SYSTEM_DAMPING_RATIO);
                        }
                    }
                    
                    //initialize the submodelpart stiffness matrix
                    auto temp_stiffness_matrix = SparseSpaceType::CreateEmptyMatrixPointer();
                    p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, 
                        temp_stiffness_matrix,
                        pDx,
                        pb,
                        r_model_part.Elements(),
                        r_model_part.Conditions(),
                        r_model_part.GetProcessInfo());

                    //build stiffness matrix for submodelpart material
                    p_builder_and_solver->BuildLHS(p_scheme, sub_model_part, *temp_stiffness_matrix);

                    //compute strain energy of the submodelpart and the effective damping ratio
                    double strain_energy = 0.5 * inner_prod( prod(modal_vector, *temp_stiffness_matrix), modal_vector );
                    down += strain_energy;
                    up += damping_coefficient * strain_energy;
                }
                if( down < std::numeric_limits<double>::epsilon() )
                    KRATOS_ERROR << "No valid effective material damping ratio could be computed. Are all elements to be damped available in the submodelparts? Are the modal vectors available? " << std::endl;
                
                mMaterialDampingRatios(i) = up / down;
            }

            if (BaseType::GetEchoLevel() > 0 && rank == 0)
            {
                std::cout << "modal_matrix_build_time : " << material_damping_build_time.elapsed() << std::endl;
                KRATOS_WATCH(mMaterialDampingRatios)
            }
        }
        
        KRATOS_CATCH("")
    }

    double Solve()
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();

        // operations to be done once
        if (this->GetIsInitialized() == false)
        {
            Initialize();
            this->SetIsInitialized(true);
        }

        this->InitializeSolutionStep();

        auto& r_process_info = r_model_part.GetProcessInfo();
        double excitation_frequency = r_process_info[TIME];

        // get eigenvalues
        DenseVectorType eigenvalues = r_process_info[EIGENVALUE_VECTOR];
        const std::size_t n_modes = eigenvalues.size();

        // DenseMatrixType eigenvectors;
        const std::size_t n_dofs = this->pGetBuilderAndSolver()->GetEquationSystemSize();
        
        auto& f = *mpForceVector;

        ComplexType mode_weight;
        ComplexVectorType modal_displacement;
        modal_displacement.resize(n_dofs, false);
        modal_displacement = ZeroVector( n_dofs );

        double modal_damping = 0.0;

        for( std::size_t i = 0; i < n_modes; ++i )
        {

            modal_damping = mSystemDamping + mRayleighAlpha / (2 * eigenvalues[i]) + mRayleighBeta * eigenvalues[i] / 2;
            
            if( mUseMaterialDamping )
            {
                modal_damping += mMaterialDampingRatios[i];
            }

            auto& r_modal_matrix = *mpModalMatrix;

            DenseVectorType modal_vector(n_dofs);
            TDenseSpace::GetColumn(i, r_modal_matrix, modal_vector);

            ComplexType factor( eigenvalues[i] - pow( excitation_frequency, 2.0 ), 2 * modal_damping * std::sqrt(eigenvalues[i]) * excitation_frequency );
            mode_weight = inner_prod( modal_vector, f ) / factor;

            // compute the modal displacement as a superposition of modal_weight * eigenvector
            for( auto& node : r_model_part.Nodes() )
            {
                auto& node_dofs = node.GetDofs();
                const std::size_t n_node_dofs = node_dofs.size();
                const Matrix& r_node_eigenvectors = node.GetValue(EIGENVECTOR_MATRIX);

                if (node_dofs.IsSorted() == false)
                {
                    node_dofs.Sort();
                }

                for (std::size_t j = 0; j < n_node_dofs; j++)
                {
                    auto it_dof = std::begin(node_dofs) + j;
                    modal_displacement[it_dof->EquationId()] = modal_displacement[it_dof->EquationId()] + mode_weight * r_node_eigenvectors(i,j);
                }
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
        auto& p_builder_and_solver = this->pGetBuilderAndSolver();
        p_builder_and_solver->GetLinearSystemSolver()->Clear();

        SparseSpaceType::Clear(mpForceVector);
        DenseSpaceType::Clear(mpModalMatrix);

        // re-setting internal flag to ensure that the dof sets are recalculated
        p_builder_and_solver->SetDofSetIsInitializedFlag(false);

        p_builder_and_solver->Clear();

        this->pGetScheme()->Clear();

        mInitializeWasPerformed = false;
        mUseMaterialDamping = false;
        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;
        mSystemDamping = 0.0;

        KRATOS_CATCH("")
    }

    /// Initialization to be performed before every solve.
    virtual void InitializeSolutionStep()
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering InitializeSolutionStep() of HarmonicAnalysisStrategy" << std::endl;

        BuilderAndSolverPointerType& p_builder_and_solver = this->pGetBuilderAndSolver();
        SchemePointerType& p_scheme = this->pGetScheme();
        auto& r_force_vector = *mpForceVector;

        // // initialize dummy vectors
        auto pA = SparseSpaceType::CreateEmptyMatrixPointer();
        auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rA = *pA;
        auto& rDx = *pDx;

        SparseSpaceType::Resize(rA,SparseSpaceType::Size(r_force_vector),SparseSpaceType::Size(r_force_vector));
        SparseSpaceType::SetToZero(rA);
        SparseSpaceType::Resize(rDx,SparseSpaceType::Size(r_force_vector));
        SparseSpaceType::Set(rDx,0.0);

        // initial operations ... things that are constant over the solution step
        p_builder_and_solver->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,r_force_vector);

        // initial operations ... things that are constant over the solution step
        p_scheme->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,r_force_vector);

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Exiting InitializeSolutionStep() of HarmonicAnalysisStrategy" << std::endl;

        KRATOS_CATCH("")
    }

    /// Check whether initial input is valid.
    virtual int Check()
    {
        KRATOS_TRY

        auto& r_model_part = BaseType::GetModelPart();
        const auto rank = r_model_part.GetCommunicator().MyPID();

        if (BaseType::GetEchoLevel() > 2 && rank == 0)
            std::cout << "Entering Check() of HarmonicAnalysisStrategy" << std::endl;

        // check the model part
        BaseType::Check();

        // check the scheme
        this->pGetScheme()->Check(r_model_part);

        // check the builder and solver
        this->pGetBuilderAndSolver()->Check(r_model_part);

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

    DenseMatrixPointerType mpModalMatrix;

    double mRayleighAlpha;

    double mRayleighBeta;

    double mSystemDamping;

    bool mUseMaterialDamping;

    vector< double > mMaterialDampingRatios;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Assign the modal displacement to the dofs and the phase angle to the reaction
    void AssignVariables(ComplexVectorType& rModalDisplacement, int step=0)
    {
        auto& r_model_part = BaseType::GetModelPart();
        for( auto& node : r_model_part.Nodes() )
        {
            ModelPart::NodeType::DofsContainerType& rNodeDofs = node.GetDofs();
            
            for( auto it_dof = std::begin(rNodeDofs); it_dof != std::end(rNodeDofs); it_dof++ )
            {
                if( !it_dof->IsFixed() )
                {
                    //absolute displacement
                    it_dof->GetSolutionStepValue(step) = std::abs(rModalDisplacement(it_dof->EquationId()));
                    //phase angle
                    it_dof->GetSolutionStepReactionValue(step) = std::abs(std::arg(rModalDisplacement(it_dof->EquationId())));
                }
                else
                {
                    it_dof->GetSolutionStepValue(step) = 0.0;
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

