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

    typedef typename TDenseSpace::MatrixPointerType DenseMatrixPointerType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename TSparseSpace::VectorPointerType SparseVectorPointerType;

    typedef typename TSparseSpace::VectorType SparseVectorType;

    typedef typename TSparseSpace::MatrixType SparseMatrixType;
    typedef typename TSparseSpace::MatrixPointerType SparseMatrixPointerType;

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

        DenseMatrixType* AuxModalMatrix = new DenseMatrixType;
        mpModalMatrix = boost::shared_ptr<DenseMatrixType>(AuxModalMatrix);

        mRayleighAlpha = 0.0;
        mRayleighBeta = 0.0;
        // KRATOS_WATCH(mMaterialDampingRatios.size())
        // mMaterialDampingRatios[0] = 0.01;
        // mMaterialDampingRatios[1] = 0.05;

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

    DenseMatrixType& GetModalMatrix()
    {
        return *mpModalMatrix;
    }

    DenseMatrixPointerType& pGetModalMatrix()
    {
        return mpModalMatrix;
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

        auto& rModelPart = BaseType::GetModelPart();
        auto& rProcessInfo = rModelPart.GetProcessInfo();
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
        rForceVector = ZeroVector( system_size );
        pBuilderAndSolver->BuildRHS(pScheme,rModelPart,rForceVector);
        
        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "force_vector_build_time : " << force_vector_build_time.elapsed() << std::endl;
        }

        // initialize the modal matrix
        auto& pModalMatrix = this->pGetModalMatrix();
        auto& rModalMatrix = *pModalMatrix;
        const std::size_t n_modes = rProcessInfo[EIGENVALUE_VECTOR].size();
        if( rModalMatrix.size1() != system_size || rModalMatrix.size2() != n_modes )
            rModalMatrix.resize( system_size, n_modes, false );
        rModalMatrix = ZeroMatrix( system_size, n_modes );

        boost::timer modal_matrix_build_time;
        for( std::size_t i = 0; i < n_modes; ++i )
        {
            for( ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); itNode++ )
            {
                ModelPart::NodeType::DofsContainerType& node_dofs = itNode->GetDofs();
                const std::size_t n_node_dofs = node_dofs.size();
                Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);

                if( node_dofs.IsSorted() == false )
                {
                    node_dofs.Sort();
                }

                for( std::size_t j = 0; j < n_node_dofs; ++j )
                {
                    auto itDof = std::begin(node_dofs) + j;
                    rModalMatrix(itDof->EquationId(), i) = rNodeEigenvectors(i, j);
                }
            }
        }

        if (BaseType::GetEchoLevel() > 0 && rank == 0)
        {
            std::cout << "modal_matrix_build_time : " << modal_matrix_build_time.elapsed() << std::endl;
        }

        // get the damping coefficients

        if( rProcessInfo.Has(RAYLEIGH_ALPHA) )
            mRayleighAlpha = rProcessInfo[RAYLEIGH_ALPHA];

        if( rProcessInfo.Has(RAYLEIGH_BETA) )
            mRayleighBeta = rProcessInfo[RAYLEIGH_BETA];

        this->SetUseMaterialDampingFlag(true);
        if( mUseMaterialDamping )
        {
            //initialize all required variables
            rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
            SparseMatrixType* AuxStiffnessMatrix;
            SparseMatrixPointerType temp_stiffness_matrix;
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
                // double modal_material_damping = 0.0;
                double up = 0.0;
                double down = 0.0;
                Vector modal_vector = column( rModalMatrix, i );
                for( ModelPart::SubModelPartIterator itSubModelPart = rModelPart.SubModelPartsBegin(); itSubModelPart!= rModelPart.SubModelPartsEnd(); itSubModelPart++ )
                {
                    // KRATOS_WATCH(itSubModelPart->Name())
                    auto current_properties = itSubModelPart->rProperties();
                    double damping_coefficient = 0.0;
                    // KRATOS_WATCH(current_properties)
                    for( ModelPart::PropertiesIterator itProperty = itSubModelPart->PropertiesBegin(); itProperty != itSubModelPart->PropertiesEnd(); itProperty++ )
                    {
                        // std::cout << "hey" << std::endl;
                        // KRATOS_WATCH(*itProperty)
                        // KRATOS_WATCH(*itProperty.Has(SYSTEM_DAMPING_RATIO))
                        if( itProperty->Has(SYSTEM_DAMPING_RATIO) )
                        {
                            // std::cout << "haha" << std::endl;
                            damping_coefficient = itProperty->GetValue(SYSTEM_DAMPING_RATIO);
                        }
                    }
                    
                    //initialize the submodelpart stiffness matrix
                    AuxStiffnessMatrix = new SparseMatrixType;
                    temp_stiffness_matrix = boost::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);
                    pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, 
                        temp_stiffness_matrix,
                        pDx,
                        pb,
                        rModelPart.Elements(),
                        rModelPart.Conditions(),
                        rModelPart.GetProcessInfo());

                    //build stiffness matrix for submodelpart material
                    pBuilderAndSolver->BuildLHS(pScheme,*itSubModelPart,*temp_stiffness_matrix);
                    double strain_energy = 0.5 * inner_prod( prod(modal_vector, *temp_stiffness_matrix), modal_vector );
                    down += strain_energy;
                    up += damping_coefficient * strain_energy;
                    // auto test = prod(*temp_stiffness_matrix,modal_vector);
                    // KRATOS_WATCH(damping_coefficient)
                    // KRATOS_WATCH(modal_vector)
                    // KRATOS_WATCH(*temp_stiffness_matrix)
                    // KRATOS_WATCH(down)
                    // KRATOS_WATCH(up)
                }

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

        auto& rModelPart = BaseType::GetModelPart();

        // operations to be done once
        if (this->GetIsInitialized() == false)
        {
            Initialize();
            this->SetIsInitialized(true);
        }

        this->InitializeSolutionStep();

        auto& rProcessInfo = rModelPart.GetProcessInfo();
        double excitation_frequency = rProcessInfo[TIME];

        // get eigenvalues
        DenseVectorType eigenvalues = rProcessInfo[EIGENVALUE_VECTOR];
        const std::size_t n_modes = eigenvalues.size();

        // DenseMatrixType eigenvectors;
        const std::size_t n_dofs = this->pGetBuilderAndSolver()->GetEquationSystemSize();
        // if( eigenvectors.size1() != n_modes || eigenvectors.size2() != n_dofs )
        // {
        //     eigenvectors.resize(n_modes, n_dofs, false);
        // }
        
        auto f = this->GetForceVector();

        ComplexType mode_weight;
        ComplexVectorType modal_displacement;
        modal_displacement.resize(n_dofs, false);
        modal_displacement = ZeroVector( n_dofs );

        double modal_damping = 0.0;

        for( std::size_t i = 0; i < n_modes; ++i )
        {
            if( mUseMaterialDamping )
            {
                modal_damping = mMaterialDampingRatios[i];
            }
            else if( rProcessInfo.Has(SYSTEM_DAMPING_RATIO) && rProcessInfo[SYSTEM_DAMPING_RATIO] != 0.0 )
            {
                modal_damping = rProcessInfo[SYSTEM_DAMPING_RATIO];
            }
            else
            {
                modal_damping = mRayleighAlpha / (2 * eigenvalues[i]) + mRayleighBeta * eigenvalues[i] / 2;
            }

            // compute the modal weight
            // DenseVectorType modal_eigenvector = ZeroVector( n_dofs ); //eigenvector for mode i
            // for( ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++ )
            // {
            //     ModelPart::NodeType::DofsContainerType& node_dofs = itNode->GetDofs();
            //     const std::size_t n_node_dofs = node_dofs.size();
            //     Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);

            //     if (node_dofs.IsSorted() == false)
            //     {
            //         node_dofs.Sort();
            //     }

            //     for (std::size_t j = 0; j < n_node_dofs; j++)
            //     {
            //         auto itDof = std::begin(node_dofs) + j;
            //         modal_eigenvector[itDof->EquationId()] = rNodeEigenvectors(i,j);
            //     }
            // }

            auto& pModalMatrix = this->pGetModalMatrix();
            auto& rModalMatrix = *pModalMatrix;

            DenseVectorType modal_vector(n_dofs);
            TDenseSpace::GetColumn(i, rModalMatrix, modal_vector);

            ComplexType factor( eigenvalues[i] - pow( excitation_frequency, 2.0 ), 2 * modal_damping * std::sqrt(eigenvalues[i]) * excitation_frequency );
            mode_weight = inner_prod( modal_vector, f ) / factor;
            // mode_weight = inner_prod( modal_eigenvector, f ) / factor;

            // compute the modal displacement as a superposition of modal_weight * eigenvector
            for( ModelPart::NodeIterator itNode = rModelPart.NodesBegin(); itNode!= rModelPart.NodesEnd(); itNode++ )
            {
                ModelPart::NodeType::DofsContainerType& node_dofs = itNode->GetDofs();
                const std::size_t n_node_dofs = node_dofs.size();
                Matrix& rNodeEigenvectors = itNode->GetValue(EIGENVECTOR_MATRIX);

                if (node_dofs.IsSorted() == false)
                {
                    node_dofs.Sort();
                }

                for (std::size_t j = 0; j < n_node_dofs; j++)
                {
                    auto itDof = std::begin(node_dofs) + j;
                    modal_displacement[itDof->EquationId()] = modal_displacement[itDof->EquationId()] + mode_weight * rNodeEigenvectors(i,j);
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

        BuilderAndSolverPointerType& pBuilderAndSolver = this->pGetBuilderAndSolver();
        SchemePointerType& pScheme = this->pGetScheme();
        auto& pForceVector = this->pGetForceVector();
        auto& rForceVector = *pForceVector;

        // // initialize dummy vectors
        auto pA = SparseSpaceType::CreateEmptyMatrixPointer();
        auto pDx = SparseSpaceType::CreateEmptyVectorPointer();
        // auto pb = SparseSpaceType::CreateEmptyVectorPointer();
        auto& rA = *pA;
        auto& rDx = *pDx;
        // auto& rb = *pb;

        SparseSpaceType::Resize(rA,SparseSpaceType::Size(rForceVector),SparseSpaceType::Size(rForceVector));
        SparseSpaceType::SetToZero(rA);
        SparseSpaceType::Resize(rDx,SparseSpaceType::Size(rForceVector));
        SparseSpaceType::Set(rDx,0.0);
        // // SparseSpaceType::Resize(rb,SparseSpaceType::Size(rForceVector));
        // // SparseSpaceType::Set(rb,0.0);


        // SparseMatrixType* AuxStiffnessMatrix = new SparseMatrixType;
        // mpStiffnessMatrix = boost::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);

        // boost::timer setup_dofs_time;
        // pBuilderAndSolver->SetUpDofSet(pScheme, rModelPart);
        // if (BaseType::GetEchoLevel() > 0 && rank == 0)
        // {
        //     std::cout << "setup_dofs_time : " << setup_dofs_time.elapsed() << std::endl;
        // }

        // // Set global equation ids
        // boost::timer setup_system_time;
        // pBuilderAndSolver->SetUpSystem(rModelPart);
        // if (BaseType::GetEchoLevel() > 0 && rank == 0)
        // {
        //     std::cout << "setup_system_time : " << setup_system_time.elapsed() << std::endl;
        // }

        // // pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, rA, rDx, rb, rModelPart.ElementsArray(), rModelPart.ConditionsArray(), rModelPart.GetProcessInfo());
        // rModelPart.GetProcessInfo()[BUILD_LEVEL] = 2;
        // KRATOS_WATCH(rModelPart.GetProcessInfo())
        // pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, 
        //     mpStiffnessMatrix,
        //     pDx,
        //     pForceVector,
        //     rModelPart.Elements(),
        //     rModelPart.Conditions(),
        //     rModelPart.GetProcessInfo());

        // // void ResizeAndInitializeVectors(
        // //     typename TSchemeType::Pointer pScheme,
        // //     TSystemMatrixPointerType& pA,
        // //     TSystemVectorPointerType& pDx,
        // //     TSystemVectorPointerType& pb,
        // //     ElementsArrayType& rElements,
        // //     ConditionsArrayType& rConditions,
        // //     ProcessInfo& CurrentProcessInfo

        // initial operations ... things that are constant over the solution step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,rForceVector);

        // initial operations ... things that are constant over the solution step
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),rA,rDx,rForceVector);

        // auto test = rModelPart.GetSubModelPartNames();
        // //!!!compute material damping ratio. we need different model parts for the different materials. 
        // //modelparts to be defined in the project parameters
        // //get stiffness matrix only from the model part -> is this then of full size??

        // ///////////// das läuft!!
        // pBuilderAndSolver->BuildLHS(pScheme,rModelPart,*mpStiffnessMatrix);
        // // KRATOS_WATCH(*mpStiffnessMatrix)
        // ////////////////
        // AuxStiffnessMatrix = new SparseMatrixType;
        // mpStiffnessMatrix = boost::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);
        // pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, 
        //     mpStiffnessMatrix,
        //     pDx,
        //     pForceVector,
        //     rModelPart.Elements(),
        //     rModelPart.Conditions(),
        //     rModelPart.GetProcessInfo());

        // for( ModelPart::SubModelPartIterator itSubModelPart = rModelPart.SubModelPartsBegin(); itSubModelPart!= rModelPart.SubModelPartsEnd(); itSubModelPart++ )
        // {
        //     AuxStiffnessMatrix = new SparseMatrixType;
        //     mpStiffnessMatrix = boost::shared_ptr<SparseMatrixType>(AuxStiffnessMatrix);
        //     pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, 
        //         mpStiffnessMatrix,
        //         pDx,
        //         pForceVector,
        //         rModelPart.Elements(),
        //         rModelPart.Conditions(),
        //         rModelPart.GetProcessInfo());
        //     // KRATOS_WATCH(itSubModelPart->Name())
        //     // KRATOS_WATCH(*itSubModelPart)
        //     pBuilderAndSolver->BuildLHS(pScheme,*itSubModelPart,*mpStiffnessMatrix);
        //     // auto index = itSubModelPart - rModelPart.SubModelPartsBegin();
        //     auto index = std::distance(rModelPart.SubModelPartsBegin(), itSubModelPart);
        //     auto damping_ratio = mMaterialDampingRatios[index];
        //     // std::cout << "size=" << rA.size1() << "/" << rA.size2() << " damping=" << damping_ratio << " i=" << index << std::endl;
        //     // KRATOS_WATCH(*mpStiffnessMatrix)
        //     //jetzt die eigenvektoren und das ganze dann in einer neuen membervariable (material damping oder so)
        //     //speichern. das kann dann oben wieder eingesetzt werden
        // }
        
        //no need to set BUILD_LEVEL to sth fancy (check eigensolver_dynamic_scheme.hpp) ??!?!?!?!

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

    DenseMatrixPointerType mpModalMatrix;

    double mRayleighAlpha;

    double mRayleighBeta;

    bool mUseMaterialDamping;

    vector< double > mMaterialDampingRatios;

    // SparseMatrixPointerType mpStiffnessMatrix;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /// Assign the modal displacement to the dofs and the phase angle to the reaction
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
                    itDof->GetSolutionStepReactionValue(step) = std::abs(std::arg(rModalDisplacement(itDof->EquationId())));
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

