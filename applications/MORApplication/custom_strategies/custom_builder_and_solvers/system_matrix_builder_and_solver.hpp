//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//
#if !defined(SYSTEM_MATRIX_BUILDER_AND_SOLVER )
#define  SYSTEM_MATRIX_BUILDER_AND_SOLVER


/* System includes */
#include <unordered_set>
// #include <iostream>
// #include <fstream>

/* External includes */
// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
    #include "sparsehash/dense_hash_set" //included in external libraries
#endif
// #define USE_LOCKS_IN_ASSEMBLY

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
#include "utilities/timer.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_flags.h"

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
 * @class SystemMatrixBuilderAndSolver
 * @ingroup KratosCore
 * @brief Current class provides an implementation for standard builder for mass and stiffness matrices
 * @details This class only builds the mass and stiffness matrices
 * @author Aditya Ghantasala
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class SystemMatrixBuilderAndSolver
    : public ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION(SystemMatrixBuilderAndSolver);


    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef Node<3> NodeType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    SystemMatrixBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >(pNewLinearSystemSolver)
    {

    }

    SystemMatrixBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver, Parameters ThisParameters)
        : ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > (pNewLinearSystemSolver)
    {
        // Validate default parameters
        Parameters default_parameters = Parameters(R"(
        {
            "build_output_structure": false,
            "output_structure_type": "vector",
            "model_part_names" : ["model"],
            "output_variable_names" : ["DISPLACEMENT_Z"]
        })" );

        ThisParameters.ValidateAndAssignDefaults(default_parameters);

        // if( ThisParameters["build_output_vector"].GetBool() && ThisParameters["build_output_matrix"].GetBool() )
        // KRATOS_ERROR_IF(ThisParameters["build_output_vector"].GetBool() && ThisParameters["build_output_matrix"].GetBool() )
        //     << "Cannot build output vector and matrix at the same time!" << std::endl;
        if( ThisParameters["build_output_structure"].GetBool() &&
            ( ThisParameters["output_structure_type"].GetString() == "vector" ) )
                mBuildOutputVectorFlag = true;

        if( ThisParameters["build_output_structure"].GetBool() &&
            ( ThisParameters["output_structure_type"].GetString() == "matrix" ) )
                KRATOS_ERROR << "not yet implemented!\n";

        mOutputModelPartNames = ThisParameters["model_part_names"].GetStringArray();
        mOutputVariableNames = ThisParameters["output_variable_names"].GetStringArray();
    }

    /** Destructor.
     */
    ~SystemMatrixBuilderAndSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void BuildEquationIdVector(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        const int nnodes = static_cast<int>(rModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator n_begin = rModelPart.NodesBegin();

        for (int k=0; k < nnodes; ++k)
        {
            ModelPart::NodesContainerType::iterator it = n_begin + k;

            unsigned int i = 0;
            for (auto const& dof : it->GetDofs())
            {
                b(3*k+i) = (*dof).EquationId();
                i++;
            }
        }
        KRATOS_WATCH(b)

        KRATOS_CATCH("")
    }

    void BuildOutputStructure(
        // typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)//,
        // std::string variableName)
    {
        KRATOS_TRY

        if( mBuildOutputVectorFlag )
        {
            size_t n_model_parts = mOutputModelPartNames.size();
            for( size_t i=0; i<n_model_parts; ++i )
            {
                auto const this_model_part = rModelPart.pGetSubModelPart( mOutputModelPartNames[i] );
                for( auto& node : this_model_part->Nodes() )
                {
                    for( auto& dof : node.GetDofs() )
                    {
                        if( (*dof).GetVariable().Name() == mOutputVariableNames[i] )
                        {
                            b( (*dof).EquationId() ) = 1;
                        }
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void ResizeAndInitializeOutput(
        typename TSchemeType::Pointer pScheme
    )
    {

    }

    /**
     * @brief Function to perform the build of the StiffnessMatrix. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The Stiffnessmatrix matrix
     * @param b The RHS vector
     */
    void BuildStiffnessMatrix(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    //pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                }

            }


            //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    //pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("SystemMatrixBuilderAndSolver - Stiffness matrix", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        KRATOS_INFO_IF("SystemMatrixBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }



    /**
     * @brief Function to perform the build of the StiffnessMatrix. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The Mass matrix
     */
    void BuildMassMatrix(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {
            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {

                ModelPart::ElementsContainerType::iterator it = el_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    //pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    // Directly asking the element to give the stiffnessmatrix
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
                    RHS_Contribution.resize(EquationId.size());

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                }

            }


            //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    //pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateMassMatrix(LHS_Contribution,CurrentProcessInfo);
                    RHS_Contribution.resize(EquationId.size());

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("SystemMatrixBuilderAndSolver - Mass matrix", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        KRATOS_INFO_IF("SystemMatrixBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }

    /*
     * @brief Function to perform the build of the DampingMatrix. The vector could be sized as the total number
     * of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The Damping matrix
     */
    void BuildDampingMatrix(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(!pScheme) << "No scheme provided!" << std::endl;

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());
        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

        #pragma omp parallel firstprivate(nelements,nconditions, LHS_Contribution, RHS_Contribution, EquationId )
        {

            # pragma omp for  schedule(guided, 512) nowait
            for (int k = 0; k < nelements; k++)
            {
                ModelPart::ElementsContainerType::iterator it = el_begin + k;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    element_is_active = (it)->Is(ACTIVE);

                if (element_is_active)
                {
                    //calculate elemental contribution
                    //pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    // Directly asking the element to give the stiffnessmatrix
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateDampingMatrix(LHS_Contribution,CurrentProcessInfo);
                    RHS_Contribution.resize(EquationId.size());


                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);

    #endif
                }
                    //std::cout<<"iteration finished"<<std::endl;
            }


            //#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
            #pragma omp for  schedule(guided, 512)
            for (int k = 0; k < nconditions; k++)
            {
                ModelPart::ConditionsContainerType::iterator it = cond_begin + k;

                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if ((it)->IsDefined(ACTIVE))
                    condition_is_active = (it)->Is(ACTIVE);

                if (condition_is_active)
                {
                    //calculate elemental contribution
                    //pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateDampingMatrix(LHS_Contribution,CurrentProcessInfo);
                    RHS_Contribution.resize(EquationId.size());

                    //assemble the elemental contribution
    #ifdef USE_LOCKS_IN_ASSEMBLY
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, mlock_array);
    #else
                    BaseType::Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
    #endif
                }
            }
        }

        const double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("SystemMatrixBuilderAndSolver - Damping matrix ", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        KRATOS_INFO_IF("SystemMatrixBuilderAndSolver", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }


    /**
     * @brief Function to perform the build of the RHS (not the residual vector).
     * @details The vector could be sized as the total number of dofs or as the number of unrestrained ones
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BuildRHSNoDirichlet(pScheme,rModelPart,b);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for (int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            const std::size_t i = dof_iterator->EquationId();

            if (dof_iterator->IsFixed())
                b[i] = 0.0f;
        }

        KRATOS_CATCH("")

    }

    // void BuildOutputStructure(ModelPart& rModelPart)


    /**
     * @brief Applies the dirichlet conditions specifically for mass matrix to output. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details the rows and colums for Dirichlet conditions are set to zero
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The Mass matrix
     */
    void ApplyDirichletConditionsForMassMatrix(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A)
    {
        std::size_t system_size = A.size1();
        std::vector<double> scaling_factors (system_size, 0.0f);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for(int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            if(dof_iterator->IsFixed())
                scaling_factors[k] = 0.0f;
            else
                scaling_factors[k] = 1.0f;

        }

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        #pragma omp parallel for firstprivate(system_size)
        for(int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            bool empty = true;
            for (std::size_t j = col_begin; j < col_end; ++j)
            {
                if(Avalues[j] != 0.0)
                {
                    empty = false;
                    break;
                }
            }

            if(empty == true)
            {
                A(k,k) = 0.0;
            }
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            double k_factor = scaling_factors[k];
            if (k_factor == 0)
            {
                // zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        }
    }

    /**
     * @brief Applies the dirichlet conditions specifically for mass matrix to output. This operation may be very heavy or completely
     * unexpensive depending on the implementation choosen and on how the System Matrix is built.
     * @details the rows and colums for Dirichlet conditions are set to zero
     * @param pScheme The integration scheme considered
     * @param rModelPart The model part of the problem to solve
     * @param A The Mass matrix
     */
    void ApplyDirichletConditionsForDampingMatrix(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A)
    {
        std::size_t system_size = A.size1();
        std::vector<double> scaling_factors (system_size, 0.0f);

        const int ndofs = static_cast<int>(BaseType::mDofSet.size());

        //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
        #pragma omp parallel for firstprivate(ndofs)
        for(int k = 0; k<ndofs; k++)
        {
            typename DofsArrayType::iterator dof_iterator = BaseType::mDofSet.begin() + k;
            if(dof_iterator->IsFixed())
                scaling_factors[k] = 0.0f;
            else
                scaling_factors[k] = 1.0f;

        }

        double* Avalues = A.value_data().begin();
        std::size_t* Arow_indices = A.index1_data().begin();
        std::size_t* Acol_indices = A.index2_data().begin();

        //detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        #pragma omp parallel for firstprivate(system_size)
        for(int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            bool empty = true;
            for (std::size_t j = col_begin; j < col_end; ++j)
            {
                if(Avalues[j] != 0.0)
                {
                    empty = false;
                    break;
                }
            }

            if(empty == true)
            {
                A(k,k) = 0.0;
            }
        }

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(system_size); ++k)
        {
            std::size_t col_begin = Arow_indices[k];
            std::size_t col_end = Arow_indices[k+1];
            double k_factor = scaling_factors[k];
            if (k_factor == 0)
            {
                // zero out the whole row, except the diagonal
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if (static_cast<int>(Acol_indices[j]) != k )
                        Avalues[j] = 0.0;
            }
            else
            {
                // zero out the column which is associated with the zero'ed row
                for (std::size_t j = col_begin; j < col_end; ++j)
                    if(scaling_factors[ Acol_indices[j] ] == 0 )
                        Avalues[j] = 0.0;
            }
        }
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
        bool mBuildOutputVectorFlag = false;
        std::vector<std::string> mOutputModelPartNames;
        std::vector<std::string> mOutputVariableNames;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void BuildRHSNoDirichlet(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = rModelPart.Conditions();

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        //for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)

        const int nelements = static_cast<int>(pElements.size());
        #pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
        {
            #pragma omp for schedule(guided, 512) nowait
            for(int i=0; i<nelements; i++)
            {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool element_is_active = true;
                if( (it)->IsDefined(ACTIVE) )
                    element_is_active = (it)->Is(ACTIVE);

                if(element_is_active)
                {
                    //calculate elemental Right Hand Side Contribution
                    //pScheme->Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

                    //assemble the elemental contribution
                    BaseType::AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            RHS_Contribution.resize(0, false);

            // assemble all conditions
            //for (typename ConditionsArrayType::ptr_iterator it = ConditionsArray.ptr_begin(); it != ConditionsArray.ptr_end(); ++it)
            const int nconditions = static_cast<int>(ConditionsArray.size());
            //#pragma omp parallel for firstprivate(nconditions, RHS_Contribution, EquationId) schedule(dynamic, 1024)
            #pragma omp for schedule(guided, 512)
            for (int i = 0; i<nconditions; i++)
            {
                auto it = ConditionsArray.begin() + i;
                //detect if the element is active or not. If the user did not make any choice the element
                //is active by default
                bool condition_is_active = true;
                if( (it)->IsDefined(ACTIVE) )
                    condition_is_active = (it)->Is(ACTIVE);

                if(condition_is_active)
                {
                    //calculate conditional contribution
                    //pScheme->Condition_Calculate_RHS_Contribution(*(it.base()), RHS_Contribution, EquationId, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);
                    it->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

                    //assemble the elemental contribution
                    BaseType::AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_CATCH("")

    }
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

    ///@}

}; /* Class SystemMatrixBuilderAndSolver */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */