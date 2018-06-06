//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Navaneeth K Narayanan
//  Collaborator:    Vicente Mataix Ferrandiz
//

#ifndef KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_
#define KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_

/* System includes */
#include <algorithm>

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "custom_utilities/multipoint_constraint_data.hpp"
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

/**
 * @class ResidualBasedBlockBuilderAndSolverWithMpc
 * @ingroup StructuralMechanicsApplication
 * @brief  Current class provides an implementation for multipoint constraints builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual)
 * Degrees of freedom are reordered putting the restrained degrees of freedom at
 * the end of the system ordered in reverse order with respect to the DofSet.
 * Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
 * this information.
 * Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @author Aditya Ghantasala
 * @author Navaneeth K Narayanan
 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ResidualBasedBlockBuilderAndSolverWithMpc
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ResidualBasedBlockBuilderAndSolverWithMpc
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedBlockBuilderAndSolverWithMpc);

    typedef ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::ElementsContainerType ElementsContainerType;

    /// MPC definitions
    typedef MpcData::Pointer MpcDataPointerType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef MpcData::VariableType VariableType;
    typedef MpcData::MasterDofWeightMapType MasterDofWeightMapType;
    typedef MpcData::MasterIdWeightMapType MasterIdWeightMapType;
    typedef MpcData::SlavePairType SlavePairType;
    typedef MpcData::VariableDataType VariableDataType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef std::vector<MpcDataPointerType> MpcDataPointerVectorType;
    typedef Kratos::shared_ptr<MpcDataPointerVectorType> MpcDataSharedPointerVectorType;

    /// Dof definitions
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;
    typedef Element::EquationIdVectorType EquationIdVectorType;

    // Process info definitions
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;

    /// Index definition
    typedef std::size_t IndexType;
    typedef std::vector<IndexType> VectorIndexType;

    /// Size definition
    typedef std::size_t SizeType;

    /// Node definitions
    typedef Node<3> NodeType;
    typedef ModelPart::NodeIterator NodeIterator;

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ResidualBasedBlockBuilderAndSolverWithMpc(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ResidualBasedBlockBuilderAndSolverWithMpc() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetUpSystem(
        ModelPart &rModelPart
        ) override
    {
        BaseType::SetUpSystem(rModelPart);
        FormulateEquationIdRelationMap(rModelPart);
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        UpdateConstraintEquationsAfterIteration(rModelPart, A, Dx, b);

        Build(pScheme, rModelPart, A, b);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithMpc", ( this->GetEchoLevel() == 3)) << "Before the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        this->SystemSolveWithPhysics(A, Dx, b, rModelPart);

        Timer::Stop("Solve");
        double stop_solve = OpenMPUtils::GetCurrentTime();

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithMpc", (this->GetEchoLevel() >=1 && rModelPart.GetCommunicator().MyPID() == 0)) << "System solve time: " << stop_solve - start_solve << std::endl;

        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithMpc", ( this->GetEchoLevel() == 3)) << "After the solution of the system" << "\nSystem Matrix = " << A << "\nUnknowns vector = " << Dx << "\nRHS vector = " << b << std::endl;

        ReconstructSlaveDofForIterationStep(rModelPart, A, Dx, b); // Reconstructing the slave dofs from master solutions

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    // This is modified to include the MPC information.
    //
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        EquationIdVectorType EquationId;

        // assemble all elements
        double start_build = OpenMPUtils::GetCurrentTime();

#pragma omp parallel firstprivate(nelements, nconditions, LHS_Contribution, RHS_Contribution, EquationId)
        {
#pragma omp for schedule(guided, 512) nowait
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
                    pScheme->CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
                    // Modifying the local contributions for MPC
                    this->Element_ApplyMultipointConstraints(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);
//assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
#else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif
                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }

//#pragma omp parallel for firstprivate(nconditions, LHS_Contribution, RHS_Contribution, EquationId ) schedule(dynamic, 1024)
#pragma omp for schedule(guided, 512)
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
                    pScheme->Condition_CalculateSystemContributions(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

                    // Modifying the local contributions for MPC
                    this->Condition_ApplyMultipointConstraints(*(it.base()), LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

//assemble the elemental contribution
#ifdef USE_LOCKS_IN_ASSEMBLY
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId, BaseType::mlock_array);
#else
                    this->Assemble(A, b, LHS_Contribution, RHS_Contribution, EquationId);
#endif

                    // clean local elemental memory
                    pScheme->CleanMemory(*(it.base()));
                }
            }
        }

        // equation_ids.close();
        // KRATOS_THROW_ERROR(std::logic_error,"i want to stop here :-D","")

        double stop_build = OpenMPUtils::GetCurrentTime();
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithMpc", (this->GetEchoLevel() >= 1 && rModelPart.GetCommunicator().MyPID() == 0)) << "Build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        KRATOS_INFO_IF("ResidualBasedBlockBuilderAndSolverWithMpc", (this->GetEchoLevel() > 2 && rModelPart.GetCommunicator().MyPID() == 0)) << "Finished parallel building" << std::endl;

        KRATOS_CATCH("")
    }
    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType& A,
        ModelPart& rModelPart
        ) override
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

        // Getting the elements from the model
        const int nelements = static_cast<int>(rModelPart.Elements().size());

        // Getting the array of the conditions
        const int nconditions = static_cast<int>(rModelPart.Conditions().size());

        ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = rModelPart.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = rModelPart.ConditionsBegin();

        const std::size_t equation_size = BaseType::mDofSet.size();

#ifdef USE_GOOGLE_HASH
        std::vector<google::dense_hash_set<std::size_t>> indices(equation_size);
        const std::size_t empty_key = 2 * equation_size + 10;
#else
        std::vector<std::unordered_set<std::size_t>> indices(equation_size);
#endif

#pragma omp parallel for firstprivate(equation_size)
        for (int iii = 0; iii < static_cast<int>(equation_size); iii++)
        {
#ifdef USE_GOOGLE_HASH
            indices[iii].set_empty_key(empty_key);
#else
            indices[iii].reserve(40);
#endif
        }
        EquationIdVectorType ids(3, 0);

#pragma omp parallel for firstprivate(nelements, ids)
        for (int iii = 0; iii < nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = el_begin + iii;
            (i_element)->EquationIdVector(ids, CurrentProcessInfo);

            // Modifying the equation IDs of this element to suit MPCs
            this->Element_ModifyEquationIdsForMPC(*(i_element.base()), ids, CurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&(BaseType::mlock_array[ids[i]]));
#endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());

#ifdef _OPENMP
                omp_unset_lock(&(BaseType::mlock_array[ids[i]]));
#endif
            }
        }

#pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii < nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = cond_begin + iii;
            (i_condition)->EquationIdVector(ids, CurrentProcessInfo);
            // Modifying the equation IDs of this element to suit MPCs
            this->Condition_ModifyEquationIdsForMPC(*(i_condition.base()), ids, CurrentProcessInfo);

            for (std::size_t i = 0; i < ids.size(); i++)
            {
#ifdef _OPENMP
                omp_set_lock(&(BaseType::mlock_array[ids[i]]));
#endif
                auto &row_indices = indices[ids[i]];
                row_indices.insert(ids.begin(), ids.end());
#ifdef _OPENMP
                omp_unset_lock(&(BaseType::mlock_array[ids[i]]));
#endif
            }
        }
        //count the row sizes
        unsigned int nnz = 0;
        for (unsigned int i = 0; i < indices.size(); i++)
            nnz += indices[i].size();

        A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);

        double *Avalues = A.value_data().begin();
        std::size_t *Arow_indices = A.index1_data().begin();
        std::size_t *Acol_indices = A.index2_data().begin();

        //filling the index1 vector - DO NOT MAKE PARALLEL THE FOLLOWING LOOP!
        Arow_indices[0] = 0;
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
            Arow_indices[i + 1] = Arow_indices[i] + indices[i].size();

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(A.size1()); i++)
        {
            const unsigned int row_begin = Arow_indices[i];
            const unsigned int row_end = Arow_indices[i + 1];
            unsigned int k = row_begin;
            for (auto it = indices[i].begin(); it != indices[i].end(); it++)
            {
                Acol_indices[k] = *it;
                Avalues[k] = 0.0;
                k++;
            }

            indices[i].clear(); //deallocating the memory

            std::sort(&Acol_indices[row_begin], &Acol_indices[row_end]);
        }

        A.set_filled(indices.size() + 1, nnz);
        Timer::Stop("MatrixStructure");
    }

    /**
     * @brief This function modifies the provided equation ID vector to accommodate MPC constraints
     */
    template<class TContainerType>
    void ModifyEquationIdsForMPC(
        typename TContainerType::Pointer pCurrentContainer,
        EquationIdVectorType &EquationId,
        ProcessInfo &CurrentProcessInfo
        )
    {
        const SizeType number_of_nodes = pCurrentContainer->GetGeometry().PointsNumber();
        MpcDataSharedPointerVectorType p_mpc_data_vector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto p_mpc_data : *p_mpc_data_vector) {
            if (p_mpc_data->IsActive()) {
                // For each node check if it is a slave or not If it is .. we change the Transformation matrix
                for (IndexType j = 0; j < number_of_nodes; j++) {
                    DofsVectorType element_dofs;
                    pCurrentContainer->GetDofList(element_dofs, CurrentProcessInfo);
                    SizeType number_dofs_per_node = element_dofs.size() / number_of_nodes;
                    if (pCurrentContainer->GetGeometry()[j].Is(SLAVE)) { //temporary, will be checked once at the beginning only
                        // Necessary data for iterating and modifying the matrix
                        IndexType slave_equation_id;
                        IndexType start_position_node_dofs = number_dofs_per_node * (j);
                        for (IndexType i = 0; i < number_dofs_per_node; i++) {
                            slave_equation_id = element_dofs[start_position_node_dofs + i]->EquationId();
                            if (p_mpc_data->mEquationIdToWeightsMap.count(slave_equation_id) > 0) {
                                MasterIdWeightMapType master_weights_map = p_mpc_data->mEquationIdToWeightsMap[slave_equation_id];
                                for (auto master : master_weights_map) {
                                    EquationId.push_back(master.first);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief This method is an specialization of ModifyEquationIdsForMPC for elements
     */
    void Element_ModifyEquationIdsForMPC(
        Element::Pointer pCurrentElement,
        EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        ModifyEquationIdsForMPC<Element>(pCurrentElement, rEquationId, rCurrentProcessInfo);
    }

    /**
     * @brief This method is an specialization of ModifyEquationIdsForMPC for conditions
     */
    void Condition_ModifyEquationIdsForMPC(
        Condition::Pointer pCurrentCondition,
        EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        ModifyEquationIdsForMPC<Condition>(pCurrentCondition, rEquationId, rCurrentProcessInfo);
    }

    /**
     * @brief This function changes/extends the element LHS and RHS to apply MPC
     */
    template<class TContainerType>
    void ApplyMultipointConstraints(
        typename TContainerType::Pointer pCurrentContainer,
        LocalSystemMatrixType &LHS_Contribution,
        LocalSystemVectorType &RHS_Contribution,
        EquationIdVectorType &EquationId,
        ProcessInfo &CurrentProcessInfo
        )
    {
        KRATOS_TRY
        bool slaveFound = false;
        auto& geometry = pCurrentContainer->GetGeometry();
        const SizeType number_of_nodes = geometry.PointsNumber();
        for (IndexType j = 0; j < number_of_nodes; j++) {
            if (geometry[j].Is(SLAVE)) { // temporary, will be checked once at the beginning only
                slaveFound = true;
                break;
            }
        }
        // If no slave is found no need of going on
        if (!slaveFound) {
            return;
        }
        MpcDataSharedPointerVectorType p_mpc_data_vector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto p_mpc_data : *p_mpc_data_vector) {
            if (p_mpc_data->IsActive()) {
                VectorIndexType local_equations_ids;
                VectorIndexType local_slave_equation_ids;
                VectorIndexType local_internal_equation_ids;
                VectorIndexType local_master_equations_ids;
                std::vector<double> weights_corresponding_to_master;
                VectorIndexType slaves_corresponding_to_masters;
                // Formulating the local slave equationId vector
                for (IndexType i = 0; i < EquationId.size(); ++i) {
                    local_equations_ids.push_back(i);
                    if (p_mpc_data->mEquationIdToWeightsMap.count(EquationId[i]) > 0) {
                        local_slave_equation_ids.push_back(i);
                    }
                }
                std::sort(local_equations_ids.begin(), local_equations_ids.end());
                std::sort(local_slave_equation_ids.begin(), local_slave_equation_ids.end());
                std::set_difference(local_equations_ids.begin(), local_equations_ids.end(), local_slave_equation_ids.begin(), local_slave_equation_ids.end(), std::back_inserter(local_internal_equation_ids));
                for (IndexType j = 0; j < number_of_nodes; ++j) { // Loop over the nodes
                    VectorIndexType slave_equation_ids;
                    SizeType total_number_of_slaves = 0;
                    SizeType total_number_of_masters = 0;
                    DofsVectorType elementDofs;
                    pCurrentContainer->GetDofList(elementDofs, CurrentProcessInfo);
                    SizeType numDofsPerNode = elementDofs.size() / number_of_nodes;

                    if (geometry[j].Is(SLAVE)) { // If the node has a slave DOF
                        IndexType start_position_node_dofs = numDofsPerNode * (j);
                        IndexType slave_equation_id;
                        for (IndexType i = 0; i < numDofsPerNode; i++) {
                            slave_equation_id = elementDofs[start_position_node_dofs + i]->EquationId();
                            if (p_mpc_data->mEquationIdToWeightsMap.count(slave_equation_id) > 0) {
                                total_number_of_slaves++;
                                slave_equation_ids.push_back(slave_equation_id);
                                MasterIdWeightMapType &master_weights_map = p_mpc_data->mEquationIdToWeightsMap[slave_equation_id];

                                total_number_of_masters += master_weights_map.size();
                            }
                        }

                        VectorIndexType::iterator it;
                        VectorIndexType local_nodel_slave_equations_ids;
                        // We resize the LHS and RHS contribution with the master sizes
                        const SizeType current_sys_size = LHS_Contribution.size1();
                        const SizeType lhs_size_1 = current_sys_size + total_number_of_masters;
                        const SizeType lhs_size_2 = current_sys_size + total_number_of_masters;
                        LHS_Contribution.resize(lhs_size_1, lhs_size_2, true); //true for Preserving the data and resizing the matrix
                        RHS_Contribution.resize(lhs_size_1, true);
                        // Making the extra part of matrx
                        for (IndexType m = current_sys_size; m < lhs_size_1; m++) {
                            for (IndexType n = 0; n < lhs_size_1; n++) {
                                LHS_Contribution(m, n) = 0.0;
                                LHS_Contribution(n, m) = 0.0;
                            }
                            RHS_Contribution(m) = 0.0;
                        }
                        // Formulating the local slave equationId vector
                        for (IndexType slaveI = 0; slaveI < total_number_of_slaves; ++slaveI) { // For each of the Slave DOF
                            // Obtaining the local dof number for the slave.
                            int local_slave_eqid = -1;
                            IndexType slaveEqId = slave_equation_ids[slaveI];
                            it = std::find(EquationId.begin(), EquationId.end(), slaveEqId);
                            if (it != EquationId.end()) {
                                std::size_t pos = std::distance(EquationId.begin(), it);
                                local_slave_eqid = pos;
                            }
                            local_nodel_slave_equations_ids.push_back(local_slave_eqid);
                        }

                        SizeType current_number_of_masters_processed = 0;
                        for (auto local_slave_eqid : local_nodel_slave_equations_ids) { // Loop over all the slaves for this node
                            it = std::find(local_nodel_slave_equations_ids.begin(), local_nodel_slave_equations_ids.end(), local_slave_eqid);
                            const IndexType slave_index = std::distance(local_nodel_slave_equations_ids.begin(), it);
                            MasterIdWeightMapType &master_weights_map = p_mpc_data->mEquationIdToWeightsMap[slave_equation_ids[slave_index]];
                            for (auto masterI : master_weights_map) { // Loop over all the masters the slave has

                                IndexType local_master_eq_id = current_number_of_masters_processed + current_sys_size;
                                ++current_number_of_masters_processed;
                                const double weight = masterI.second;
                                const double constant = p_mpc_data->mSlaveEquationIdConstantsUpdate[slave_equation_ids[slave_index]];
                                for (auto local_intern_eqid : local_internal_equation_ids) {
                                    RHS_Contribution(local_intern_eqid) += -LHS_Contribution(local_intern_eqid, local_slave_eqid) * constant;
                                }

                                // For K(m,u) and K(u,m)
                                for (auto local_intern_eqid : local_internal_equation_ids) { // Loop over all the local equation ids
                                    LHS_Contribution(local_intern_eqid, local_master_eq_id) += LHS_Contribution(local_intern_eqid, local_slave_eqid) * weight;
                                    LHS_Contribution(local_master_eq_id, local_intern_eqid) += LHS_Contribution(local_slave_eqid, local_intern_eqid) * weight;
                                } // Loop over all the local equation ids

                                // For RHS(m) += A'*LHS(s,s)*B
                                for (auto local_slave_eqid_other : local_nodel_slave_equations_ids) {
                                    //VectorIndexType::iterator itOther = std::find(local_nodel_slave_equations_ids.begin(), local_nodel_slave_equations_ids.end(), local_slave_eqid_other);
                                    IndexType slave_index_other = std::distance(local_nodel_slave_equations_ids.begin(), it);
                                    double constantOther = p_mpc_data->mSlaveEquationIdConstantsUpdate[slave_equation_ids[slave_index_other]];
                                    RHS_Contribution(local_master_eq_id) += LHS_Contribution(local_slave_eqid, local_slave_eqid_other) * weight * constantOther;
                                }

                                EquationId.push_back(masterI.first);
                                // Changing the RHS side of the equation
                                RHS_Contribution(local_master_eq_id) += weight * RHS_Contribution(local_slave_eqid);

                                local_master_equations_ids.push_back(local_master_eq_id);
                                weights_corresponding_to_master.push_back(weight);
                                slaves_corresponding_to_masters.push_back(local_slave_eqid);

                            } // Loop over all the masters the slave has

                            RHS_Contribution(local_slave_eqid) = 0.0;
                        } // Loop over all the slaves for this node

                        //Adding contribution from slave to Kmm
                        for (IndexType local_master_index = 0; local_master_index < local_master_equations_ids.size(); local_master_index++) {
                            for (IndexType local_master_index_other = 0; local_master_index_other < local_master_equations_ids.size(); local_master_index_other++) {
                                LHS_Contribution(local_master_equations_ids[local_master_index], local_master_equations_ids[local_master_index_other]) += weights_corresponding_to_master[local_master_index] *
                                                                                                                                             LHS_Contribution(slaves_corresponding_to_masters[local_master_index], slaves_corresponding_to_masters[local_master_index_other]) * weights_corresponding_to_master[local_master_index_other];
                            }
                        }
                    } // If the node has a slave DOF
                }     // Loop over the nodes

                // For K(u,s) and K(s,u)
                for (auto local_slave_eqid : local_slave_equation_ids) { // Loop over all the slaves for this node
                    for (auto local_intern_eqid : local_internal_equation_ids) { // Loop over all the local equation ids
                        LHS_Contribution(local_slave_eqid, local_intern_eqid) = 0.0;
                        LHS_Contribution(local_intern_eqid, local_slave_eqid) = 0.0;
                    }
                } // Loop over all the slaves for this node
            }
        }
        KRATOS_CATCH("Applying Multipoint constraints failed ..");
    } // End of function

    /**
     * @brief This is a specialization of ApplyMultipointConstraints for elements
     */
    void Element_ApplyMultipointConstraints(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        ApplyMultipointConstraints<Element>(pCurrentElement, rLHS_Contribution, rRHS_Contribution, rEquationId, rCurrentProcessInfo);
    }

    /**
     * @brief This is a specialization of ApplyMultipointConstraints for conditions
     */
    void Condition_ApplyMultipointConstraints(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        ApplyMultipointConstraints<Condition>(pCurrentCondition, rLHS_Contribution, rRHS_Contribution, rEquationId, rCurrentProcessInfo);
    }

    /**
     * @brief fThis function Formulates the MPC data in equation ID terms
     * @param rModelPart Reference to the ModelPart containing the problem.
     */
    void FormulateEquationIdRelationMap(ModelPart &rModelPart)
    {
        ProcessInfoType info = rModelPart.GetProcessInfo();

        if (info.Has(MPC_DATA_CONTAINER)) {
            MpcDataSharedPointerVectorType p_p_mpc_data_vector = info.GetValue(MPC_DATA_CONTAINER);
            for (auto p_mpc_data : *p_p_mpc_data_vector) {
                if (p_mpc_data->IsActive()) {
                    for (auto& slave_master_dof_map : p_mpc_data->mDofConstraints) {
                        const SlavePairType slave_dof_map = slave_master_dof_map.first;
                        const MasterDofWeightMapType &master_dof_map = slave_master_dof_map.second;
                        const IndexType slave_node_id = slave_dof_map.first;
                        const IndexType slave_dof_key = slave_dof_map.second;
                        NodeType& slave_node = rModelPart.Nodes()[slave_node_id];
                        NodeType::DofsContainerType::iterator it_slave = slave_node.GetDofs().find(slave_dof_key);
                        if (it_slave != slave_node.GetDofs().end()) {
                            const IndexType slave_equation_id = it_slave->EquationId();

                            for (auto& master_dof_map_elem : master_dof_map) {
                                IndexType master_nodeId;
                                double constant;
                                IndexType master_equation_id;
                                IndexType master_dof_key;
                                const double weight = master_dof_map_elem.second;
                                std::tie(master_nodeId, master_dof_key, constant) = master_dof_map_elem.first;
                                NodeType &master_node = rModelPart.Nodes()[master_nodeId];
                                NodeType::DofsContainerType::iterator it_master = master_node.GetDofs().find(master_dof_key);
                                if (it_master != master_node.GetDofs().end()) {
                                    master_equation_id = it_master->EquationId();

                                    p_mpc_data->AddConstraint(slave_equation_id, master_equation_id, weight, constant);
                                } else {
                                    KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithMpc") << "WARNING:: Could not find dof " << master_dof_key << " in node " << master_node.Id() << std::endl;
                                }
                            }
                        } else {
                            KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithMpc") << "WARNING:: Could not find dof " << slave_dof_key << " in node " << slave_node.Id() << std::endl;
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief This method reconstructs the slave dof for iteration step
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void ReconstructSlaveDofForIterationStep(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b
        )
    {
        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        MpcDataSharedPointerVectorType p_mpc_data_vector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto p_mpc_data : (*p_mpc_data_vector)) {
            if (p_mpc_data->IsActive()) {
                for (auto slave_master_dof_map : p_mpc_data->mDofConstraints) {
                    const SlavePairType slave_dof_map = slave_master_dof_map.first;
                    const MasterDofWeightMapType &master_dof_map = slave_master_dof_map.second;
                    const IndexType slave_node_id = slave_dof_map.first;
                    const IndexType slave_dof_key = slave_dof_map.second;
                    NodeType& slave_node = rModelPart.Nodes()[slave_node_id];
                    NodeType::DofsContainerType::iterator it_slave = slave_node.GetDofs().find(slave_dof_key);
                    if (it_slave != slave_node.GetDofs().end()) {
                        const IndexType slave_equation_id = it_slave->EquationId();
                        double slave_dx_value = 0.0;

                        for (auto master_dof_map_elem : master_dof_map) {
                            IndexType master_nodeId;
                            double constant;
                            IndexType master_dof_key;
                            const double weight = master_dof_map_elem.second;
                            std::tie(master_nodeId, master_dof_key, constant) = master_dof_map_elem.first;
                            NodeType &master_node = rModelPart.Nodes()[master_nodeId];
                            NodeType::DofsContainerType::iterator it_master = master_node.GetDofs().find(master_dof_key);
                            if (it_master != master_node.GetDofs().end()) {
                                const IndexType master_equation_id = it_master->EquationId();
                                slave_dx_value = slave_dx_value + TSparseSpace::GetValue(Dx, master_equation_id) * weight;
                            } else {
                                KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithMpc") << "WARNING:: Could not find dof " << master_dof_key << " in node " << master_node.Id() << std::endl;
                            }
                        }
                        slave_dx_value = slave_dx_value + p_mpc_data->mSlaveEquationIdConstantsUpdate[slave_equation_id];

                        Dx[slave_equation_id] = slave_dx_value;
                        p_mpc_data->mSlaveEquationIdConstantsUpdate[slave_equation_id] = 0.0;
                    } else {
                        KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithMpc") << "WARNING:: Could not find dof " << slave_dof_key << " in node " << slave_node.Id() << std::endl;
                    }
                }
            }
        }
    }

    /**
     * @brief This method is used to update the constraint equation after each iteration
     * @param rModelPart Reference to the ModelPart containing the problem.
     * @param A System matrix
     * @param Dx Vector of results (variations on nodal variables)
     * @param b RHS vector (residual)
     */
    void UpdateConstraintEquationsAfterIteration(
        ModelPart &rModelPart,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b
        )
    {
        ProcessInfo &CurrentProcessInfo = rModelPart.GetProcessInfo();
        MpcDataSharedPointerVectorType p_mpc_data_vector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto p_mpc_data : (*p_mpc_data_vector)){
            if (p_mpc_data->IsActive()) {
                for (auto slave_master_dof_map : p_mpc_data->mDofConstraints) {
                    const SlavePairType slave_dof_map = slave_master_dof_map.first;
                    const MasterDofWeightMapType &master_dof_map = slave_master_dof_map.second;
                    const IndexType slave_node_id = slave_dof_map.first;
                    const IndexType slave_dof_key = slave_dof_map.second;
                    NodeType& slave_node = rModelPart.Nodes()[slave_node_id];
                    NodeType::DofsContainerType::iterator it_slave = slave_node.GetDofs().find(slave_dof_key);
                    if (it_slave != slave_node.GetDofs().end()) {
                        const IndexType slave_equation_id = it_slave->EquationId();
                        const double slave_dof_value = it_slave->GetSolutionStepValue();
                        double slave_dof_value_calc = 0.0;

                        for (auto master_dof_map_elem : master_dof_map) {
                            IndexType master_nodeId;
                            double constant;
                            IndexType master_dof_key;
                            const double weight = master_dof_map_elem.second;
                            std::tie(master_nodeId, master_dof_key, constant) = master_dof_map_elem.first;
                            NodeType& master_node = rModelPart.Nodes()[master_nodeId];
                            NodeType::DofsContainerType::iterator it_master = master_node.GetDofs().find(master_dof_key);
                            if (it_master != master_node.GetDofs().end()) {
                                slave_dof_value_calc += it_master->GetSolutionStepValue() * weight;
                            } else {
                                KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithMpc") << "WARNING:: Could not find dof " << master_dof_key << " in node " << master_node.Id() << std::endl;
                            }
                        }

                        slave_dof_value_calc += p_mpc_data->mSlaveEquationIdConstantsMap[slave_equation_id];

                        const double d_constant = slave_dof_value_calc - slave_dof_value;
                        p_mpc_data->mSlaveEquationIdConstantsUpdate[slave_equation_id] = d_constant;
                    } else {
                        KRATOS_WARNING("ResidualBasedBlockBuilderAndSolverWithMpc") << "WARNING:: Could not find dof " << slave_dof_key << " in node " << slave_node.Id() << std::endl;
                    }
                }
            }
        }
    }
};
}

#endif /* KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_ */
