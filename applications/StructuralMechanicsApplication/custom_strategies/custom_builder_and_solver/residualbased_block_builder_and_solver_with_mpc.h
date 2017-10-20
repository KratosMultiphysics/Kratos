//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala / Navaneeth K Narayanan
//
//

#ifndef KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_
#define KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_

/* System includes */

#include "utilities/openmp_utils.h"
#include <unordered_set>
#include <algorithm>
/* External includes */
#include "boost/smart_ptr.hpp"

#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "includes/model_part.h"
#include "custom_utilities/multipoint_constraint_data.hpp"

// #define USE_GOOGLE_HASH
#ifdef USE_GOOGLE_HASH
#include "sparsehash/dense_hash_set" //included in external libraries
#endif
// #define USE_LOCKS_IN_ASSEMBLY
// #include <iostream>
// #include <fstream>

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */

/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.

Detail class definition.

Current class provides an implementation for standard builder and solving operations.

the RHS is constituted by the unbalanced loads (residual)

Degrees of freedom are reordered putting the restrained degrees of freedom at
the end of the system ordered in reverse order with respect to the DofSet.

Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information.

Calculation of the reactions involves a cost very similiar to the calculation of the total residual

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


 */
template <class TSparseSpace,
          class TDenseSpace,  //= DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ResidualBasedBlockBuilderAndSolverWithMpc
    : public ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    /**@name Type Definitions */
    /*@{ */
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
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;

    typedef Kratos::MpcData::MasterIdWeightMapType MasterIdWeightMapType;
    typedef Kratos::MpcData::SlavePairType SlavePairType;
    typedef Kratos::MpcData::VariableDataType VariableDataType;
    typedef Kratos::MpcData::MasterDofWeightMapType MasterDofWeightMapType;
    typedef Node<3> NodeType;
    typedef MpcData::Pointer MpcDataPointerType;
    //typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
    typedef boost::shared_ptr<std::vector<MpcDataPointerType>> MpcDataPointerVectorType;

    typedef ProcessInfo ProcessInfoType;

    /*@} */
    /**@name Life Cycle
	 */
    /*@{ */

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

    void SetUpSystem(
        ModelPart &r_model_part) override
    {
        BaseType::SetUpSystem(r_model_part);
        FormulateEquationIdRelationMap(r_model_part);
    }

    void BuildAndSolve(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b) override
    {
        KRATOS_TRY

        Timer::Start("Build");

        UpdateConstraintEquationsAfterIteration(r_model_part, A, Dx, b);

        Build(pScheme, r_model_part, A, b);

        Timer::Stop("Build");

        this->ApplyDirichletConditions(pScheme, r_model_part, A, Dx, b);

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "before the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        double start_solve = OpenMPUtils::GetCurrentTime();
        Timer::Start("Solve");

        this->SystemSolveWithPhysics(A, Dx, b, r_model_part);

        Timer::Stop("Solve");
        double stop_solve = OpenMPUtils::GetCurrentTime();
        if (this->GetEchoLevel() >= 1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "system solve time: " << stop_solve - start_solve << std::endl;

        if (this->GetEchoLevel() == 3)
        {
            std::cout << "after the solution of the system" << std::endl;
            std::cout << "System Matrix = " << A << std::endl;
            std::cout << "unknowns vector = " << Dx << std::endl;
            std::cout << "RHS vector = " << b << std::endl;
        }

        ReconstructSlaveDofForIterationStep(r_model_part, A, Dx, b); // Reconstructing the slave dofs from master solutions

        KRATOS_CATCH("")
    }

    //**************************************************************************
    //**************************************************************************
    // This is modified to include the MPC information.
    //
    void Build(
        typename TSchemeType::Pointer pScheme,
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &b) override
    {
        KRATOS_TRY
        if (!pScheme)
            KRATOS_THROW_ERROR(std::runtime_error, "No scheme provided!", "");

        //getting the elements from the model
        const int nelements = static_cast<int>(r_model_part.Elements().size());

        //getting the array of the conditions
        const int nconditions = static_cast<int>(r_model_part.Conditions().size());

        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();
        ModelPart::ElementsContainerType::iterator el_begin = r_model_part.ElementsBegin();
        ModelPart::ConditionsContainerType::iterator cond_begin = r_model_part.ConditionsBegin();

        //contributions to the system
        LocalSystemMatrixType LHS_Contribution = LocalSystemMatrixType(0, 0);
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

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
        if (this->GetEchoLevel() >= 1 && r_model_part.GetCommunicator().MyPID() == 0)
            std::cout << "build time: " << stop_build - start_build << std::endl;

        //for (int i = 0; i < A_size; i++)
        //    omp_destroy_lock(&lock_array[i]);
        if (this->GetEchoLevel() > 2 && r_model_part.GetCommunicator().MyPID() == 0)
        {
            KRATOS_WATCH("finished parallel building");
        }

        KRATOS_CATCH("")
    }

  protected:
    void ConstructMatrixStructure(
        typename TSchemeType::Pointer pScheme,
        TSystemMatrixType &A,
        ElementsContainerType &rElements,
        ConditionsArrayType &rConditions,
        ProcessInfo &CurrentProcessInfo) override
    {
        //filling with zero the matrix (creating the structure)
        Timer::Start("MatrixStructure");

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
        Element::EquationIdVectorType ids(3, 0);

        const int nelements = static_cast<int>(rElements.size());
#pragma omp parallel for firstprivate(nelements, ids)
        for (int iii = 0; iii < nelements; iii++)
        {
            typename ElementsContainerType::iterator i_element = rElements.begin() + iii;
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
        const int nconditions = static_cast<int>(rConditions.size());
#pragma omp parallel for firstprivate(nconditions, ids)
        for (int iii = 0; iii < nconditions; iii++)
        {
            typename ConditionsArrayType::iterator i_condition = rConditions.begin() + iii;
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

  protected:
    /*@} */
    /**@name Protected member Variables */
    /*@{ */

    /*
     * This function modifies the provided equation ID vector to accommodate MPC constraints
     */
    void Element_ModifyEquationIdsForMPC(Element::Pointer rCurrentElement,
                                         Element::EquationIdVectorType &EquationId,
                                         ProcessInfo &CurrentProcessInfo)
    {
        const unsigned int number_of_nodes = rCurrentElement->GetGeometry().PointsNumber();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->IsActive())
            {
                // For each node check if it is a slave or not If it is .. we change the Transformation matrix
                for (unsigned int j = 0; j < number_of_nodes; j++)
                {
                    DofsVectorType elementDofs;
                    rCurrentElement->GetDofList(elementDofs, CurrentProcessInfo);
                    int numDofsPerNode = elementDofs.size() / number_of_nodes;
                    if (rCurrentElement->GetGeometry()[j].Is(SLAVE))
                    { //temporary, will be checked once at the beginning only
                        // Necessary data for iterating and modifying the matrix
                        unsigned int slaveEquationId;
                        int startPositionNodeDofs = numDofsPerNode * (j);
                        for (int i = 0; i < numDofsPerNode; i++)
                        {
                            slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                            if (mpcData->mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                            {
                                MasterIdWeightMapType masterWeightsMap = mpcData->mEquationIdToWeightsMap[slaveEquationId];
                                for (auto master : masterWeightsMap)
                                {
                                    EquationId.push_back(master.first);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     * This function modifies the provided equation ID vector to accommodate MPC constraints
     */
    void Condition_ModifyEquationIdsForMPC(Condition::Pointer rCurrentCondition,
                                           Element::EquationIdVectorType &EquationId,
                                           ProcessInfo &CurrentProcessInfo)
    {

        const unsigned int number_of_nodes = rCurrentCondition->GetGeometry().PointsNumber();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->IsActive())
            {
                // For each node check if it is a slave or not If it is .. we change the Transformation matrix
                for (unsigned int j = 0; j < number_of_nodes; j++)
                {
                    DofsVectorType conditionDofs;
                    rCurrentCondition->GetDofList(conditionDofs, CurrentProcessInfo);
                    int numDofsPerNode = conditionDofs.size() / number_of_nodes;
                    if (rCurrentCondition->GetGeometry()[j].Is(SLAVE))
                    { //temporary, will be checked once at the beginning only
                        // Necessary data for iterating and modifying the matrix
                        unsigned int slaveEquationId;
                        int startPositionNodeDofs = numDofsPerNode * (j);
                        for (int i = 0; i < numDofsPerNode; i++)
                        {
                            slaveEquationId = conditionDofs[startPositionNodeDofs + i]->EquationId();
                            if (mpcData->mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                            {
                                MasterIdWeightMapType masterWeightsMap = mpcData->mEquationIdToWeightsMap[slaveEquationId];
                                for (auto master : masterWeightsMap)
                                {
                                    EquationId.push_back(master.first);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /*
     * This function changes/extends the element LHS and RHS to apply MPC
     */

    void Element_ApplyMultipointConstraints(Element::Pointer rCurrentElement,
                                            LocalSystemMatrixType &LHS_Contribution,
                                            LocalSystemVectorType &RHS_Contribution,
                                            Element::EquationIdVectorType &EquationId,
                                            ProcessInfo &CurrentProcessInfo)
    {

        KRATOS_TRY
        bool slaveFound = false;
        Element::NodesArrayType nodesArray = rCurrentElement->GetGeometry();
        const unsigned int number_of_nodes = rCurrentElement->GetGeometry().PointsNumber();
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            if (rCurrentElement->GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                slaveFound = true;
                break;
            }
        }
        // If no slave is found no need of going on
        if (!slaveFound)
        {
            return;
        }
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->IsActive())
            {
                std::vector<std::size_t> localEquationIds;
                std::vector<std::size_t> localSlaveEquationIds;
                std::vector<std::size_t> localInternEquationIds;
                std::vector<std::size_t> localMasterEquationIds;
                std::vector<double> WeightsCorrespondingToMasters;
                std::vector<std::size_t> SlavesCorrespondingToMasters;
                // Formulating the local slave equationId vector
                for (unsigned int i = 0; i < EquationId.size(); ++i)
                {
                    localEquationIds.push_back(i);
                    if (mpcData->mEquationIdToWeightsMap.count(EquationId[i]) > 0)
                    {
                        localSlaveEquationIds.push_back(i);
                    }
                }
                std::sort(localEquationIds.begin(), localEquationIds.end());
                std::sort(localSlaveEquationIds.begin(), localSlaveEquationIds.end());
                std::set_difference(localEquationIds.begin(), localEquationIds.end(), localSlaveEquationIds.begin(), localSlaveEquationIds.end(), std::back_inserter(localInternEquationIds));
                for (unsigned int j = 0; j < number_of_nodes; ++j)
                { // Loop over the nodes
                    std::vector<int> slaveEquationIds;
                    int totalNumberOfSlaves = 0;
                    int totalNumberOfMasters = 0;
                    DofsVectorType elementDofs;
                    rCurrentElement->GetDofList(elementDofs, CurrentProcessInfo);
                    int numDofsPerNode = elementDofs.size() / number_of_nodes;

                    if (rCurrentElement->GetGeometry()[j].Is(SLAVE))
                    { // If the node has a slave DOF
                        int startPositionNodeDofs = numDofsPerNode * (j);
                        unsigned int slaveEquationId;
                        for (int i = 0; i < numDofsPerNode; i++)
                        {
                            slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                            if (mpcData->mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                            {
                                totalNumberOfSlaves++;
                                slaveEquationIds.push_back(slaveEquationId);
                                MasterIdWeightMapType &masterWeightsMap = mpcData->mEquationIdToWeightsMap[slaveEquationId];

                                totalNumberOfMasters += masterWeightsMap.size();
                            }
                        }

                        std::vector<std::size_t>::iterator it;
                        std::vector<std::size_t> localNodalSlaveEquationIds;
                        // We resize the LHS and RHS contribution with the master sizes
                        int currentSysSize = LHS_Contribution.size1();
                        int lhsSize1 = currentSysSize + totalNumberOfMasters;
                        int lhsSize2 = currentSysSize + totalNumberOfMasters;
                        LHS_Contribution.resize(lhsSize1, lhsSize2, true); //true for Preserving the data and resizing the matrix
                        RHS_Contribution.resize(lhsSize1, true);
                        // Making the extra part of matrx
                        for (int m = currentSysSize; m < lhsSize1; m++)
                        {
                            for (int n = 0; n < lhsSize1; n++)
                            {
                                LHS_Contribution(m, n) = 0.0;
                                LHS_Contribution(n, m) = 0.0;
                            }
                            RHS_Contribution(m) = 0.0;
                        }
                        // Formulating the local slave equationId vector
                        for (int slaveI = 0; slaveI < totalNumberOfSlaves; ++slaveI)
                        { // For each of the Slave DOF
                            // Obtaining the local dof number for the slave.
                            int localSlaveEqId = -1;
                            int slaveEqId = slaveEquationIds[slaveI];
                            it = std::find(EquationId.begin(), EquationId.end(), slaveEqId);
                            if (it != EquationId.end())
                            {
                                std::size_t pos = std::distance(EquationId.begin(), it);
                                localSlaveEqId = pos;
                            }
                            localNodalSlaveEquationIds.push_back(localSlaveEqId);
                        }

                        int currentNumberOfMastersProcessed = 0;
                        for (auto localSlaveEqId : localNodalSlaveEquationIds)
                        { // Loop over all the slaves for this node
                            it = std::find(localNodalSlaveEquationIds.begin(), localNodalSlaveEquationIds.end(), localSlaveEqId);
                            int slaveIndex = std::distance(localNodalSlaveEquationIds.begin(), it);
                            MasterIdWeightMapType &masterWeightsMap = mpcData->mEquationIdToWeightsMap[slaveEquationIds[slaveIndex]];
                            for (auto masterI : masterWeightsMap)
                            { // Loop over all the masters the slave has

                                int localMasterEqId = currentNumberOfMastersProcessed + currentSysSize;
                                ++currentNumberOfMastersProcessed;
                                double weight = masterI.second;
                                double constant = mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndex]];
                                for (auto localInternEqId : localInternEquationIds)
                                {
                                    RHS_Contribution(localInternEqId) += -LHS_Contribution(localInternEqId, localSlaveEqId) * constant;
                                }

                                // For K(m,u) and K(u,m)
                                for (auto localInternEqId : localInternEquationIds)
                                { // Loop over all the local equation ids
                                    LHS_Contribution(localInternEqId, localMasterEqId) += LHS_Contribution(localInternEqId, localSlaveEqId) * weight;
                                    LHS_Contribution(localMasterEqId, localInternEqId) += LHS_Contribution(localSlaveEqId, localInternEqId) * weight;
                                } // Loop over all the local equation ids

                                // For RHS(m) += A'*LHS(s,s)*B
                                for (auto localSlaveEqIdOther : localNodalSlaveEquationIds)
                                {
                                    //std::vector<std::size_t>::iterator itOther = std::find(localNodalSlaveEquationIds.begin(), localNodalSlaveEquationIds.end(), localSlaveEqIdOther);
                                    int slaveIndexOther = std::distance(localNodalSlaveEquationIds.begin(), it);
                                    double constantOther = mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndexOther]];
                                    RHS_Contribution(localMasterEqId) += LHS_Contribution(localSlaveEqId, localSlaveEqIdOther) * weight * constantOther;
                                }

                                EquationId.push_back(masterI.first);
                                // Changing the RHS side of the equation
                                RHS_Contribution(localMasterEqId) += weight * RHS_Contribution(localSlaveEqId);

                                localMasterEquationIds.push_back(localMasterEqId);
                                WeightsCorrespondingToMasters.push_back(weight);
                                SlavesCorrespondingToMasters.push_back(localSlaveEqId);

                            } // Loop over all the masters the slave has

                            RHS_Contribution(localSlaveEqId) = 0.0;
                        } // Loop over all the slaves for this node

                        //Adding contribution from slave to Kmm
                        for (unsigned int localMasterIndex = 0; localMasterIndex < localMasterEquationIds.size(); localMasterIndex++)
                        {
                            for (unsigned int localMasterIndexOther = 0; localMasterIndexOther < localMasterEquationIds.size(); localMasterIndexOther++)
                            {
                                LHS_Contribution(localMasterEquationIds[localMasterIndex], localMasterEquationIds[localMasterIndexOther]) += WeightsCorrespondingToMasters[localMasterIndex] *
                                                                                                                                             LHS_Contribution(SlavesCorrespondingToMasters[localMasterIndex], SlavesCorrespondingToMasters[localMasterIndexOther]) * WeightsCorrespondingToMasters[localMasterIndexOther];
                            }
                        }
                    } // If the node has a slave DOF
                }     // Loop over the nodes

                // For K(u,s) and K(s,u)
                for (auto localSlaveEqId : localSlaveEquationIds)
                { // Loop over all the slaves for this node
                    for (auto localInternEqId : localInternEquationIds)
                    { // Loop over all the local equation ids
                        LHS_Contribution(localSlaveEqId, localInternEqId) = 0.0;
                        LHS_Contribution(localInternEqId, localSlaveEqId) = 0.0;
                    }
                } // Loop over all the slaves for this node
            }
        }
        KRATOS_CATCH("Applying Multipoint constraints failed ..");
    } // End of function

    void Condition_ApplyMultipointConstraints(Condition::Pointer rCurrentElement,
                                              LocalSystemMatrixType &LHS_Contribution,
                                              LocalSystemVectorType &RHS_Contribution,
                                              Element::EquationIdVectorType &EquationId,
                                              ProcessInfo &CurrentProcessInfo)
    {

        KRATOS_TRY
        bool slaveFound = false;
        Element::NodesArrayType nodesArray = rCurrentElement->GetGeometry();
        const unsigned int number_of_nodes = rCurrentElement->GetGeometry().PointsNumber();
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            if (rCurrentElement->GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                slaveFound = true;
                break;
            }
        }
        // If no slave is found no need of going on
        if (!slaveFound)
        {
            return;
        }
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);
        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->IsActive())
            {
                std::vector<std::size_t> localEquationIds;
                std::vector<std::size_t> localSlaveEquationIds;
                std::vector<std::size_t> localInternEquationIds;
                std::vector<std::size_t> localMasterEquationIds;
                std::vector<double> WeightsCorrespondingToMasters;
                std::vector<std::size_t> SlavesCorrespondingToMasters;
                // Formulating the local slave equationId vector
                for (unsigned int i = 0; i < EquationId.size(); ++i)
                {
                    localEquationIds.push_back(i);
                    if (mpcData->mEquationIdToWeightsMap.count(EquationId[i]) > 0)
                    {
                        localSlaveEquationIds.push_back(i);
                    }
                }
                std::sort(localEquationIds.begin(), localEquationIds.end());
                std::sort(localSlaveEquationIds.begin(), localSlaveEquationIds.end());
                std::set_difference(localEquationIds.begin(), localEquationIds.end(), localSlaveEquationIds.begin(), localSlaveEquationIds.end(), std::back_inserter(localInternEquationIds));
                for (unsigned int j = 0; j < number_of_nodes; ++j)
                { // Loop over the nodes
                    std::vector<int> slaveEquationIds;
                    int totalNumberOfSlaves = 0;
                    int totalNumberOfMasters = 0;
                    DofsVectorType elementDofs;
                    rCurrentElement->GetDofList(elementDofs, CurrentProcessInfo);
                    int numDofsPerNode = elementDofs.size() / number_of_nodes;

                    if (rCurrentElement->GetGeometry()[j].Is(SLAVE))
                    { // If the node has a slave DOF
                        int startPositionNodeDofs = numDofsPerNode * (j);
                        unsigned int slaveEquationId;
                        for (int i = 0; i < numDofsPerNode; i++)
                        {
                            slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                            if (mpcData->mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                                if (mpcData->mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                                {
                                    totalNumberOfSlaves++;
                                    slaveEquationIds.push_back(slaveEquationId);
                                    MasterIdWeightMapType &masterWeightsMap = mpcData->mEquationIdToWeightsMap[slaveEquationId];

                                    totalNumberOfMasters += masterWeightsMap.size();
                                }
                        }

                        std::vector<std::size_t>::iterator it;
                        std::vector<std::size_t> localNodalSlaveEquationIds;
                        // We resize the LHS and RHS contribution with the master sizes
                        int currentSysSize = LHS_Contribution.size1();
                        int lhsSize1 = currentSysSize + totalNumberOfMasters;
                        int lhsSize2 = currentSysSize + totalNumberOfMasters;
                        LHS_Contribution.resize(lhsSize1, lhsSize2, true); //true for Preserving the data and resizing the matrix
                        RHS_Contribution.resize(lhsSize1, true);
                        // Making the extra part of matrx
                        for (int m = currentSysSize; m < lhsSize1; m++)
                        {
                            for (int n = 0; n < lhsSize1; n++)
                            {
                                LHS_Contribution(m, n) = 0.0;
                                LHS_Contribution(n, m) = 0.0;
                            }
                            RHS_Contribution(m) = 0.0;
                        }
                        // Formulating the local slave equationId vector
                        for (int slaveI = 0; slaveI < totalNumberOfSlaves; ++slaveI)
                        { // For each of the Slave DOF
                            // Obtaining the local dof number for the slave.
                            int localSlaveEqId = -1;
                            int slaveEqId = slaveEquationIds[slaveI];
                            it = std::find(EquationId.begin(), EquationId.end(), slaveEqId);
                            if (it != EquationId.end())
                            {
                                std::size_t pos = std::distance(EquationId.begin(), it);
                                localSlaveEqId = pos;
                            }
                            localNodalSlaveEquationIds.push_back(localSlaveEqId);
                        }

                        int currentNumberOfMastersProcessed = 0;
                        for (auto localSlaveEqId : localNodalSlaveEquationIds)
                        { // Loop over all the slaves for this node
                            it = std::find(localNodalSlaveEquationIds.begin(), localNodalSlaveEquationIds.end(), localSlaveEqId);
                            int slaveIndex = std::distance(localNodalSlaveEquationIds.begin(), it);
                            MasterIdWeightMapType &masterWeightsMap = mpcData->mEquationIdToWeightsMap[slaveEquationIds[slaveIndex]];
                            for (auto masterI : masterWeightsMap)
                            { // Loop over all the masters the slave has

                                int localMasterEqId = currentNumberOfMastersProcessed + currentSysSize;
                                ++currentNumberOfMastersProcessed;
                                double weight = masterI.second;
                                double constant = mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndex]];
                                for (auto localInternEqId : localInternEquationIds)
                                {
                                    RHS_Contribution(localInternEqId) += -LHS_Contribution(localInternEqId, localSlaveEqId) * constant;
                                }

                                // For K(m,u) and K(u,m)
                                for (auto localInternEqId : localInternEquationIds)
                                { // Loop over all the local equation ids
                                    LHS_Contribution(localInternEqId, localMasterEqId) += LHS_Contribution(localInternEqId, localSlaveEqId) * weight;
                                    LHS_Contribution(localMasterEqId, localInternEqId) += LHS_Contribution(localSlaveEqId, localInternEqId) * weight;
                                } // Loop over all the local equation ids

                                // For RHS(m) += A'*LHS(s,s)*B
                                for (auto localSlaveEqIdOther : localNodalSlaveEquationIds)
                                {
                                    //std::vector<std::size_t>::iterator itOther = std::find(localNodalSlaveEquationIds.begin(), localNodalSlaveEquationIds.end(), localSlaveEqIdOther);
                                    int slaveIndexOther = std::distance(localNodalSlaveEquationIds.begin(), it);
                                    double constantOther = mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndexOther]];
                                    RHS_Contribution(localMasterEqId) += LHS_Contribution(localSlaveEqId, localSlaveEqIdOther) * weight * constantOther;
                                }

                                EquationId.push_back(masterI.first);
                                // Changing the RHS side of the equation
                                RHS_Contribution(localMasterEqId) = RHS_Contribution(localMasterEqId) + weight * RHS_Contribution(localSlaveEqId);

                                localMasterEquationIds.push_back(localMasterEqId);
                                WeightsCorrespondingToMasters.push_back(weight);
                                SlavesCorrespondingToMasters.push_back(localSlaveEqId);

                            } // Loop over all the masters the slave has

                            RHS_Contribution(localSlaveEqId) = 0.0;
                        } // Loop over all the slaves for this node

                        //Adding contribution from slave to Kmm
                        for (unsigned int localMasterIndex = 0; localMasterIndex < localMasterEquationIds.size(); localMasterIndex++)
                        {
                            for (unsigned int localMasterIndexOther = 0; localMasterIndexOther < localMasterEquationIds.size(); localMasterIndexOther++)
                            {
                                LHS_Contribution(localMasterEquationIds[localMasterIndex], localMasterEquationIds[localMasterIndexOther]) += WeightsCorrespondingToMasters[localMasterIndex] *
                                                                                                                                             LHS_Contribution(SlavesCorrespondingToMasters[localMasterIndex], SlavesCorrespondingToMasters[localMasterIndexOther]) * WeightsCorrespondingToMasters[localMasterIndexOther];
                            }
                        }
                    } // If the node has a slave DOF
                }     // Loop over the nodes

                // For K(u,s) and K(s,u)
                for (auto localSlaveEqId : localSlaveEquationIds)
                { // Loop over all the slaves for this node
                    for (auto localInternEqId : localInternEquationIds)
                    { // Loop over all the local equation ids
                        LHS_Contribution(localSlaveEqId, localInternEqId) = 0.0;
                        LHS_Contribution(localInternEqId, localSlaveEqId) = 0.0;
                    }
                } // Loop over all the slaves for this node
            }
        }
        KRATOS_CATCH("Applying Multipoint constraints failed ..");
    } // End of the function

    /*
     * This function Formulates the MPC data in equation ID terms
     */
    void FormulateEquationIdRelationMap(ModelPart &r_model_part)
    {

        ProcessInfoType info = r_model_part.GetProcessInfo();

        if (info.Has(MPC_DATA_CONTAINER))
        {
            MpcDataPointerVectorType mpcDataVector = info.GetValue(MPC_DATA_CONTAINER);
            for (auto mpcData : (*mpcDataVector))
            {
                if (mpcData->IsActive())
                {
                    for (auto slaveMasterDofMap : mpcData->mDofConstraints)
                    {
                        SlavePairType slaveDofMap = slaveMasterDofMap.first;
                        MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                        unsigned int slaveNodeId = slaveDofMap.first;
                        unsigned int slaveDofKey = slaveDofMap.second;
                        NodeType &node = r_model_part.Nodes()[slaveNodeId];
                        Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slaveDofKey);
                        unsigned int slaveEquationId = it->EquationId();

                        for (auto masterDofMapElem : masterDofMap)
                        {
                            unsigned int masterNodeId;
                            double constant;
                            unsigned int masterEquationId;
                            unsigned int masterDofKey;
                            double weight = masterDofMapElem.second;
                            std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                            NodeType &masterNode = r_model_part.Nodes()[masterNodeId];
                            Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);
                            masterEquationId = itMaster->EquationId();
                            //
                            mpcData->AddConstraint(slaveEquationId, masterEquationId, weight, constant);
                        }
                    }
                }
            }
        }
    }

    void ReconstructSlaveDofForIterationStep(
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b)
    {
        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->IsActive())
            {
                for (auto slaveMasterDofMap : mpcData->mDofConstraints)
                {
                    SlavePairType slaveDofMap = slaveMasterDofMap.first;
                    MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                    unsigned int slaveNodeId = slaveDofMap.first;
                    unsigned int slaveDofKey = slaveDofMap.second;
                    NodeType &node = r_model_part.Nodes()[slaveNodeId];
                    Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slaveDofKey);
                    unsigned int slaveEquationId = it->EquationId();

                    for (auto masterDofMapElem : masterDofMap)
                    {
                        unsigned int masterNodeId;
                        double constant;
                        unsigned int masterEquationId;
                        unsigned int masterDofKey;
                        double weight = masterDofMapElem.second;
                        std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                        NodeType &masterNode = r_model_part.Nodes()[masterNodeId];
                        Node<3>::DofsContainerType::iterator it = masterNode.GetDofs().find(masterDofKey);
                        masterEquationId = it->EquationId();

                        Dx[slaveEquationId] = TSparseSpace::GetValue(Dx, slaveEquationId) + TSparseSpace::GetValue(Dx, masterEquationId) * weight;
                    }

                    Dx[slaveEquationId] = TSparseSpace::GetValue(Dx, slaveEquationId) + mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationId];
                    mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationId] = 0.0;
                }
            }
        }
    }

    void UpdateConstraintEquationsAfterIteration(
        ModelPart &r_model_part,
        TSystemMatrixType &A,
        TSystemVectorType &Dx,
        TSystemVectorType &b)
    {

        ProcessInfo &CurrentProcessInfo = r_model_part.GetProcessInfo();
        MpcDataPointerVectorType mpcDataVector = CurrentProcessInfo.GetValue(MPC_DATA_CONTAINER);

        for (auto mpcData : (*mpcDataVector))
        {
            if (mpcData->IsActive())
            {
                for (auto slaveMasterDofMap : mpcData->mDofConstraints)
                {
                    SlavePairType slaveDofMap = slaveMasterDofMap.first;
                    MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
                    unsigned int slaveNodeId = slaveDofMap.first;
                    unsigned int slaveDofKey = slaveDofMap.second;
                    NodeType &node = r_model_part.Nodes()[slaveNodeId];
                    Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slaveDofKey);
                    unsigned int slaveEquationId = it->EquationId();
                    double slaveDofValue = it->GetSolutionStepValue();
                    double slaveDofValueCalc = 0.0;

                    for (auto masterDofMapElem : masterDofMap)
                    {
                        unsigned int masterNodeId;
                        double constant;
                        unsigned int masterDofKey;
                        double weight = masterDofMapElem.second;
                        std::tie(masterNodeId, masterDofKey, constant) = masterDofMapElem.first;
                        NodeType &masterNode = r_model_part.Nodes()[masterNodeId];
                        Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);

                        slaveDofValueCalc += itMaster->GetSolutionStepValue() * weight;
                    }

                    slaveDofValueCalc += mpcData->mSlaveEquationIdConstantsMap[slaveEquationId];

                    double dConstant = slaveDofValueCalc - slaveDofValue;
                    mpcData->mSlaveEquationIdConstantsUpdate[slaveEquationId] = dConstant;
                }
            }
        }
    }
};
}

#endif /* KRATOS_SOLVING_STRATEGIES_BUILDER_AND_SOLVERS_RESIDUALBASED_BLOCK_BUILDER_AND_SOLVER_WITH_MPC_H_ */
