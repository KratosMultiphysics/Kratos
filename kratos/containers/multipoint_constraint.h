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
//
//

#if !defined(MULTIPOINT_CONSTRAINT_H)
#define MULTIPOINT_CONSTRAINT_H
// System includes
#include <vector>
#include <unordered_map>
#include <iostream>
#include <tuple>
#include <utility>
#include <assert.h>

// project includes
#include <boost/functional/hash.hpp>
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"

namespace Kratos
{
/** \brief MpcData
	* A class that implements the data structure needed for applying Multipoint constraints.
    */
template <class TDenseSpace>
class MultipointConstraint : public Constraint<TDenseSpace>
{

  public:
    /// Pointer definition of MultipointConstraint
    KRATOS_CLASS_POINTER_DEFINITION(MultipointConstraint);

    typedef Dof<double> DofType;
    typedef Constraint<TDenseSpace> BaseType;
    typedef std::unordered_map<unsigned int, double> MasterIdWeightMapType;
    typedef std::pair<unsigned int, unsigned int> SlavePairType;
    typedef Node<3> NodeType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef typename BaseType::MasterDofWeightMapType MasterDofWeightMapType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@name Life Cycle
    ///@{

    /**
		Creates a MultipointConstraint object
		*/
    MultipointConstraint() : BaseType()
    {
    }
    /// Destructor.
    virtual ~MultipointConstraint(){};

    ///@}

    ///@name Access
    ///@{

    virtual void FormulateEquationIdRelationMap(NodesContainerType &Nodes) override
    {
        for (auto slaveMasterDofMap : this->GetData().mDofConstraints)
        {
            SlavePairType slaveDofMap = slaveMasterDofMap.first;
            MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
            unsigned int slaveNodeId = slaveDofMap.first;
            unsigned int slaveDofKey = slaveDofMap.second;
            NodeType &node = Nodes[slaveNodeId];
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
                NodeType &masterNode = Nodes[masterNodeId];
                Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);
                masterEquationId = itMaster->EquationId();
                //
                this->AddConstraint(slaveEquationId, masterEquationId, weight, constant);
            }
        }
    }

    virtual void UpdateConstraintEquationsAfterIteration(NodesContainerType &Nodes) override
    {

        for (auto slaveMasterDofMap : this->GetData().mDofConstraints)
        {
            SlavePairType slaveDofMap = slaveMasterDofMap.first;
            MasterDofWeightMapType &masterDofMap = slaveMasterDofMap.second;
            unsigned int slaveNodeId = slaveDofMap.first;
            unsigned int slaveDofKey = slaveDofMap.second;
            NodeType &node = Nodes[slaveNodeId];
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
                NodeType &masterNode = Nodes[masterNodeId];
                Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);

                slaveDofValueCalc += itMaster->GetSolutionStepValue() * weight;
            }
            slaveDofValueCalc += (this->GetData()).mSlaveEquationIdConstantsMap[slaveEquationId];

            double dConstant = slaveDofValueCalc - slaveDofValue;
            (this->GetData()).mSlaveEquationIdConstantsUpdate[slaveEquationId] = dConstant;
        }
    }

    virtual void Element_ModifyEquationIdsForMPC(Element &rCurrentElement,
                                                 EquationIdVectorType &EquationId,
                                                 ProcessInfo &CurrentProcessInfo)
    {
        const unsigned int number_of_nodes = rCurrentElement.GetGeometry().PointsNumber();
        // For each node check if it is a slave or not If it is .. we change the Transformation matrix
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            DofsVectorType elementDofs;
            rCurrentElement.GetDofList(elementDofs, CurrentProcessInfo);
            int numDofsPerNode = elementDofs.size() / number_of_nodes;
            if (rCurrentElement.GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                // Necessary data for iterating and modifying the matrix
                unsigned int slaveEquationId;
                int startPositionNodeDofs = numDofsPerNode * (j);
                for (int i = 0; i < numDofsPerNode; i++)
                {
                    slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                    if (BaseType::GetData().mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                    {
                        MasterIdWeightMapType &masterWeightsMap = BaseType::GetData().mEquationIdToWeightsMap[slaveEquationId];
                        for (auto &master : masterWeightsMap)
                        {
                            EquationId.push_back(master.first);
                        }
                    }
                }
            }
        }
    }

    virtual void Condition_ModifyEquationIdsForMPC(Condition &rCurrentCondition,
                                                   EquationIdVectorType &EquationId,
                                                   ProcessInfo &CurrentProcessInfo)
    {
        const unsigned int number_of_nodes = rCurrentCondition.GetGeometry().PointsNumber();
        // For each node check if it is a slave or not If it is .. we change the Transformation matrix
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            DofsVectorType elementDofs;
            rCurrentCondition.GetDofList(elementDofs, CurrentProcessInfo);
            int numDofsPerNode = elementDofs.size() / number_of_nodes;
            if (rCurrentCondition.GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                // Necessary data for iterating and modifying the matrix
                unsigned int slaveEquationId;
                int startPositionNodeDofs = numDofsPerNode * (j);
                for (int i = 0; i < numDofsPerNode; i++)
                {
                    slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                    if (BaseType::GetData().mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                    {
                        MasterIdWeightMapType &masterWeightsMap = BaseType::GetData().mEquationIdToWeightsMap[slaveEquationId];
                        for (auto &master : masterWeightsMap)
                        {
                            EquationId.push_back(master.first);
                        }
                    }
                }
            }
        }
    }

    virtual void Element_ApplyConstraints(Element &rCurrentElement,
                                          LocalSystemMatrixType &LHS_Contribution,
                                          LocalSystemVectorType &RHS_Contribution,
                                          EquationIdVectorType &EquationId,
                                          ProcessInfo &CurrentProcessInfo)
    {

        KRATOS_TRY

        bool slaveFound = false;
        Element::NodesArrayType nodesArray = rCurrentElement.GetGeometry();
        const unsigned int number_of_nodes = rCurrentElement.GetGeometry().PointsNumber();
        for (unsigned int j = 0; j < number_of_nodes; j++)
        {
            if (rCurrentElement.GetGeometry()[j].Is(SLAVE))
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
            if (BaseType::GetData().mEquationIdToWeightsMap.count(EquationId[i]) > 0)
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
            rCurrentElement.GetDofList(elementDofs, CurrentProcessInfo);
            int numDofsPerNode = elementDofs.size() / number_of_nodes;

            if (rCurrentElement.GetGeometry()[j].Is(SLAVE))
            { // If the node has a slave DOF
                int startPositionNodeDofs = numDofsPerNode * (j);
                unsigned int slaveEquationId;
                for (int i = 0; i < numDofsPerNode; i++)
                {
                    slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                    if (BaseType::GetData().mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                    {
                        totalNumberOfSlaves++;
                        slaveEquationIds.push_back(slaveEquationId);
                        MasterIdWeightMapType &masterWeightsMap = BaseType::GetData().mEquationIdToWeightsMap[slaveEquationId];

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
                    MasterIdWeightMapType &masterWeightsMap = BaseType::GetData().mEquationIdToWeightsMap[slaveEquationIds[slaveIndex]];
                    for (auto masterI : masterWeightsMap)
                    { // Loop over all the masters the slave has

                        int localMasterEqId = currentNumberOfMastersProcessed + currentSysSize;
                        ++currentNumberOfMastersProcessed;
                        double weight = masterI.second;
                        double constant = BaseType::GetData().mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndex]];
                        for (auto &localInternEqId : localInternEquationIds)
                        {
                            RHS_Contribution(localInternEqId) += -LHS_Contribution(localInternEqId, localSlaveEqId) * constant;
                        }

                        // For K(m,u) and K(u,m)
                        for (auto &localInternEqId : localInternEquationIds)
                        { // Loop over all the local equation ids
                            LHS_Contribution(localInternEqId, localMasterEqId) += LHS_Contribution(localInternEqId, localSlaveEqId) * weight;
                            LHS_Contribution(localMasterEqId, localInternEqId) += LHS_Contribution(localSlaveEqId, localInternEqId) * weight;
                        } // Loop over all the local equation ids

                        // For RHS(m) += A'*LHS(s,s)*B
                        for (auto &localSlaveEqIdOther : localNodalSlaveEquationIds)
                        {
                            //std::vector<std::size_t>::iterator itOther = std::find(localNodalSlaveEquationIds.begin(), localNodalSlaveEquationIds.end(), localSlaveEqIdOther);
                            int slaveIndexOther = std::distance(localNodalSlaveEquationIds.begin(), it);
                            double constantOther = BaseType::GetData().mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndexOther]];
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
        for (auto &localSlaveEqId : localSlaveEquationIds)
        { // Loop over all the slaves for this node
            for (auto &localInternEqId : localInternEquationIds)
            { // Loop over all the local equation ids
                LHS_Contribution(localSlaveEqId, localInternEqId) = 0.0;
                LHS_Contribution(localInternEqId, localSlaveEqId) = 0.0;
            }
        } // Loop over all the slaves for this node

        KRATOS_CATCH("Applying Multipoint constraints failed ..");
    }

    void Condition_ApplyMultipointConstraints(Condition &rCurrentCondition,
                                              LocalSystemMatrixType &LHS_Contribution,
                                              LocalSystemVectorType &RHS_Contribution,
                                              EquationIdVectorType &EquationId,
                                              ProcessInfo &CurrentProcessInfo)
    {
        KRATOS_TRY

        if (CurrentProcessInfo.Has(CONSTRAINTS_CONTAINER))
        {
            bool slaveFound = false;
            Element::NodesArrayType nodesArray = rCurrentCondition.GetGeometry();
            const unsigned int number_of_nodes = rCurrentCondition.GetGeometry().PointsNumber();
            for (unsigned int j = 0; j < number_of_nodes; j++)
            {
                if (rCurrentCondition.GetGeometry()[j].Is(SLAVE))
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
            ConstraintSharedPointerVectorType constraintVector = CurrentProcessInfo.GetValue(CONSTRAINTS_CONTAINER);
            for (auto constraint : (*constraintVector))
            {
                if (constraint->IsActive())
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
                        if (constraint->GetData().mEquationIdToWeightsMap.count(EquationId[i]) > 0)
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
                        rCurrentCondition.GetDofList(elementDofs, CurrentProcessInfo);
                        int numDofsPerNode = elementDofs.size() / number_of_nodes;

                        if (rCurrentCondition.GetGeometry()[j].Is(SLAVE))
                        { // If the node has a slave DOF
                            int startPositionNodeDofs = numDofsPerNode * (j);
                            unsigned int slaveEquationId;
                            for (int i = 0; i < numDofsPerNode; i++)
                            {
                                slaveEquationId = elementDofs[startPositionNodeDofs + i]->EquationId();
                                if (constraint->GetData().mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                                    if (constraint->GetData().mEquationIdToWeightsMap.count(slaveEquationId) > 0)
                                    {
                                        totalNumberOfSlaves++;
                                        slaveEquationIds.push_back(slaveEquationId);
                                        MasterIdWeightMapType &masterWeightsMap = constraint->GetData().mEquationIdToWeightsMap[slaveEquationId];

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
                                MasterIdWeightMapType &masterWeightsMap = constraint->GetData().mEquationIdToWeightsMap[slaveEquationIds[slaveIndex]];
                                for (auto& masterI : masterWeightsMap)
                                { // Loop over all the masters the slave has

                                    int localMasterEqId = currentNumberOfMastersProcessed + currentSysSize;
                                    ++currentNumberOfMastersProcessed;
                                    double weight = masterI.second;
                                    double constant = constraint->GetData().mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndex]];
                                    for (auto& localInternEqId : localInternEquationIds)
                                    {
                                        RHS_Contribution(localInternEqId) += -LHS_Contribution(localInternEqId, localSlaveEqId) * constant;
                                    }

                                    // For K(m,u) and K(u,m)
                                    for (auto& localInternEqId : localInternEquationIds)
                                    { // Loop over all the local equation ids
                                        LHS_Contribution(localInternEqId, localMasterEqId) += LHS_Contribution(localInternEqId, localSlaveEqId) * weight;
                                        LHS_Contribution(localMasterEqId, localInternEqId) += LHS_Contribution(localSlaveEqId, localInternEqId) * weight;
                                    } // Loop over all the local equation ids

                                    // For RHS(m) += A'*LHS(s,s)*B
                                    for (auto& localSlaveEqIdOther : localNodalSlaveEquationIds)
                                    {
                                        //std::vector<std::size_t>::iterator itOther = std::find(localNodalSlaveEquationIds.begin(), localNodalSlaveEquationIds.end(), localSlaveEqIdOther);
                                        int slaveIndexOther = std::distance(localNodalSlaveEquationIds.begin(), it);
                                        double constantOther = constraint->GetData().mSlaveEquationIdConstantsUpdate[slaveEquationIds[slaveIndexOther]];
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
                    for (auto& localSlaveEqId : localSlaveEquationIds)
                    { // Loop over all the slaves for this node
                        for (auto& localInternEqId : localInternEquationIds)
                        { // Loop over all the local equation ids
                            LHS_Contribution(localSlaveEqId, localInternEqId) = 0.0;
                            LHS_Contribution(localInternEqId, localSlaveEqId) = 0.0;
                        }
                    } // Loop over all the slaves for this node
                }
            }
        }
        KRATOS_CATCH("Applying Multipoint constraints failed ..");
        
    }

    /**
		Get the Total number of MasterDOFs for a given slave dof
		@return Total number of MasterDOFs for a given slave dof
		 */
    unsigned int
    GetNumbeOfMasterDofsForSlave(const DofType &SlaveDof)
    {
        return (this->GetData()).mDofConstraints[std::make_pair(SlaveDof.Id(), SlaveDof.GetVariable().Key())].size();
    }

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " MultipointConstraint object " << std::endl;
    }

    ///@name Member Variables
    ///@{

    ///@}

    ///@name Serialization
    ///@{

    ///@}
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif //
