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
//  The modification of the element matrices follows the algorithm described in
//  "AN ALGQRITHM FOR MULTIPOINT CONSTRAINTS IN FINITE ELEMENT ANALYSIS"
//   by John F. Abel and Mark S. Shephard
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
template <class TSparseSpace,
          class TDenseSpace //= DenseSpace<double>
          >
class MultipointConstraint : public Constraint<TSparseSpace, TDenseSpace>
{

  public:
    /// Pointer definition of MultipointConstraint
    KRATOS_CLASS_POINTER_DEFINITION(MultipointConstraint);

    typedef Constraint<TSparseSpace, TDenseSpace> BaseType;
    typedef Node<3> NodeType;
    typedef std::vector<Dof<double>::Pointer> DofsVectorType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

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

    void UpdateConstraintEquations(NodesContainerType &Nodes)
    {
        for (auto &slaveData : this->GetData())
        {
            unsigned int slaveNodeId = slaveData->dofId;
            unsigned int slaveDofKey = slaveData->dofKey;
            NodeType &node = Nodes[slaveNodeId];
            Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slaveDofKey);
            double slaveDofValue = it->GetSolutionStepValue();
            double slaveDofValueCalc = 0.0;
            double constant = slaveData->constant;

            int index = 0;
            for (auto masterDofId : slaveData->masterDofIds)
            {
                unsigned int masterDofKey = slaveData->masterDofKeys[index];
                double weight = slaveData->masterWeights[index];
                NodeType &masterNode = Nodes[masterDofId]; // DofId and nodeId are same
                Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);

                slaveDofValueCalc += itMaster->GetSolutionStepValue() * weight;
            }
            slaveDofValueCalc += constant;

            double dConstant = slaveDofValueCalc - slaveDofValue;
            slaveData->constantUpdate = dConstant;
        }
    }

    virtual void ExecuteBeforeBuilding(NodesContainerType &Nodes) override
    {
        UpdateConstraintEquations(Nodes);
    }


    virtual void SetUp(NodesContainerType &Nodes) override
    {
        for (const auto &slaveData : this->GetData())
        {
            unsigned int slaveNodeId = slaveData->dofId;
            unsigned int slaveDofKey = slaveData->dofKey;
            NodeType &node = Nodes[slaveNodeId];
            Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slaveDofKey);
            unsigned int slaveEquationId = it->EquationId();
            this->GetData().AddEquationIdToSlave(*it, slaveEquationId);

            int index = 0;
            for (auto masterDofId : slaveData->masterDofIds)
            {
                unsigned int masterEquationId;
                unsigned int masterDofKey = slaveData->masterDofKeys[index];

                NodeType &masterNode = Nodes[masterDofId];
                Node<3>::DofsContainerType::iterator itMaster = masterNode.GetDofs().find(masterDofKey);
                masterEquationId = itMaster->EquationId();
                slaveData->SetMasterEquationId(masterDofId, masterDofKey, masterEquationId);
                index++;
            }
        }
        
    }    

    virtual void ExecuteAfterSolving(TSystemMatrixType &A,
                                     TSystemVectorType &Dx,
                                     TSystemVectorType &b) override
    {
        for (auto &slaveData : this->GetData())
        {
            unsigned int slaveEquationId = slaveData->equationId;

            int index = 0;
            for (auto masterEquationId : slaveData->masterEquationIds)
            {
                double weight = slaveData->masterWeights[index];
                Dx[slaveEquationId] = TSparseSpace::GetValue(Dx, slaveEquationId) + TSparseSpace::GetValue(Dx, masterEquationId) * weight;
                index++;
            }

            Dx[slaveEquationId] = TSparseSpace::GetValue(Dx, slaveEquationId) + slaveData->constantUpdate;
            slaveData->constantUpdate = 0.0;
        }
    }

    virtual void Element_ModifyEquationIdsForConstraints(Element &rCurrentElement,
                                                         EquationIdVectorType &EquationId,
                                                         ProcessInfo &CurrentProcessInfo) override
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
                int startPositionNodeDofs = numDofsPerNode * (j);
                for (int i = 0; i < numDofsPerNode; i++)
                {
                    if (this->GetData().GetNumbeOfMasterDofsForSlave(*(elementDofs[startPositionNodeDofs + i])) > 0)
                    {
                        auto &slaveData = this->GetData().GetConstraintEquation(*(elementDofs[startPositionNodeDofs + i]));

                        for (auto masterEqId : slaveData.masterEquationIds)
                        {
                            EquationId.push_back(masterEqId);
                        }
                    }
                }
            }
        }
    }

    virtual void Condition_ModifyEquationIdsForConstraints(Condition &rCurrentCondition,
                                                           EquationIdVectorType &EquationId,
                                                           ProcessInfo &CurrentProcessInfo) override
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
                int startPositionNodeDofs = numDofsPerNode * (j);
                for (int i = 0; i < numDofsPerNode; i++)
                {
                    if (this->GetData().GetNumbeOfMasterDofsForSlave(*(elementDofs[startPositionNodeDofs + i])) > 0)
                    {
                        auto &slaveData = this->GetData().GetConstraintEquation(*(elementDofs[startPositionNodeDofs + i]));

                        for (auto masterEqId : slaveData.masterEquationIds)
                        {
                            EquationId.push_back(masterEqId);
                        }
                    }
                }
            }
        }
    }

    /**
     * (i -> internal, s -> slave, m -> master)
     * We modify the element LHS matrix  K and RHS vector as 
     *             i       s
     *         [ K(i,i)  K(i,s) ]   i
     * K =     [                ]
     *         [ K(s,i)  K(s,s) ]   s
     * 
     * to 
     * 
     *             i                  s              m
     * 
     *         [ K(i,i)             0         K(i,m) + K(i,s)*A     ]      i
     *         [                                                    ]
     * K =     [ 0                  K(s,s)    0                     ]      s 
     *         [                                                    ]
     *         [ K(m,i)+A'*K(s,i)   0         K(m,m) + A*K(s,s)*A'  ]      m
     * 
     * 
     * 
     *       [ RHS(i) ]  i
     * RHS = [        ] 
     *       [ RHS(s) ]  s
     * 
     * to
     * 
     *       [ RHS(i) - K(i,s)*B               ]   i
     *       [                                 ]
     * RHS = [ 0                               ]   s
     *       [                                 ]
     *       [ RHS(m) + A'*RHS(s) + w*K(s,s)*A ]   m
     * 
     */

    virtual void Element_ApplyConstraints(Element &rCurrentElement,
                                          LocalSystemMatrixType &LHS_Contribution,
                                          LocalSystemVectorType &RHS_Contribution,
                                          EquationIdVectorType &EquationId,
                                          ProcessInfo &CurrentProcessInfo) override
    {

        KRATOS_TRY
        ////
        //// Check if there exists a slave in the current element. If not return without any modifications
        //// NOTE : further in the comments indices (written in small) for matrix K and vector RHS : i -> internal, s -> slave, m -> master
        ////
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
        ////
        //// Formulate the list of uneffected  equations for this element. There will be some as all the DOFs cannot be slaves
        ////
        DofsVectorType elementDofs;
        rCurrentElement.GetDofList(elementDofs, CurrentProcessInfo);

        int totalNumberOfMasters = 0;
        int localMasterIndex = 0;
        std::vector<std::size_t> localIndices;
        std::vector<std::size_t> localSlaveIndices;
        std::vector<std::size_t> localInternIndices;

        std::vector<std::size_t> localMasterIndices;
        std::vector<double> localMasterWeights;
        std::vector<std::size_t> slavesIndexToMasterIndex;

        for (unsigned int i = 0; i < elementDofs.size(); ++i)
        {
            localIndices.push_back(i);
            int numMasters = this->GetData().GetNumbeOfMasterDofsForSlave(*elementDofs[i]);

            if (numMasters > 0)
            {
                localSlaveIndices.push_back(i);
                totalNumberOfMasters += numMasters;
            }
        }
        std::sort(localIndices.begin(), localIndices.end());
        std::sort(localSlaveIndices.begin(), localSlaveIndices.end());
        // Get the NOT intersection between localIndices and localSlaveIndices to find the uneffected equations
        std::set_difference(localIndices.begin(), localIndices.end(), localSlaveIndices.begin(), localSlaveIndices.end(), std::back_inserter(localInternIndices));
        ////
        //// Once we know the final size of the modified element system, resizing the element matrix to accomodate the master dofs of the found slave dofs
        ////
        int currentSysSize = LHS_Contribution.size1();
        localMasterIndex = currentSysSize;
        int lhsSize1 = currentSysSize + totalNumberOfMasters;
        int lhsSize2 = currentSysSize + totalNumberOfMasters;
        LHS_Contribution.resize(lhsSize1, lhsSize2, true); //true for Preserving the data and resizing the matrix
        RHS_Contribution.resize(lhsSize1, true);
        // Making the newly added part of matrx, that is K(m,m), K(i,m), k(s,m), k(m,s)
        for (int m = currentSysSize; m < lhsSize1; m++)
        {
            for (int n = 0; n < lhsSize1; n++)
            {
                LHS_Contribution(m, n) = 0.0;
                LHS_Contribution(n, m) = 0.0;
            }
            RHS_Contribution(m) = 0.0;
        }
        ////
        //// Now change the element matrix for each of the slave dof it has for corresponding master dofs
        ////
        for (auto &slaveIndex : localSlaveIndices)
        { // Loop over all the slaves DOFs for this element
            auto &slaveData = this->GetData().GetConstraintEquation(*elementDofs[slaveIndex]);

            int index = 0;
            for (auto masterEqId : slaveData.masterEquationIds)
            { // Loop over all the masters the slave has
                double weight = slaveData.masterWeights[index];
                double constant = slaveData.constant;
                for (auto &internIndex : localInternIndices)
                {
                    RHS_Contribution(internIndex) += -LHS_Contribution(internIndex, slaveIndex) * constant;
                }

                // For K(m,i) and K(i,m)
                for (auto &internIndex : localInternIndices)
                { // Loop over all the local equation ids
                    LHS_Contribution(internIndex, localMasterIndex) += LHS_Contribution(internIndex, slaveIndex) * weight;
                    LHS_Contribution(localMasterIndex, internIndex) += LHS_Contribution(slaveIndex, internIndex) * weight;
                } // Loop over all the local equation ids

                // For RHS(m) += A'*K(s,s)*B
                for (auto &slaveIndexOther : localSlaveIndices)
                {
                    auto slaveDataOther = this->GetData().GetConstraintEquation(*elementDofs[slaveIndexOther]);
                    double constantOther = slaveDataOther.constant;
                    RHS_Contribution(localMasterIndex) += LHS_Contribution(slaveIndex, slaveIndexOther) * weight * constantOther;
                }

                EquationId.push_back(masterEqId);
                // Changing the RHS side of the equation
                RHS_Contribution(localMasterIndex) += weight * RHS_Contribution(slaveIndex);

                localMasterIndices.push_back(localMasterIndex);
                localMasterWeights.push_back(weight);
                slavesIndexToMasterIndex.push_back(slaveIndex);

                index++;
                localMasterIndex++;
            } // Loop over all the masters the slave has

            RHS_Contribution(slaveIndex) = 0.0;
        } // Loop over all the slaves for this node

        // For K(m,m) = K(m,m) + A'*K(s,s)*A
        int index(0);
        for (auto masterIndex : localMasterIndices)
        {
            int indexOther(0);
            for (auto masterIndexOther : localMasterIndices)
            {
                LHS_Contribution(masterIndex, masterIndexOther) += localMasterWeights[index] *
                                                                   LHS_Contribution(slavesIndexToMasterIndex[index], slavesIndexToMasterIndex[indexOther]) *
                                                                   localMasterWeights[indexOther];
                indexOther++;
            }
            index++;
        }

        // For K(i,s) and K(s,i)
        for (auto slaveIndex : localSlaveIndices)
        { // Loop over all the slaves for this node
            for (auto &internIndex : localInternIndices)
            { // Loop over all the local equation ids
                LHS_Contribution(slaveIndex, internIndex) = 0.0;
                LHS_Contribution(internIndex, slaveIndex) = 0.0;
            }
        } // Loop over all the slaves for this node

        KRATOS_CATCH("Applying Multipoint constraints failed ..")
    }

    void Condition_ApplyConstraints(Condition &rCurrentElement,
                                    LocalSystemMatrixType &LHS_Contribution,
                                    LocalSystemVectorType &RHS_Contribution,
                                    EquationIdVectorType &EquationId,
                                    ProcessInfo &CurrentProcessInfo) override
    {
        KRATOS_TRY
        ////
        //// Check if there exists a slave in the current element. If not return without any modifications
        //// NOTE : further in the comments indices (written in small) for matrix K and vector RHS : i -> internal, s -> slave, m -> master
        ////
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
        ////
        //// Formulate the list of uneffected  equations for this element. There will be some as all the DOFs cannot be slaves
        ////
        DofsVectorType elementDofs;
        rCurrentElement.GetDofList(elementDofs, CurrentProcessInfo);

        int totalNumberOfMasters = 0;
        int localMasterIndex = 0;
        std::vector<std::size_t> localIndices;
        std::vector<std::size_t> localSlaveIndices;
        std::vector<std::size_t> localInternIndices;

        std::vector<std::size_t> localMasterIndices;
        std::vector<double> localMasterWeights;
        std::vector<std::size_t> slavesIndexToMasterIndex;

        for (unsigned int i = 0; i < elementDofs.size(); ++i)
        {
            localIndices.push_back(i);
            int numMasters = this->GetData().GetNumbeOfMasterDofsForSlave(*elementDofs[i]);

            if (numMasters > 0)
            {
                localSlaveIndices.push_back(i);
                totalNumberOfMasters += numMasters;
            }
        }
        std::sort(localIndices.begin(), localIndices.end());
        std::sort(localSlaveIndices.begin(), localSlaveIndices.end());
        // Get the NOT intersection between localIndices and localSlaveIndices to find the uneffected equations
        std::set_difference(localIndices.begin(), localIndices.end(), localSlaveIndices.begin(), localSlaveIndices.end(), std::back_inserter(localInternIndices));
        ////
        //// Once we know the final size of the modified element system, resizing the element matrix to accomodate the master dofs of the found slave dofs
        ////
        int currentSysSize = LHS_Contribution.size1();
        localMasterIndex = currentSysSize;
        int lhsSize1 = currentSysSize + totalNumberOfMasters;
        int lhsSize2 = currentSysSize + totalNumberOfMasters;
        LHS_Contribution.resize(lhsSize1, lhsSize2, true); //true for Preserving the data and resizing the matrix
        RHS_Contribution.resize(lhsSize1, true);
        // Making the newly added part of matrx, that is K(m,m), K(i,m), k(s,m), k(m,s)
        for (int m = currentSysSize; m < lhsSize1; m++)
        {
            for (int n = 0; n < lhsSize1; n++)
            {
                LHS_Contribution(m, n) = 0.0;
                LHS_Contribution(n, m) = 0.0;
            }
            RHS_Contribution(m) = 0.0;
        }
        ////
        //// Now change the element matrix for each of the slave dof it has for corresponding master dofs
        ////
        for (auto &slaveIndex : localSlaveIndices)
        { // Loop over all the slaves DOFs for this element
            auto &slaveData = this->GetData().GetConstraintEquation(*elementDofs[slaveIndex]);

            int index = 0;
            for (auto masterEqId : slaveData.masterEquationIds)
            { // Loop over all the masters the slave has
                double weight = slaveData.masterWeights[index];
                double constant = slaveData.constant;
                for (auto &internIndex : localInternIndices)
                {
                    RHS_Contribution(internIndex) += -LHS_Contribution(internIndex, slaveIndex) * constant;
                }

                // For K(m,i) and K(i,m)
                for (auto &internIndex : localInternIndices)
                { // Loop over all the local equation ids
                    LHS_Contribution(internIndex, localMasterIndex) += LHS_Contribution(internIndex, slaveIndex) * weight;
                    LHS_Contribution(localMasterIndex, internIndex) += LHS_Contribution(slaveIndex, internIndex) * weight;
                } // Loop over all the local equation ids

                // For RHS(m) += A'*K(s,s)*B
                for (auto &slaveIndexOther : localSlaveIndices)
                {
                    auto slaveDataOther = this->GetData().GetConstraintEquation(*elementDofs[slaveIndexOther]);
                    double constantOther = slaveDataOther.constant;
                    RHS_Contribution(localMasterIndex) += LHS_Contribution(slaveIndex, slaveIndexOther) * weight * constantOther;
                }

                EquationId.push_back(masterEqId);
                // Changing the RHS side of the equation
                RHS_Contribution(localMasterIndex) += weight * RHS_Contribution(slaveIndex);

                localMasterIndices.push_back(localMasterIndex);
                localMasterWeights.push_back(weight);
                slavesIndexToMasterIndex.push_back(slaveIndex);

                index++;
                localMasterIndex++;
            } // Loop over all the masters the slave has

            RHS_Contribution(slaveIndex) = 0.0;
        } // Loop over all the slaves for this node

        // For K(m,m) = K(m,m) + A'*K(s,s)*A
        int index(0);
        for (auto masterIndex : localMasterIndices)
        {
            int indexOther(0);
            for (auto masterIndexOther : localMasterIndices)
            {
                LHS_Contribution(masterIndex, masterIndexOther) += localMasterWeights[index] *
                                                                   LHS_Contribution(slavesIndexToMasterIndex[index], slavesIndexToMasterIndex[indexOther]) *
                                                                   localMasterWeights[indexOther];
                indexOther++;
            }
            index++;
        }

        // For K(i,s) and K(s,i)
        for (auto slaveIndex : localSlaveIndices)
        { // Loop over all the slaves for this node
            for (auto &internIndex : localInternIndices)
            { // Loop over all the local equation ids
                LHS_Contribution(slaveIndex, internIndex) = 0.0;
                LHS_Contribution(internIndex, slaveIndex) = 0.0;
            }
        } // Loop over all the slaves for this node

        KRATOS_CATCH("Applying Multipoint constraints failed ..")
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

    ///@}*/
};

///@name Input/Output funcitons
///@{

///@}

} // namespace Kratos

#endif //
