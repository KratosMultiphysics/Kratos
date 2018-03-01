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
//  The modification of the element matrix follows the algorithm described in
//  "AN ALGORITHM FOR MULTIPOINT CONSTRAINTS IN FINITE ELEMENT ANALYSIS"
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
    typedef Dof<double> DofType;    
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

    void UpdateConstraintEquations(NodesContainerType &rNodes)
    {
        for (auto &constraint_eq_data : this->GetData())
        {
            size_t slave_node_id = constraint_eq_data->SlaveDofId();
            size_t slave_dof_key = constraint_eq_data->SlaveDofKey();
            NodeType &node = rNodes[slave_node_id];
            Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slave_dof_key);
            double slave_dof_value = it->GetSolutionStepValue();
            double slave_dof_value_calc = 0.0;
            double constant = constraint_eq_data->Constant();

            for (auto &master_data : *constraint_eq_data)
            {
                size_t master_dof_key = master_data->MasterDofKey();
                double weight = master_data->MasterWeight();
                NodeType &r_master_node = rNodes[master_data->MasterDofId()]; // DofId and nodeId are same
                Node<3>::DofsContainerType::iterator it_master = r_master_node.GetDofs().find(master_dof_key);

                slave_dof_value_calc += it_master->GetSolutionStepValue() * weight;
            }
            slave_dof_value_calc += constant;

            double d_constant = slave_dof_value_calc - slave_dof_value;
            constraint_eq_data->SetConstantUpdate(d_constant);
        }
    }

    virtual void ExecuteBeforeBuilding(NodesContainerType &rNodes) override
    {
        UpdateConstraintEquations(rNodes);
    }

    virtual void SetUp(NodesContainerType &rNodes) override
    {
        for (const auto &constraint_eq_data : this->GetData())
        {
            size_t slave_node_id = constraint_eq_data->SlaveDofId();
            size_t slave_dof_key = constraint_eq_data->SlaveDofKey();
            NodeType &node = rNodes[slave_node_id];
            Node<3>::DofsContainerType::iterator it = node.GetDofs().find(slave_dof_key);
            constraint_eq_data->SetSlaveEquationId(it->EquationId());
            for (auto &master_data : *constraint_eq_data)
            {
                size_t master_dof_key = master_data->MasterDofKey();
                size_t master_dof_id = master_data->MasterDofId();

                NodeType &master_node = rNodes[master_dof_id];
                Node<3>::DofsContainerType::iterator it_master = master_node.GetDofs().find(master_dof_key);
                master_data->SetMasterEqId(it_master->EquationId());
            }
        } 
    }

    virtual void ExecuteAfterSolving(TSystemMatrixType &rA,
                                     TSystemVectorType &rDx,
                                     TSystemVectorType &rb) override
    {
        for (auto &constraint_eq_data : this->GetData())
        {
            size_t slave_equation_id = constraint_eq_data->SlaveEquationId();
            for (auto &master_data : *constraint_eq_data)
            {
                rDx[slave_equation_id] = TSparseSpace::GetValue(rDx, slave_equation_id) + TSparseSpace::GetValue(rDx, master_data->MasterEqId()) * master_data->MasterWeight();
            }

            rDx[slave_equation_id] = TSparseSpace::GetValue(rDx, slave_equation_id) + constraint_eq_data->ConstantUpdate();
            constraint_eq_data->SetConstantUpdate(0.0);
        }        
    }

    virtual void ModifyEquationIdsForConstraints(Element &rCurrentElement,
                                                         EquationIdVectorType &rEquationId,
                                                         ProcessInfo &rCurrentProcessInfo) override
    {
        const size_t number_of_nodes = rCurrentElement.GetGeometry().PointsNumber();
        // For each node check if it is a slave or not If it is .. we change the Transformation matrix
        for (size_t j = 0; j < number_of_nodes; j++)
        {            
            DofsVectorType element_dofs;
            rCurrentElement.GetDofList(element_dofs, rCurrentProcessInfo);
            int num_dofs_per_node = element_dofs.size() / number_of_nodes;
            if (rCurrentElement.GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                // Necessary data for iterating and modifying the matrix
                int start_position_node_dofs = num_dofs_per_node * (j);
                for (int i = 0; i < num_dofs_per_node; i++)
                {
                    if (this->GetData().GetNumberOfMasterDofsForSlave(*(element_dofs[start_position_node_dofs + i])) > 0)
                    {
                        auto &constraint_eq_data = this->GetData().GetConstraintEquation(*(element_dofs[start_position_node_dofs + i]));
                        for (auto &master_data : constraint_eq_data)
                        {
                            rEquationId.push_back(master_data->MasterEqId());
                        }
                    }
                }
            }
        }    
    }

    virtual void ModifyEquationIdsForConstraints(Condition &rCurrentCondition,
                                                           EquationIdVectorType &rEquationId,
                                                           ProcessInfo &rCurrentProcessInfo) override
    {
        const size_t number_of_nodes = rCurrentCondition.GetGeometry().PointsNumber();
        // For each node check if it is a slave or not If it is .. we change the Transformation matrix
        for (size_t j = 0; j < number_of_nodes; j++)
        {
            DofsVectorType element_dofs;
            rCurrentCondition.GetDofList(element_dofs, rCurrentProcessInfo);
            int num_dofs_per_node = element_dofs.size() / number_of_nodes;
            if (rCurrentCondition.GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                // Necessary data for iterating and modifying the matrix
                int start_position_node_dofs = num_dofs_per_node * (j);
                for (int i = 0; i < num_dofs_per_node; i++)
                {
                    if (this->GetData().GetNumberOfMasterDofsForSlave(*(element_dofs[start_position_node_dofs + i])) > 0)
                    {
                        auto &constraint_eq_data = this->GetData().GetConstraintEquation(*(element_dofs[start_position_node_dofs + i]));

                        for (auto &master_data : constraint_eq_data)
                        {
                            rEquationId.push_back(master_data->MasterEqId());
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

    virtual void ApplyConstraints(Element &rCurrentElement,
                                          LocalSystemMatrixType &rLHS_Contribution,
                                          LocalSystemVectorType &rRHS_Contribution,
                                          EquationIdVectorType &rEquationId,
                                          ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY
        ////
        //// Check if there exists a slave in the current element. If not return without any modifications
        //// NOTE : further in the comments indices (written in small) for matrix K and vector RHS : i -> internal, s -> slave, m -> master
        ////
        bool slave_found = false;
        const size_t number_of_nodes = rCurrentElement.GetGeometry().PointsNumber();
        for (size_t j = 0; j < number_of_nodes; j++)
        {
            if (rCurrentElement.GetGeometry()[j].Is(SLAVE))
            { //temporary, will be checked once at the beginning only
                slave_found = true;
                break;
            }
        }
        // If no slave is found no need of going on
        if (!slave_found)
        {
            return;
        }
        ////
        //// Formulate the list of uneffected  equations for this element. There will be some as all the DOFs cannot be slaves
        ////
        DofsVectorType element_dofs;
        rCurrentElement.GetDofList(element_dofs, rCurrentProcessInfo);

        int total_number_of_masters = 0;
        int local_master_index = 0;
        std::vector<std::size_t> local_indices;
        std::vector<std::size_t> local_slave_indices;
        std::vector<std::size_t> local_intern_indices;

        std::vector<std::size_t> local_master_indices;
        std::vector<double> local_master_weights;
        std::vector<std::size_t> slaves_index_to_master_index;

        for (size_t i = 0; i < element_dofs.size(); ++i)
        {
            local_indices.push_back(i);
            int num_masters = this->GetData().GetNumberOfMasterDofsForSlave(*element_dofs[i]);

            if (num_masters > 0)
            {
                local_slave_indices.push_back(i);
                total_number_of_masters += num_masters;
            }
        }
        std::sort(local_indices.begin(), local_indices.end());
        std::sort(local_slave_indices.begin(), local_slave_indices.end());
        // Get the NOT intersection between local_indices and local_slave_indices to find the uneffected equations
        std::set_difference(local_indices.begin(), local_indices.end(), local_slave_indices.begin(), local_slave_indices.end(), std::back_inserter(local_intern_indices));
        ////
        //// Once we know the final size of the modified element system, resizing the element matrix to accomodate the master dofs of the found slave dofs
        ////
        int current_system_size = rLHS_Contribution.size1();
        local_master_index = current_system_size;
        int lhs_size1 = current_system_size + total_number_of_masters;
        int lhs_size2 = current_system_size + total_number_of_masters;
        rLHS_Contribution.resize(lhs_size1, lhs_size2, true); //true for Preserving the data and resizing the matrix
        rRHS_Contribution.resize(lhs_size1, true);
        // Making the newly added part of matrx, that is K(m,m), K(i,m), k(s,m), k(m,s)
        for (int m = current_system_size; m < lhs_size1; m++)
        {
            for (int n = 0; n < lhs_size1; n++)
            {
                rLHS_Contribution(m, n) = 0.0;
                rLHS_Contribution(n, m) = 0.0;
            }
            rRHS_Contribution(m) = 0.0;
        }
        ////
        //// Now change the element matrix for each of the slave dof it has for corresponding master dofs
        ////
        for (auto &slave_index : local_slave_indices)
        { // Loop over all the slaves DOFs for this element
            auto &constraint_eq_data = this->GetData().GetConstraintEquation(*element_dofs[slave_index]);

            for (auto &master_data : constraint_eq_data)
            { // Loop over all the masters the slave has
                double weight = master_data->MasterWeight();
                const double constant_update = constraint_eq_data.ConstantUpdate();
                for (auto &intern_index : local_intern_indices)
                {
                    rRHS_Contribution(intern_index) += -rLHS_Contribution(intern_index, slave_index) * constant_update;
                }

                // For K(m,i) and K(i,m)
                for (auto &intern_index : local_intern_indices)
                { // Loop over all the local equation ids
                    rLHS_Contribution(intern_index, local_master_index) += rLHS_Contribution(intern_index, slave_index) * weight;
                    rLHS_Contribution(local_master_index, intern_index) += rLHS_Contribution(slave_index, intern_index) * weight;
                } // Loop over all the local equation ids

                // For RHS(m) += A'*K(s,s)*B
                for (auto &slave_index_other : local_slave_indices)
                {
                    auto& slave_data_other = this->GetData().GetConstraintEquation(*element_dofs[slave_index_other]);
                    double constant_update_other = slave_data_other.ConstantUpdate();
                    rRHS_Contribution(local_master_index) += rLHS_Contribution(slave_index, slave_index_other) * weight * constant_update_other;
                }

                rEquationId.push_back(master_data->MasterEqId());
                // Changing the RHS side of the equation
                rRHS_Contribution(local_master_index) += weight * rRHS_Contribution(slave_index);

                local_master_indices.push_back(local_master_index);
                local_master_weights.push_back(weight);
                slaves_index_to_master_index.push_back(slave_index);

                local_master_index++;
            } // Loop over all the masters the slave has

            rRHS_Contribution(slave_index) = 0.0;
        } // Loop over all the slaves for this node

        // For K(m,m) = K(m,m) + A'*K(s,s)*A
        int index(0);
        for (auto masterIndex : local_master_indices)
        {
            int index_other(0);
            for (auto masterIndexOther : local_master_indices)
            {
                rLHS_Contribution(masterIndex, masterIndexOther) += local_master_weights[index] *
                                                                   rLHS_Contribution(slaves_index_to_master_index[index], slaves_index_to_master_index[index_other]) *
                                                                   local_master_weights[index_other];
                index_other++;
            }
            index++;
        }

        // For K(i,s) and K(s,i)
        for (auto& slave_index : local_slave_indices)
        { // Loop over all the slaves for this node
            for (auto &intern_index : local_intern_indices)
            { // Loop over all the local equation ids
                rLHS_Contribution(slave_index, intern_index) = 0.0;
                rLHS_Contribution(intern_index, slave_index) = 0.0;
            }
        } // Loop over all the slaves for this node


        // For K(s,s) = I
        for (auto& slave_index : local_slave_indices)
        { // Loop over all the slaves for this node
            for (auto& slave_index_other : local_slave_indices)
            { // Loop over all the local equation ids
                rLHS_Contribution(slave_index, slave_index_other) = 0.0;
            }
            rLHS_Contribution(slave_index, slave_index) = 1.0;
        } // Loop over all the slaves for this node        

        KRATOS_CATCH("Applying Multipoint constraints failed ..")
    }

    void ApplyConstraints(Condition &rCurrentElement,
                                    LocalSystemMatrixType &rLHS_Contribution,
                                    LocalSystemVectorType &rRHS_Contribution,
                                    EquationIdVectorType &rEquationId,
                                    ProcessInfo &rCurrentProcessInfo) override
    {
        KRATOS_TRY
        ////
        //// Check if there exists a slave in the current element. If not return without any modifications
        //// NOTE : further in the comments indices (written in small) for matrix K and vector RHS : i -> internal, s -> slave, m -> master
        ////
        bool slaveFound = false;
        const size_t number_of_nodes = rCurrentElement.GetGeometry().PointsNumber();
        for (size_t j = 0; j < number_of_nodes; j++)
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
        DofsVectorType element_dofs;
        rCurrentElement.GetDofList(element_dofs, rCurrentProcessInfo);

        int total_number_of_masters = 0;
        int local_master_index = 0;
        std::vector<std::size_t> local_indices;
        std::vector<std::size_t> local_slave_indices;
        std::vector<std::size_t> local_intern_indices;

        std::vector<std::size_t> local_master_indices;
        std::vector<double> local_master_weights;
        std::vector<std::size_t> slaves_index_to_master_index;

        for (size_t i = 0; i < element_dofs.size(); ++i)
        {
            local_indices.push_back(i);
            int num_masters = this->GetData().GetNumberOfMasterDofsForSlave(*element_dofs[i]);

            if (num_masters > 0)
            {
                local_slave_indices.push_back(i);
                total_number_of_masters += num_masters;
            }
        }
        std::sort(local_indices.begin(), local_indices.end());
        std::sort(local_slave_indices.begin(), local_slave_indices.end());
        // Get the NOT intersection between local_indices and local_slave_indices to find the uneffected equations
        std::set_difference(local_indices.begin(), local_indices.end(), local_slave_indices.begin(), local_slave_indices.end(), std::back_inserter(local_intern_indices));
        ////
        //// Once we know the final size of the modified element system, resizing the element matrix to accomodate the master dofs of the found slave dofs
        ////
        int current_system_size = rLHS_Contribution.size1();
        local_master_index = current_system_size;
        int lhsSize1 = current_system_size + total_number_of_masters;
        int lhsSize2 = current_system_size + total_number_of_masters;
        rLHS_Contribution.resize(lhsSize1, lhsSize2, true); //true for Preserving the data and resizing the matrix
        rRHS_Contribution.resize(lhsSize1, true);
        // Making the newly added part of matrx, that is K(m,m), K(i,m), k(s,m), k(m,s)
        for (int m = current_system_size; m < lhsSize1; m++)
        {
            for (int n = 0; n < lhsSize1; n++)
            {
                rLHS_Contribution(m, n) = 0.0;
                rLHS_Contribution(n, m) = 0.0;
            }
            rRHS_Contribution(m) = 0.0;
        }
        ////
        //// Now change the element matrix for each of the slave dof it has for corresponding master dofs
        ////
        for (auto &slave_index : local_slave_indices)
        { // Loop over all the slaves DOFs for this element
            auto &constraint_eq_data = this->GetData().GetConstraintEquation(*element_dofs[slave_index]);

            for (auto &master_data : constraint_eq_data)
            { // Loop over all the masters the slave has
                double weight = master_data->MasterWeight();
                const double constant_update = constraint_eq_data.ConstantUpdate();
                for (auto &intern_index : local_intern_indices)
                {
                    rRHS_Contribution(intern_index) += -rLHS_Contribution(intern_index, slave_index) * constant_update;
                }

                // For K(m,i) and K(i,m)
                for (auto &intern_index : local_intern_indices)
                { // Loop over all the local equation ids
                    rLHS_Contribution(intern_index, local_master_index) += rLHS_Contribution(intern_index, slave_index) * weight;
                    rLHS_Contribution(local_master_index, intern_index) += rLHS_Contribution(slave_index, intern_index) * weight;
                } // Loop over all the local equation ids

                // For RHS(m) += A'*K(s,s)*B
                for (auto &slave_index_other : local_slave_indices)
                {
                    auto slave_data_other = this->GetData().GetConstraintEquation(*element_dofs[slave_index_other]);
                    double constant_update_other = slave_data_other.ConstantUpdate();
                    rRHS_Contribution(local_master_index) += rLHS_Contribution(slave_index, slave_index_other) * weight * constant_update_other;
                }

                rEquationId.push_back(master_data->MasterEqId());
                // Changing the RHS side of the equation
                rRHS_Contribution(local_master_index) += weight * rRHS_Contribution(slave_index);

                local_master_indices.push_back(local_master_index);
                local_master_weights.push_back(weight);
                slaves_index_to_master_index.push_back(slave_index);

                local_master_index++;
            } // Loop over all the masters the slave has

            rRHS_Contribution(slave_index) = 0.0;
        } // Loop over all the slaves for this node

        // For K(m,m) = K(m,m) + A'*K(s,s)*A
        int index(0);
        for (auto masterIndex : local_master_indices)
        {
            int index_other(0);
            for (auto masterIndexOther : local_master_indices)
            {
                rLHS_Contribution(masterIndex, masterIndexOther) += local_master_weights[index] *
                                                                   rLHS_Contribution(slaves_index_to_master_index[index], slaves_index_to_master_index[index_other]) *
                                                                   local_master_weights[index_other];
                index_other++;
            }
            index++;
        }

        // For K(i,s) and K(s,i)
        for (auto slave_index : local_slave_indices)
        { // Loop over all the slaves for this node
            for (auto &intern_index : local_intern_indices)
            { // Loop over all the local equation ids
                rLHS_Contribution(slave_index, intern_index) = 0.0;
                rLHS_Contribution(intern_index, slave_index) = 0.0;
            }
        } // Loop over all the slaves for this node

        // For K(s,s) = I
        for (auto& slave_index : local_slave_indices)
        { // Loop over all the slaves for this node
            for (auto& slave_index_other : local_slave_indices)
            { // Loop over all the local equation ids
                rLHS_Contribution(slave_index, slave_index_other) = 0.0;
            }
            rLHS_Contribution(slave_index, slave_index) = 1.0;
        } // Loop over all the slaves for this node        

        KRATOS_CATCH("Applying Multipoint constraints failed ..")
    }

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << " MultipointConstraint object " << std::endl;
    }


    /**
	* Get the constraint equation for this slave
	* @return MasterDOFs vector for this slave
	*/
    const ConstraintEquation &GetConstraintEquation (DofType &rSlaveDof)
    {
        return this->GetData().GetConstraintEquation(rSlaveDof);
    } 

    /**
    * Adds a constraints between the given slave and master with a weight. 		
	*/
    // Takes in a slave dof equationId and a master dof equationId
    void AddConstraint(DofType &rSlaveDof, DofType &rMasterDof, double Weight, double Constant = 0.0)
    {
        this->GetData().AddConstraint(rSlaveDof, rMasterDof, Weight, Constant);
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
