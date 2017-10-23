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

#ifndef APPLY_MULTI_POINT_CONSTRAINTS_PROCESS_H
#define APPLY_MULTI_POINT_CONSTRAINTS_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"

// Application includes
#include "custom_utilities/multipoint_constraint_data.hpp"

namespace Kratos
{

class ApplyMultipointConstraintsProcess : public Process
{
  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyMultipointConstraintsProcess);

    typedef MpcData::Pointer MpcDataPointerType;
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef unsigned int IndexType;
    typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
    typedef MpcData::VariableType VariableType;
    typedef boost::shared_ptr<std::vector<MpcDataPointerType>> MpcDataSharedPointerVectorType;
    typedef ModelPart::NodeIterator NodeIterator;

    /// Constructor.
    ApplyMultipointConstraintsProcess(ModelPart &model_part,
                                      Parameters rParameters) : Process(Flags()), mr_model_part(model_part), m_parameters(rParameters)
    {

        Parameters default_parameters(R"(
            {
                "constraint_set_name":"default",
                "master_sub_model_part_name":"default_master",
                "slave_sub_model_part_name":"default_slave",                
                "variable_names":[""],
                "reform_every_step":false   
            }  )");

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if (info->GetValue(MPC_DATA_CONTAINER) == nullptr)
            info->SetValue(MPC_DATA_CONTAINER, MpcDataSharedPointerVectorType(new std::vector<MpcDataPointerType>()));

        pMpc = MpcDataPointerType(new MpcData());
        std::string name = rParameters["constraint_set_name"].GetString();
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataSharedPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        (*mpcDataVector).push_back(pMpc);

        if (!m_parameters["reform_every_step"].GetBool())
            // Adding the master slave relation between the master and slave sub model parts
            AddMasterSlaveRelation();
    }

    ApplyMultipointConstraintsProcess(ModelPart &model_part, std::string name = "default") : Process(Flags()), mr_model_part(model_part), m_parameters("{}")
    {

        // IMPORTANT : This constructor is not to be used when using this process in the normal KRATOS process_list of python script
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if (info->GetValue(MPC_DATA_CONTAINER) == nullptr)
            info->SetValue(MPC_DATA_CONTAINER, MpcDataSharedPointerVectorType(new std::vector<MpcDataPointerType>()));

        pMpc = MpcDataPointerType(new MpcData());
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataSharedPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        (*mpcDataVector).push_back(pMpc);
    }

    /**
		Applies the MPC condition using two model parts, one as master and other as slave.
        Here a nearest element interpolation is used by default to get the relation between master and slave
		*/
    void AddMasterSlaveRelation()
    {
        ModelPart &master_model_part = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part = mr_model_part.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        int &dim = info->GetValue(DOMAIN_SIZE);

        std::string interpolationType = m_parameters["interpolation_type"].GetString();
        Parameters mapper_parameters = m_parameters["interpolation_settings"];

        if (dim == 2)
        {
            ApplyConstraints<2>(master_model_part, slave_model_part);
        }
        else if (dim == 3)
        {
            ApplyConstraints<3>(master_model_part, slave_model_part);
        }
    }

    // Functions which use two variable components
    template <int TDim>
    void ApplyConstraints(ModelPart &master_model_part, ModelPart &slave_model_part)
    {
        BinBasedFastPointLocator<TDim> *p_point_locator = new BinBasedFastPointLocator<TDim>(master_model_part);
        int numVars = m_parameters["variable_names"].size();
        // iterating over slave nodes to find the corresponding masters
        const int n_slave_nodes = slave_model_part.Nodes().size();
        array_1d<double, TDim + 1> N; // This is only for triangular meshes
        const int max_results = 100;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        for (int i = 0; i < n_slave_nodes; i++)
        {
            ModelPart::NodesContainerType::iterator iparticle = slave_model_part.NodesBegin() + i;
            Node<3>::Pointer p_slave_node = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            Element::Pointer pMasterElement;
            bool is_found = false;
            is_found = p_point_locator->FindPointOnMesh(p_slave_node->Coordinates(), N, pMasterElement, result_begin, max_results);
            if (is_found == true)
            {
                for (int i = 0; i < numVars; i++)
                {
                    std::string varName = m_parameters["variable_names"][i].GetString();
                    Geometry<Node<3>> &geom = pMasterElement->GetGeometry();
                    for (unsigned int i = 0; i < geom.size(); i++)
                    {
                        if (KratosComponents<Variable<double>>::Has(varName))
                        { //case of double variable
                            VariableType rVar = KratosComponents<Variable<double>>::Get(m_parameters["variable_names"][i].GetString());
                            this->AddMasterSlaveRelationWithNodesAndVariable(geom[i], rVar, *p_slave_node, rVar, N[i], 0.0);
                        }
                        else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(varName))
                        {
                            VariableComponentType rVar = KratosComponents<VariableComponentType>::Get(m_parameters["variable_names"][i].GetString());
                            this->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], rVar, *p_slave_node, rVar, N[i], 0.0);
                        }
                    }
                }
            }
        }
        delete p_point_locator;
    }

    /**
		Applies the MPC condition using two nodes, one as master and other as slave, and with the given weight
		@arg MasterNode 
        @arg MasterVariable 
        @arg SlaveNode 
        @arg SlaveVariable
        @arg weight
		*/
    void AddMasterSlaveRelationWithNodesAndVariableComponents(Node<3> &MasterNode, VariableComponentType &MasterVariable, Node<3> &SlaveNode, VariableComponentType &SlaveVariable, double weight, double constant = 0.0)
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariableComponents(IndexType MasterNodeId, VariableComponentType &MasterVariable, IndexType SlaveNodeId, VariableComponentType &SlaveVariable, double weight, double constant = 0.0)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    // Functions with use two variables
    void AddMasterSlaveRelationWithNodesAndVariable(Node<3> &MasterNode, VariableType &MasterVariable, Node<3> &SlaveNode, VariableType &SlaveVariable, double weight, double constant = 0.0)
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariable(IndexType MasterNodeId, VariableType &MasterVariable, IndexType SlaveNodeId, VariableType &SlaveVariable, double weight, double constant = 0.0)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    // Default functions
    /**
		Applies the MPC condition using DOFs, one as master and other as slave, and with the given weight
		@arg slaveDOF 
        @arg masterDOF 
        @arg weight
		*/
    void AddMasterSlaveRelationWithDofs(DofType slaveDOF, DofType masterDOF, double masterWeight, double constant = 0.0)
    {
        pMpc->AddConstraint(slaveDOF, masterDOF, masterWeight, constant);
    }

    /**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/
    void SetActive(bool isActive = true)
    {
        pMpc->SetActive(isActive);
    }

    /**
		Sets the name of the constraint set
		@arg isActive true/false
		*/
    void SetName(std::string name)
    {
        pMpc->SetName(name);
    }

    /// Destructor.
    ~ApplyMultipointConstraintsProcess() override
    {
    }

    void ExecuteBeforeSolutionLoop() override
    {
        KRATOS_TRY;
        //// Use the master and slave sub model parts to formulate the constraints.
        //// Parallel implementation can be taken care here as we can define the partition id from the mapper here.

        KRATOS_CATCH("");
    }

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;

        if (m_parameters["reform_every_step"].GetBool())
            // Adding the master slave relation between the master and slave sub model parts
            AddMasterSlaveRelation();

        KRATOS_CATCH("");
    }

    void ExecuteAfterOutputStep() override
    {
        Clear();
    }

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ApplyMultipointConstraintsProcess";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override { rOStream << "ApplyMultipointConstraintsProcess"; }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
        std::cout << "Number of slave nodes :: " << std::endl;
        pMpc->GetInfo();
    }

    /// Print object's data.
    void Clear()
    {
        pMpc->Clear();
    }

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    ModelPart &mr_model_part;
    MpcDataPointerType pMpc;
    Parameters m_parameters;

  private:
    /// Assignment operator.
    ApplyMultipointConstraintsProcess &operator=(ApplyMultipointConstraintsProcess const &rOther) { return *this; }

}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
