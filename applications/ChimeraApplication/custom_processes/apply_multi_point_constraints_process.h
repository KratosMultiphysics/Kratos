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

// Application includes
#include "custom_utilities/multipoint_constraint_data.hpp"
#include "chimera_application_variables.h"

namespace Kratos
{

class ApplyMultipointConstraintsProcess : public Process
{
  public:
    /// Pointer definition of MoveRotorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyMultipointConstraintsProcess);

    typedef MpcData::Pointer MpcDataPointerType;
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef std::size_t IndexType;
    typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
    typedef MpcData::VariableType VariableType;
    typedef ModelPart::NodeIterator NodeIterator;

    /// Constructor.
    ApplyMultipointConstraintsProcess(std::string type, ModelPart &model_part, ModelPart &rSubModelPart,
                                      Parameters rParameters) : Process(Flags()), mtype(type), mr_model_part(model_part), mrSubModelPart(rSubModelPart), m_parameters(rParameters)
    {

        Parameters default_parameters(R"(
            {
                "constraint_set_name":"default",
                "master_sub_model_part_name":"default_master",
                "slave_sub_model_part_name":"default_slave",
                "variable_names":[""],
                "interpolation_type":"nearest_node",
                "reform_every_step":false
            }  )");

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if (info->GetValue(MPC_DATA_CONTAINER) == NULL)
            info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

        std::string SubModelPartName = mrSubModelPart.Name();
        pMpc = MpcDataPointerType(new MpcData(type, SubModelPartName));
        std::string name = rParameters["constraint_set_name"].GetString();
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        (*mpcDataVector).push_back(pMpc);

        std::string interpolationType = rParameters["interpolation_type"].GetString();
        if (interpolationType != "nearest_element" && interpolationType != "nearest_node")
        {
            KRATOS_THROW_ERROR(std::runtime_error, "No valid interpolation type provided !", "");
        }

        //AddMasterSlaveRelation();
    }

    ApplyMultipointConstraintsProcess(std::string type, ModelPart &model_part, ModelPart &rSubModelPart, std::string name = "default") : Process(Flags()), mtype(type), mr_model_part(model_part), mrSubModelPart(rSubModelPart), m_parameters("{}")
    {
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if (info->GetValue(MPC_DATA_CONTAINER) == NULL)
            info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

        std::string SubModelPartName = mrSubModelPart.Name();
        pMpc = MpcDataPointerType(new MpcData(type, SubModelPartName));
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        (*mpcDataVector).push_back(pMpc);
    }

    // Functions which use two variable components

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

    void AddMasterSlaveRelationWithNodeIdsAndVariable(IndexType MasterNodeId, VariableType &MasterVariable, IndexType SlaveNodeId, VariableType &SlaveVariable, double weight,double constant = 0.0)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }


    // Remove constraints
    void RemoveMasterSlaveRelationWithNodesAndVariableComponents(Node<3> &MasterNode, VariableComponentType &MasterVariable, Node<3> &SlaveNode, VariableComponentType &SlaveVariable)
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        RemoveMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF);
    }

    void RemoveMasterSlaveRelationWithNodeIdsAndVariableComponents(IndexType MasterNodeId, VariableComponentType &MasterVariable, IndexType SlaveNodeId, VariableComponentType &SlaveVariable)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        RemoveMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF);
    }

    // Functions with use two variables
    void RemoveMasterSlaveRelationWithNodesAndVariable(Node<3> &MasterNode, VariableType &MasterVariable, Node<3> &SlaveNode, VariableType &SlaveVariable)
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        RemoveMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF);
    }

    void RemoveMasterSlaveRelationWithNodeIdsAndVariable(IndexType MasterNodeId, VariableType &MasterVariable, IndexType SlaveNodeId, VariableType &SlaveVariable)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        RemoveMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF);
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

    void RemoveMasterSlaveRelationWithDofs(DofType slaveDOF, DofType masterDOF)
    {
        pMpc->RemoveConstraint(slaveDOF);
    }

    void AddNodalNormalSlaveRelationWithDofs(DofType slaveDOF, double nodalNormalComponent = 0.0)
    {
        pMpc->AddNodalNormalToSlaveDof(slaveDOF, nodalNormalComponent);
    }

    /**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/
    void SetActive(bool isActive = true)
    {
        pMpc->SetActive(isActive);
    }

    void SetRtMinvR(double value)
    {

        pMpc->RtMinvR = value;
    }

    /**
		Sets the name of the constraint set
		@arg isActive true/false
		*/
    void SetName(std::string name)
    {
        pMpc->SetName(name);
    }

    void SetWeak(bool value = true)
    {

        pMpc->SetIsWeak(value);
    }

    /// Destructor.
    virtual ~ApplyMultipointConstraintsProcess()
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

        /*if (m_parameters["reform_every_step"].GetBool())
            // Adding the master slave relation between the master and slave sub model parts
            AddMasterSlaveRelation();*/

        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ApplyMultipointConstraintsProcess";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override { rOStream << "ApplyMultipointConstraintsProcess"; }

    /// Print object's data.
    void PrintData()
    {
        pMpc->GetInfo();
    }

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    std::string mtype;
    ModelPart &mr_model_part;
    ModelPart &mrSubModelPart;
    MpcDataPointerType pMpc;
    Parameters m_parameters;

  private:
    /// Assignment operator.
    ApplyMultipointConstraintsProcess &operator=(ApplyMultipointConstraintsProcess const &rOther) { return *this; }


}; // Class MoveRotorProcess

} // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
