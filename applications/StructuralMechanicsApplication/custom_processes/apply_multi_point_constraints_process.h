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

namespace Kratos
{

class ApplyMultipointConstraintsProcess : public Process
{
public:

    /// Pointer definition of MoveRotorProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyMultipointConstraintsProcess);

    typedef MpcData::Pointer MpcDataPointerType;
    typedef Dof<double>* DofPointerType;
    typedef Dof<double> DofType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef ProcessInfo      ProcessInfoType;
    typedef ProcessInfo::Pointer      ProcessInfoPointerType;
    typedef unsigned int IndexType;
    typedef std::vector<MpcDataPointerType>*  MpcDataPointerVectorType;
    typedef MpcData::VariableType VariableType;
    typedef ModePart::NodeIterator NodeIterator;

    /// Constructor.
    ApplyMultipointConstraintsProcess(  ModelPart& model_part,
                                        Parameters rParameters
                                        ) : Process(Flags()) , mr_model_part(model_part)
    {

         Parameters default_parameters( R"(
            {
                "constraint_set_name":"default",
                "master_sub_model_part_name":"default_master",
                "slave_sub_model_part_name":"default_slave",                
                "variable_names":[],
                "interpolation_type":"nearest_node",
                "reform_every_step":false,   
            }  )" );

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if(info->GetValue(MPC_DATA_CONTAINER) == NULL)
            info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

        pMpc = MpcDataPointerType( new MpcData() );
        std::string name = rParameters["constraint_set_name"].GetString();
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        (*mpcDataVector).push_back(pMpc);

        std::string interpolationType = rParameters["interpolation_type"].GetString();
        if(interpolationType != "nearest_element" && interpolationType != "nearest_node" )
            {
                KRATOS_THROW_ERROR(std::runtime_error, "No valid interpolation type provided !", "");
            }
    }

    ApplyMultipointConstraintsProcess(  ModelPart& model_part, std::string name="default"
                                        ) : Process(Flags()) , mr_model_part(model_part)
    {
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if(info->GetValue(MPC_DATA_CONTAINER) == NULL)
            info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

        pMpc = MpcDataPointerType( new MpcData() );
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        (*mpcDataVector).push_back(pMpc);

        // Adding the master slave relation between the master and slave sub model parts
        AddMasterSlaveRelation();
    }


    /**
		Applies the MPC condition using two model parts, one as master and other as slave.
        Here a nearest element interpolation is used by default to get the relation between master and slave
		*/    
    void AddMasterSlaveRelation()
    {
        ModelPart& master_model_part = mr_model_part.GetSubModelPart(rParameters["master_sub_model_part_name"].GetString());
        ModelPart& slave_model_part  = mr_model_part.GetSubModelPart(rParameters["slave_sub_model_part_name"].GetString());
        std::string interpolationType = rParameters["interpolation_type"].GetString();

        // Create the mapper based on the type of interpolation 


        // Use the mapper to get the weights and masters for any given slave
        for(NodeIterator itSlaveNode = slave_model_part.NodesBegin(); itSlaveNode<slave_model_part.NodesEnd(); itSlaveNode++)
        {
            // Find the masters and their weights here and add them to the data using AddMasterSlaveRelationWithDofs function
        }
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
    void AddMasterSlaveRelationWithNodesAndVariableComponents(Node<3> &MasterNode, VariableComponentType& MasterVariable, Node<3> &SlaveNode, VariableComponentType& SlaveVariable, double weight)
    {
        SlaveNode.Set(SLAVE);        
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
    	DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }


    void AddMasterSlaveRelationWithNodeIdsAndVariableComponents(IndexType MasterNodeId, VariableComponentType& MasterVariable, IndexType SlaveNodeId, VariableComponentType& SlaveVariable, double weight)
    {
        Node<3>& SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3>& MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
    	DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }

    // Functions with use two variables
    void AddMasterSlaveRelationWithNodesAndVariable(Node<3> &MasterNode, VariableType& MasterVariable, Node<3> &SlaveNode, VariableType& SlaveVariable, double weight)
    {
        SlaveNode.Set(SLAVE);        
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
    	DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }


    void AddMasterSlaveRelationWithNodeIdsAndVariable(IndexType MasterNodeId, VariableType& MasterVariable, IndexType SlaveNodeId, VariableType& SlaveVariable, double weight)
    {
        Node<3>& SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3>& MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
    	DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }



    // Default functions
    /**
		Applies the MPC condition using DOFs, one as master and other as slave, and with the given weight
		@arg slaveDOF 
        @arg masterDOF 
        @arg weight
		*/     
    void AddMasterSlaveRelationWithDofs(DofType slaveDOF, DofType masterDOF, double masterWeight, int PartitionId=0 )
    {
        pMpc->AddConstraint(slaveDOF, masterDOF,  masterWeight, PartitionId);
    }

    /**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/  
    void SetActive(bool isActive=true)
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
    virtual ~ApplyMultipointConstraintsProcess(){

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
        
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ApplyMultipointConstraintsProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override {rOStream << "ApplyMultipointConstraintsProcess";}

    /// Print object's data.
    void PrintData() {
        std::cout<<"Number of slave nodes :: "<<std::endl;
        pMpc->GetInfo();
    }



protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    ModelPart&                                 mr_model_part;
    MpcDataPointerType                         pMpc;

private:

    /// Assignment operator.
    ApplyMultipointConstraintsProcess& operator=(ApplyMultipointConstraintsProcess const& rOther){return *this;}

}; // Class MoveRotorProcess

};  // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
