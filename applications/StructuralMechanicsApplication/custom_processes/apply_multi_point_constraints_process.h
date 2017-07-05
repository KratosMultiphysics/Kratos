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

//#ifdef MAPPING_APPLICATION
#include "../../MappingApplication/custom_utilities/mapper.h"
#include "../../MappingApplication/custom_utilities/mapper_communicator.h"
//#endif

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
    typedef unsigned int IndexType;
    typedef std::vector<MpcDataPointerType> *MpcDataPointerVectorType;
    typedef MpcData::VariableType VariableType;
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
                "interpolation_type":"nearest_node",
                "reform_every_step":false   
            }  )");

        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if (info->GetValue(MPC_DATA_CONTAINER) == NULL)
            info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

        pMpc = MpcDataPointerType(new MpcData());
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
        
        AddMasterSlaveRelation();        
    }

    ApplyMultipointConstraintsProcess(ModelPart &model_part, std::string name = "default") : Process(Flags()), mr_model_part(model_part), m_parameters("{}")
    {
        ProcessInfoPointerType info = mr_model_part.pGetProcessInfo();
        if (info->GetValue(MPC_DATA_CONTAINER) == NULL)
            info->SetValue(MPC_DATA_CONTAINER, new std::vector<MpcDataPointerType>());

        pMpc = MpcDataPointerType(new MpcData());
        pMpc->SetName(name);
        pMpc->SetActive(true);

        MpcDataPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
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
        std::string interpolationType = m_parameters["interpolation_type"].GetString();
        int numVars = m_parameters["variable_names"].size();
        Parameters mapper_parameters = m_parameters["interpolation_settings"];

        MapperCommunicator::Pointer mpMapperCommunicator = MapperCommunicator::Pointer(
            new MapperCommunicator(master_model_part,
                                   slave_model_part,
                                   mapper_parameters));

        // KratosComponents< Variable<double> >::Get( mvariable_name )
        if (interpolationType == "nearest_node")
        {
            mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
            mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
            mpMapperCommunicator->Initialize();
        }
        else if (interpolationType == "nearest_element")
        {
            mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Center);
            mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
            mpMapperCommunicator->Initialize();
        }

        for (int i = 0; i < numVars; i++)
        {

            VariableComponentType rVar = KratosComponents<VariableComponentType>::Get(m_parameters["variable_names"][i].GetString());

            // Create the mapper based on the type of interpolation
            // Creating the function pointers for the InterfaceObjects
            if (interpolationType == "nearest_node")
            {
                auto function_pointer_origin = std::bind(&GetMasterRelationInformationFromNode,
                                                         std::placeholders::_1,
                                                         rVar,
                                                         std::placeholders::_2);

                auto function_pointer_destination = std::bind(&SetMpcDataAtNode<void *>,
                                                              std::placeholders::_1,
                                                              rVar,
                                                              std::placeholders::_2,
                                                              pMpc);

                mpMapperCommunicator->TransferVariableData(function_pointer_origin,
                                                           function_pointer_destination);
            }
            else if (interpolationType == "nearest_element")
            {
                auto function_pointer_origin = std::bind(&GetMasterRelationInformationFromElement,
                                                         std::placeholders::_1,
                                                         rVar,
                                                         std::placeholders::_2);
                auto function_pointer_destination = std::bind(&SetMpcDataAtNode<void *>,
                                                              std::placeholders::_1,
                                                              rVar,
                                                              std::placeholders::_2,
                                                              pMpc);

                mpMapperCommunicator->TransferVariableData(function_pointer_origin,
                                                           function_pointer_destination);
            }
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
    void AddMasterSlaveRelationWithNodesAndVariableComponents(Node<3> &MasterNode, VariableComponentType &MasterVariable, Node<3> &SlaveNode, VariableComponentType &SlaveVariable, double weight)
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariableComponents(IndexType MasterNodeId, VariableComponentType &MasterVariable, IndexType SlaveNodeId, VariableComponentType &SlaveVariable, double weight)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }

    // Functions with use two variables
    void AddMasterSlaveRelationWithNodesAndVariable(Node<3> &MasterNode, VariableType &MasterVariable, Node<3> &SlaveNode, VariableType &SlaveVariable, double weight)
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, 0);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariable(IndexType MasterNodeId, VariableType &MasterVariable, IndexType SlaveNodeId, VariableType &SlaveVariable, double weight)
    {
        Node<3> &SlaveNode = mr_model_part.Nodes()[SlaveNodeId];
        Node<3> &MasterNode = mr_model_part.Nodes()[MasterNodeId];
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
    void AddMasterSlaveRelationWithDofs(DofType slaveDOF, DofType masterDOF, double masterWeight, int PartitionId = 0)
    {
        pMpc->AddConstraint(slaveDOF, masterDOF, masterWeight, PartitionId);
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

        if (m_parameters["reform_every_step"].GetBool())
            // Adding the master slave relation between the master and slave sub model parts
            AddMasterSlaveRelation();

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
        std::cout << "Number of slave nodes :: " << std::endl;
        pMpc->GetInfo();
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

    /*
    *   Structrue which contain the slave master information. This is used in conjection with the mapper communicator
    *   
    */
    struct MasterSlaveRelation
    {
        std::vector<int> MastersDOFIds;
        std::vector<double> MastersDOFWeights;
    };

    /*
    * Function to be used in realation with the nearest node mapper. Master side 
    */
    static MasterSlaveRelation *GetMasterRelationInformationFromNode(InterfaceObject *pInterfaceObject, const VariableComponentType &rVariable,
                                                                     const std::vector<double> &rShapeFunctionValues)
    {
        MasterSlaveRelation *pMasterSlaveRelation = new MasterSlaveRelation();
        Node<3> *p_base_node = static_cast<InterfaceNode *>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;

        unsigned int dofId = p_base_node->GetDof(rVariable).EquationId();
        pMasterSlaveRelation->MastersDOFIds.push_back(dofId);
        pMasterSlaveRelation->MastersDOFWeights.push_back(1.0);

        return pMasterSlaveRelation;
    }

    /*
    * Function to be used in realation with the nearest element mapper. Master side 
    */
    static MasterSlaveRelation *GetMasterRelationInformationFromElement(InterfaceObject *pInterfaceObject, const VariableComponentType &rVariable,
                                                                        const std::vector<double> &rShapeFunctionValues)
    {
        MasterSlaveRelation *pMasterSlaveRelation = new MasterSlaveRelation();
        Geometry<Node<3>> *p_base_geometry = static_cast<InterfaceGeometryObject *>(pInterfaceObject)->pGetBase();
        KRATOS_ERROR_IF_NOT(p_base_geometry) << "Base Pointer is nullptr!!!" << std::endl;

        for (std::size_t i = 0; i < p_base_geometry->PointsNumber(); ++i)
        {
            unsigned int dofId = p_base_geometry->GetPoint(i).GetDof(rVariable).EquationId();
            pMasterSlaveRelation->MastersDOFIds.push_back(dofId);
            pMasterSlaveRelation->MastersDOFWeights.push_back(rShapeFunctionValues[i]);
        }

        return pMasterSlaveRelation;
    }

    /*
    * Function to be used in realation with the mapper. Slave side 
    */
    template <typename T>
    static void SetMpcDataAtNode(InterfaceObject *pInterfaceObject, VariableComponentType &rVariable, T rValue, MpcDataPointerType pMpc)
    {

        Node<3> *p_base_node = static_cast<InterfaceNode *>(pInterfaceObject)->pGetBase();
        MasterSlaveRelation *mMasterSlaveRelation = static_cast<MasterSlaveRelation *>(rValue);
        KRATOS_ERROR_IF_NOT(p_base_node) << "Base Pointer is nullptr!!!" << std::endl;
        // Marking the node as a slave
        p_base_node->Set(SLAVE);

        unsigned int slaveDofId = p_base_node->GetDof(rVariable).EquationId();
        for (int i = 0; i < mMasterSlaveRelation->MastersDOFIds.size(); i++)
        {
            pMpc->AddConstraint(slaveDofId, mMasterSlaveRelation->MastersDOFIds[i], mMasterSlaveRelation->MastersDOFWeights[i], 0);
        }

        delete mMasterSlaveRelation;
    }

}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
