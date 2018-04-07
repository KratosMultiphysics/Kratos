//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#ifndef APPLY_MULTI_POINT_CONSTRAINTS_PROCESS_H
#define APPLY_MULTI_POINT_CONSTRAINTS_PROCESS_H

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "utilities/binbased_fast_point_locator.h"

// Application includes
#include "custom_utilities/multipoint_constraint_data.hpp"

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
class ApplyMultipointConstraintsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyMultipointConstraintsProcess);

    /// MPC definitions
    typedef MpcData::Pointer MpcDataPointerType;
    typedef MpcData::VariableComponentType VariableComponentType;
    typedef MpcData::VariableType VariableType;
    typedef std::map<std::string, MpcDataPointerType> MpcDataMapType;
    typedef std::vector<MpcDataPointerType> MpcDataPointerVectorType;
    typedef Kratos::shared_ptr<MpcDataPointerVectorType> MpcDataSharedPointerVectorType;

    /// Dof definitions
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;

    // Process info definitions
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;

    /// Index definition
    typedef std::size_t IndexType;

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
    /**
     * @brief Parameters constructor. In this constructor we pass the parameters
     * @param rModelPart The model part to be computed
     * @param rParameters The configuration parameters
     */
    ApplyMultipointConstraintsProcess(
        ModelPart &rModelPart,
        Parameters rParameters
        ) : Process(Flags()),
            mrModelPart(rModelPart),
            mParameters(rParameters)
    {
        // We assign the default parameters if not defined
        Parameters default_parameters(R"(
        {
            "constraint_set_name":"default",
            "master_sub_rModelPart_name":"default_master",
            "slave_sub_rModelPart_name":"default_slave",
            "variable_names":[""],
            "reform_every_step":false
        }  )");

        mParameters.ValidateAndAssignDefaults(default_parameters);

        // We create the MPC data container in case is not defined or replace it in case is a null pointer
        ProcessInfoPointerType info = mrModelPart.pGetProcessInfo();
        if (info->Has(MPC_DATA_CONTAINER) ) {
            if ( info->GetValue(MPC_DATA_CONTAINER) == nullptr)
                info->SetValue(MPC_DATA_CONTAINER, MpcDataSharedPointerVectorType(Kratos::make_shared<MpcDataPointerVectorType>()));
        } else {
            info->SetValue(MPC_DATA_CONTAINER, MpcDataSharedPointerVectorType(Kratos::make_shared<MpcDataPointerVectorType>()));
        }

        mpMpc = MpcDataPointerType(Kratos::make_shared<MpcData>());
        std::string name = rParameters["constraint_set_name"].GetString();
        mpMpc->SetName(name);
        mpMpc->SetActive(true);

        MpcDataSharedPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        mpcDataVector->push_back(mpMpc);

        // Adding the master slave relation between the master and slave sub model parts
        if (!mParameters["reform_every_step"].GetBool()) {
            AddMasterSlaveRelation();
        }
    }

    /**
     * @brief Name constructor. In this constructor we pass the name of the MPC
     * @param rModelPart The model part to be computed
     * @param Name The name of the MPC data container
     */
    ApplyMultipointConstraintsProcess(
        ModelPart &rModelPart,
        std::string Name = "default"
        ) : Process(Flags()),
            mrModelPart(rModelPart),
            mParameters("{}")
    {
        // IMPORTANT : This constructor is not to be used when using this process in the normal KRATOS process_list of python script
        ProcessInfoPointerType info = mrModelPart.pGetProcessInfo();
        if (info->Has(MPC_DATA_CONTAINER)) {
            if( info->GetValue(MPC_DATA_CONTAINER) == nullptr)
                info->SetValue(MPC_DATA_CONTAINER, MpcDataSharedPointerVectorType(Kratos::make_shared<MpcDataPointerVectorType>()));
        } else {
            info->SetValue(MPC_DATA_CONTAINER, MpcDataSharedPointerVectorType(Kratos::make_shared<MpcDataPointerVectorType>()));
        }

        mpMpc = MpcDataPointerType(Kratos::make_shared<MpcData>());
        mpMpc->SetName(Name);
        mpMpc->SetActive(true);

        MpcDataSharedPointerVectorType mpcDataVector = info->GetValue(MPC_DATA_CONTAINER);
        mpcDataVector->push_back(mpMpc);
    }

    /// Destructor.
    ~ApplyMultipointConstraintsProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{
    /**
     * Applies the MPC condition using two model parts, one as master and other as slave.
     * Here a nearest element interpolation is used by default to get the relation between master and slave
    */
    void AddMasterSlaveRelation()
    {
        ModelPart& r_master_model_part = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart& r_slave_model_part = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());

        ProcessInfoPointerType info = mrModelPart.pGetProcessInfo();
        const SizeType dim = info->GetValue(DOMAIN_SIZE);

        std::string interpolationType = mParameters["interpolation_type"].GetString();
        Parameters mapper_parameters = mParameters["interpolation_settings"];

        if (dim == 2) {
            ApplyConstraints<2>(r_master_model_part, r_slave_model_part);
        } else if (dim == 3) {
            ApplyConstraints<3>(r_master_model_part, r_slave_model_part);
        }
    }

    // Functions which use two variable components
    template <int TDim>
    void ApplyConstraints(
        ModelPart &rMasterModelPart,
        ModelPart &rSlaveModelPart
        )
    {
        BinBasedFastPointLocator<TDim> *p_point_locator = new BinBasedFastPointLocator<TDim>(rMasterModelPart);
        const SizeType num_vars = mParameters["variable_names"].size();
        // iterating over slave nodes to find the corresponding masters
        const SizeType n_slave_nodes = rSlaveModelPart.Nodes().size();
        array_1d<double, TDim + 1> N; // This is only for triangular meshes
        const int max_results = 100;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

        for (IndexType i = 0; i < n_slave_nodes; i++) {
            ModelPart::NodesContainerType::iterator iparticle = rSlaveModelPart.NodesBegin() + i;
            NodeType::Pointer p_slave_node = *(iparticle.base());
            typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
            Element::Pointer p_master_element;
            const bool is_found = p_point_locator->FindPointOnMesh(p_slave_node->Coordinates(), N, p_master_element, result_begin, max_results);
            if (is_found == true) {
                for (IndexType i = 0; i < num_vars; i++) {
                    std::string varName = mParameters["variable_names"][i].GetString();
                    Geometry<NodeType> &geom = p_master_element->GetGeometry();
                    for (IndexType i = 0; i < geom.size(); i++) {
                        if (KratosComponents<Variable<double>>::Has(varName)) { // Case of double variable
                            VariableType rVar = KratosComponents<Variable<double>>::Get(mParameters["variable_names"][i].GetString());
                            this->AddMasterSlaveRelationWithNodesAndVariable(geom[i], rVar, *p_slave_node, rVar, N[i], 0.0);
                        } else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(varName)) { // Case of an array variable
                            VariableComponentType rVar = KratosComponents<VariableComponentType>::Get(mParameters["variable_names"][i].GetString());
                            this->AddMasterSlaveRelationWithNodesAndVariableComponents(geom[i], rVar, *p_slave_node, rVar, N[i], 0.0);
                        }
                    }
                }
            }
        }
        delete p_point_locator;
    }

    /**
     * Applies the MPC condition using two nodes, one as master and other as slave, and with the given weight
     * @param MasterNode
     * @param MasterVariable
     * @param SlaveNode
     * @param SlaveVariable
     * @param weight
    */
    void AddMasterSlaveRelationWithNodesAndVariableComponents(
        NodeType &MasterNode,
        VariableComponentType &MasterVariable,
        NodeType &SlaveNode,
        VariableComponentType &SlaveVariable,
        const double weight,
        const double constant = 0.0
        )
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariableComponents(
        IndexType MasterNodeId,
        VariableComponentType &MasterVariable,
        IndexType SlaveNodeId,
        VariableComponentType &SlaveVariable,
        const double weight,
        const double constant = 0.0
        )
    {
        NodeType &SlaveNode = mrModelPart.Nodes()[SlaveNodeId];
        NodeType &MasterNode = mrModelPart.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    // Functions with use two variables
    void AddMasterSlaveRelationWithNodesAndVariable(
        NodeType &MasterNode,
        VariableType &MasterVariable,
        NodeType &SlaveNode,
        VariableType &SlaveVariable,
        const double weight,
        const double constant = 0.0
        )
    {
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariable(
        IndexType MasterNodeId,
        VariableType &MasterVariable,
        IndexType SlaveNodeId,
        VariableType &SlaveVariable,
        const double weight,
        const double constant = 0.0
        )
    {
        NodeType &SlaveNode = mrModelPart.Nodes()[SlaveNodeId];
        NodeType &MasterNode = mrModelPart.Nodes()[MasterNodeId];
        SlaveNode.Set(SLAVE);
        DofType &pointerSlaveDOF = SlaveNode.GetDof(SlaveVariable);
        DofType &pointerMasterDOF = MasterNode.GetDof(MasterVariable);
        AddMasterSlaveRelationWithDofs(pointerSlaveDOF, pointerMasterDOF, weight, constant);
    }

    // Default functions
    /**
     * Applies the MPC condition using DOFs, one as master and other as slave, and with the given weight
     * @param slaveDOF
     * @param masterDOF
     * @param masterWeight
    */
    void AddMasterSlaveRelationWithDofs(
        DofType slaveDOF,
        DofType masterDOF,
        const double masterWeight,
        const double constant = 0.0
        )
    {
        mpMpc->AddConstraint(slaveDOF, masterDOF, masterWeight, constant);
    }

    /**
     * Activates the constraint set or deactivates
     * @param isActive true/false
    */
    void SetActive(bool isActive = true)
    {
        mpMpc->SetActive(isActive);
    }

    /**
     * Sets the name of the constraint set
     * @param isActive true/false
    */
    void SetName(std::string name)
    {
        mpMpc->SetName(name);
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

        if (mParameters["reform_every_step"].GetBool())
            // Adding the master slave relation between the master and slave sub model parts
            AddMasterSlaveRelation();

        KRATOS_CATCH("");
    }

    void ExecuteAfterOutputStep() override
    {
        Clear();
    }

    ///@}
    ///@name Input and output
    ///@{

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
        mpMpc->GetInfo();
    }

    /// Print object's data.
    void Clear()
    {
        mpMpc->Clear();
    }

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart &mrModelPart;   /// The main model part where the MPC is computed
    MpcDataPointerType mpMpc; /// The MPC data container
    Parameters mParameters;   /// The parameters of the problem

    ///@}

private:
    ///@name Private Operators
    ///@{

    /// Assignment operator.
    ApplyMultipointConstraintsProcess &operator=(ApplyMultipointConstraintsProcess const &rOther) { return *this; }

    ///@}

}; // Class MoveRotorProcess
///@}
///@name Input and output
///@{

}; // namespace Kratos.

#endif // KRATOS_MOVE_ROTOR_PROCESS_H
