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
#include "spaces/ublas_space.h"
#include "utilities/constraint.h"
#include "utilities/multipoint_constraint.h"

namespace Kratos
{

class ApplyMultipointConstraintsProcess : public Process
{
  public:
    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyMultipointConstraintsProcess);

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType; // TODO: pass them also as templates
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef Constraint<SparseSpaceType, LocalSpaceType>::Pointer ConstraintPointerType;
    typedef MultipointConstraint<SparseSpaceType, LocalSpaceType>::Pointer MpcPointerType;
    typedef Dof<double> *DofPointerType;
    typedef Dof<double> DofType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef ProcessInfo ProcessInfoType;
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef unsigned int IndexType;
    typedef Kratos::Variable<double> VariableType;
    typedef boost::shared_ptr<std::vector<ConstraintPointerType>> ConstraintSharedPointerVectorType;
    typedef ModelPart::NodeIterator NodeIterator;

    /// Constructor.
    ApplyMultipointConstraintsProcess(ModelPart &rModelPart,
                                      Parameters parameters) : Process(Flags()), mrModelPart(rModelPart), mParameters(parameters)
    {

        Parameters default_parameters(R"(
            {
                "constraint_set_name":"default",
                "master_sub_model_part_name":"default_master",
                "slave_sub_model_part_name":"default_slave",                
                "variable_names":[],
                "reform_every_step":false   
            }  )");
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        if(mParameters["variable_names"].size() == 0)
            KRATOS_THROW_ERROR("","In ApplyMultipointConstraintsProcess class constructor :: No DOFs specified for applying Multipoint constraints."," ")


        ProcessInfoPointerType info = mrModelPart.pGetProcessInfo();
        if (!info->Has(CONSTRAINTS_CONTAINER))    
           info->SetValue(CONSTRAINTS_CONTAINER, Kratos::make_shared<std::vector<ConstraintPointerType>>());


        pMpc = Kratos::make_shared<MultipointConstraint<SparseSpaceType,LocalSpaceType>>();
        std::string name = mParameters["constraint_set_name"].GetString();
        pMpc->SetName(name);
        pMpc->SetActive(true);

        ConstraintSharedPointerVectorType mpc_data_vector = info->GetValue(CONSTRAINTS_CONTAINER);
        (*mpc_data_vector).push_back(pMpc);

        if (!mParameters["reform_every_step"].GetBool())
            // Adding the master slave relation between the master and slave sub model parts
            AddMasterSlaveRelation();
    }

    ApplyMultipointConstraintsProcess(ModelPart &rModelPart, std::string rName = "default") : Process(Flags()), mrModelPart(rModelPart), mParameters("{}")
    {
        Parameters default_parameters(R"(
            {
                "constraint_set_name":"default",
                "master_sub_model_part_name":"default_master",
                "slave_sub_model_part_name":"default_slave",                
                "variable_names":[],
                "reform_every_step":false   
            }  )");
        mParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

        // IMPORTANT : This constructor is not to be used when using this process in the normal KRATOS process_list of python script
        ProcessInfoPointerType info = mrModelPart.pGetProcessInfo();
        if (!info->Has(CONSTRAINTS_CONTAINER))
             info->SetValue(CONSTRAINTS_CONTAINER, Kratos::make_shared<std::vector<ConstraintPointerType>>());

        pMpc = Kratos::make_shared<MultipointConstraint<SparseSpaceType,LocalSpaceType>>();
        pMpc->SetName(rName);
        pMpc->SetActive(true);

        ConstraintSharedPointerVectorType mpc_data_vector = info->GetValue(CONSTRAINTS_CONTAINER);
        (*mpc_data_vector).push_back(pMpc);
    }

    /**
		Applies the MPC condition using two model parts, one as master and other as slave.
        Here a nearest element interpolation is used by default to get the relation between master and slave
		*/
    void AddMasterSlaveRelation()
    {
        ModelPart &master_model_part = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());

        ProcessInfoPointerType info = mrModelPart.pGetProcessInfo();
        int &dim = info->GetValue(DOMAIN_SIZE);

        //std::string interpolationType = mParameters["interpolation_type"].GetString();
        Parameters mapper_parameters = mParameters["interpolation_settings"];

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
        std::unique_ptr< BinBasedFastPointLocator<TDim> > p_point_locator (new BinBasedFastPointLocator<TDim>(master_model_part));
        int num_vars = mParameters["variable_names"].size();
        // iterating over slave nodes to find the corresponding masters
        const int n_slave_nodes = slave_model_part.Nodes().size();
        array_1d<double, TDim + 1> N; // This is only for triangular meshes
        const int max_results = 100;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        ModelPart::NodesContainerType::iterator iparticle = slave_model_part.NodesBegin();
        typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();        

        for (int i = 0; i < n_slave_nodes; i++)
        {            
            iparticle += i;
            Node<3>::Pointer p_slave_node = *(iparticle.base());
            Element::Pointer p_master_element;
            bool is_found = p_point_locator->FindPointOnMesh(p_slave_node->Coordinates(), N, p_master_element, result_begin, max_results);
            if (is_found == true)
            {
                for (int i = 0; i < num_vars; i++)
                {
                    std::string var_name = mParameters["variable_names"][i].GetString();
                    Geometry<Node<3>> &r_geometry = p_master_element->GetGeometry();
                    for (unsigned int i = 0; i < r_geometry.size(); i++)
                    {
                        if (KratosComponents<Variable<double>>::Has(var_name))
                        { //case of double variable
                            VariableType r_var = KratosComponents<Variable<double>>::Get(mParameters["variable_names"][i].GetString());
                            this->AddMasterSlaveRelationWithNodesAndVariable(r_geometry[i], r_var, *p_slave_node, r_var, N[i], 0.0);
                        }
                        else if (KratosComponents<VariableComponent<VectorComponentAdaptor<array_1d<double, 3>>>>::Has(var_name))
                        {
                            VariableComponentType r_var = KratosComponents<VariableComponentType>::Get(mParameters["variable_names"][i].GetString());
                            this->AddMasterSlaveRelationWithNodesAndVariableComponents(r_geometry[i], r_var, *p_slave_node, r_var, N[i], 0.0);
                        }
                    }
                }
            }
        }
    }

    /**
		Applies the MPC condition using two nodes, one as master and other as slave, and with the given Weight
		@arg rMasterNode 
        @arg rMasterVariable 
        @arg rSlaveNode 
        @arg rSlaveVariable
        @arg Weight
		*/
    void AddMasterSlaveRelationWithNodesAndVariableComponents(Node<3> &rMasterNode, VariableComponentType &rMasterVariable, Node<3> &rSlaveNode, VariableComponentType &rSlaveVariable, double Weight, double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        DofType &r_slave_dof = rSlaveNode.GetDof(rSlaveVariable);
        DofType &r_master_dof = rMasterNode.GetDof(rMasterVariable);
        AddMasterSlaveRelationWithDofs(r_slave_dof, r_master_dof, Weight, Constant);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariableComponents(IndexType MasterNodeId, VariableComponentType &rMasterVariable, IndexType SlaveNodeId, VariableComponentType &rSlaveVariable, double Weight, double Constant = 0.0)
    {
        Node<3> &r_slave_node = mrModelPart.Nodes()[SlaveNodeId];
        Node<3> &r_master_node = mrModelPart.Nodes()[MasterNodeId];
        r_slave_node.Set(SLAVE);
        DofType &r_slave_dof = r_slave_node.GetDof(rSlaveVariable);
        DofType &r_master_dof = r_master_node.GetDof(rMasterVariable);
        AddMasterSlaveRelationWithDofs(r_slave_dof, r_master_dof, Weight, Constant);
    }

    // Functions with use two variables
    void AddMasterSlaveRelationWithNodesAndVariable(Node<3> &rMasterNode, VariableType &rMasterVariable, Node<3> &rSlaveNode, VariableType &rSlaveVariable, double Weight, double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        DofType &r_slave_dof = rSlaveNode.GetDof(rSlaveVariable);
        DofType &r_master_dof = rMasterNode.GetDof(rMasterVariable);
        AddMasterSlaveRelationWithDofs(r_slave_dof, r_master_dof, Weight, Constant);
    }

    void AddMasterSlaveRelationWithNodeIdsAndVariable(IndexType MasterNodeId, VariableType &rMasterVariable, IndexType SlaveNodeId, VariableType &rSlaveVariable, double Weight, double Constant = 0.0)
    {
        Node<3> &r_slave_node = mrModelPart.Nodes()[SlaveNodeId];
        Node<3> &r_master_node = mrModelPart.Nodes()[MasterNodeId];
        r_slave_node.Set(SLAVE);
        DofType &r_slave_dof = r_slave_node.GetDof(rSlaveVariable);
        DofType &r_master_dof = r_master_node.GetDof(rMasterVariable);
        AddMasterSlaveRelationWithDofs(r_slave_dof, r_master_dof, Weight, Constant);
    }

    // Default functions
    /**
		Applies the MPC condition using DOFs, one as master and other as slave, and with the given Weight
		@arg slaveDOF 
        @arg masterDOF 
        @arg Weight
		*/
    void AddMasterSlaveRelationWithDofs(DofType& rSlaveDOF, DofType& rMasterDOF, double MasterWeight, double Constant = 0.0)
    {
        pMpc->AddConstraint(rSlaveDOF, rMasterDOF, MasterWeight, Constant);
    }

    /**
		Activates the constraint set or deactivates
		@arg isActive true/false
		*/
    void SetActive(bool IsActive = true)
    {
        pMpc->SetActive(IsActive);
    }

    /**
		Sets the name of the constraint set
		@arg Name Name of the constraint set
		*/
    void SetName(std::string Name)
    {
        pMpc->SetName(Name);
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

        if (mParameters["reform_every_step"].GetBool())
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
    ModelPart &mrModelPart;
    MpcPointerType pMpc;
    Parameters mParameters;

  private:
    /// Assignment operator.
    ApplyMultipointConstraintsProcess &operator=(ApplyMultipointConstraintsProcess const &rOther) { return *this; }

}; // Class MoveRotorProcess

}; // namespace Kratos.

#endif // 
