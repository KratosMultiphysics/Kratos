//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Miguel Angel Celigueta
//

#if !defined(KRATOS_SPLIT_INTERNAL_INTERFACE_PROCESS_H_INCLUDED )
#define  KRATOS_SPLIT_INTERNAL_INTERFACE_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "processes/find_elements_neighbours_process.h"

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

/**
 * @class SplitInternalInterfaceProcess
 * @ingroup KratosCore
 * @brief Computes NODAL_AREA
 * @details splits a domain across changes of property and generates a condition at the splitting positions
 * @author Riccardo Rossi
 * @author Miguel Angel Celigueta
 */
class KRATOS_API(KRATOS_CORE) SplitInternalInterfaceProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Index type definition
    typedef std::size_t IndexType;

    /// Size type definition
    typedef std::size_t SizeType;

    /// The definition of the node
    typedef Node<3> NodeType;

    /// Pointer definition of SplitInternalInterfaceProcess
    KRATOS_CLASS_POINTER_DEFINITION(SplitInternalInterfaceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param DomainSize The size of the space, if the value is not provided will compute from the model part
     */
    SplitInternalInterfaceProcess(Model& rModel,
                                  Parameters rParameters
                                 )
        : Process(Flags()), mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
    {
        KRATOS_TRY


        const Parameters default_parameters = GetDefaultParameters();
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mConditionName = rParameters["condition_name"].GetString();
KRATOS_WATCH("constructed ----------------------------------------------------")
        KRATOS_CATCH("");
    }

    /// Destructor.
    virtual ~SplitInternalInterfaceProcess() override
    {
    }


    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters( R"(
        {
            "model_part_name" :"MODEL_PART_NAME",
            "condition_name"   : "CONDITION_NAME"
        }  )" );
        return default_parameters;
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExecuteInitialize() override
    {
        std::set< std::size_t> property_ids;




        for(auto& rElem : mrModelPart.Elements())
        {
            property_ids.insert( rElem.GetProperties().Id() );
        }

        std::size_t max_id = 0;
        std::size_t min_id = 1e6;
        for(auto id : property_ids)
        {
            max_id = std::max(id,max_id);
            min_id = std::min(id,min_id);
        }

        for(std::size_t id=min_id; id<max_id; id++)
        {
            std::cout << "processing property Id "  << id << std::endl;
            SplitBoundary(id, mrModelPart);
        }


    }






    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SplitInternalInterfaceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SplitInternalInterfaceProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{
    void SplitBoundary(const std::size_t PropertyIdBeingProcessed, ModelPart& rModelPart)
    {
        std::size_t domain_size = rModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension(); //TODO: this may not be a very good solution.
        FindElementalNeighboursProcess(rModelPart, domain_size).Execute();

        //construct list of faces on the interface
        std::vector< Geometry<Node<3> > > interface_faces;
        std::vector<
        std::pair<
        Geometry<Node<3>>::Pointer,
                 Geometry<Node<3>>::Pointer
                 >
                 > neighbouring_elements;

        for(auto& rElem : mrModelPart.Elements())
        {
            const auto& neighb = rElem.GetValue(NEIGHBOUR_ELEMENTS);

            for(unsigned int i=0; i<rElem.GetGeometry().size(); ++i)
            {

                if(rElem.GetProperties().Id() == PropertyIdBeingProcessed && neighb[i].GetProperties().Id() > PropertyIdBeingProcessed)
                {
                    auto boundary_entities = rElem.GetGeometry().GenerateBoundariesEntities();
                    interface_faces.push_back(boundary_entities[i]);
                    neighbouring_elements.push_back( std::make_pair(rElem.pGetGeometry()  ,neighb[i].pGetGeometry()  ) );
                }
            }
        }

        //construct list of nodes on the interface
        std::set< std::size_t > ids_on_interface;
        for(auto& geom : interface_faces)
        {
            for(auto& rNode : geom)
                ids_on_interface.insert(rNode.Id());
        }
        //create duplicated nodes list
        std::size_t max_node_id = 0;
        if(mrModelPart.GetRootModelPart().Nodes().size() != 0)
            max_node_id = (mrModelPart.GetRootModelPart().Nodes().end()-1)->Id() + 1;

        std::map<std::size_t, Node<3>::Pointer> new_nodes_map;
        for(auto& id : ids_on_interface)
        {
            auto& rOrigNode = rModelPart.Nodes()[id];
            auto pNode = mrModelPart.CreateNewNode(max_node_id++, rOrigNode );
            auto& origin_dofs = rOrigNode.GetDofs();
            for (auto it_dof = origin_dofs.begin(); it_dof != origin_dofs.end(); it_dof++)
            {
                pNode->pAddDof(**it_dof);
            }

            new_nodes_map[id] = pNode;
        }

        //now change the nodes to make the split and generate the new conditions
        std::size_t max_cond_id = 0;
        if(mrModelPart.GetRootModelPart().Conditions().size() != 0)
            max_cond_id = (mrModelPart.GetRootModelPart().Conditions().end()-1)->Id() + 1;

        Properties::Pointer pInterfaceProp = mrModelPart.pGetProperties(1); //TODO: understand if the property 1 is what we want
        for(std::size_t i=0; i<interface_faces.size(); ++i)
        {
            //do the split
            auto& pgeom = neighbouring_elements[i].second;
            for(std::size_t k=0; k<pgeom->size(); ++k)
            {
                auto it = new_nodes_map.find((*pgeom)[k].Id());
                if( it != new_nodes_map.end() )
                    (*pgeom)(i) = it->second;
            }
            //create prism
            std::vector<std::size_t> prism_ids;
            for(std::size_t k=0; k<interface_faces[i].size(); ++k)
                prism_ids.push_back(interface_faces[i][k].Id());
            for(std::size_t k=0; k<interface_faces[i].size(); ++k)
                prism_ids.push_back(new_nodes_map[interface_faces[i][k].Id()]->Id());

            rModelPart.CreateNewCondition(mConditionName, max_cond_id++, prism_ids, pInterfaceProp ); //TODO: understand if the property 1 is what we want
        }


for(auto& cond : mrModelPart.Elements())
    std::cout << cond << std::endl;


    }

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;  /// The model part where the nodal area is computed
    std::string mConditionName;    /// The dimension of the space

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{



    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SplitInternalInterfaceProcess& operator=(SplitInternalInterfaceProcess const& rOther);

    /// Copy constructor.
    //SplitInternalInterfaceProcess(SplitInternalInterfaceProcess const& rOther);


    ///@}

}; // Class SplitInternalInterfaceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SplitInternalInterfaceProcess& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SplitInternalInterfaceProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SPLIT_INTERNAL_INTERFACE_PROCESS_H_INCLUDED  defined 


