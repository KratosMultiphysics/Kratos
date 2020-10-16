//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt

#if !defined(KRATOS_PML_DIRECTION_PROCESS_H_INCLUDED)
#define KRATOS_MONOLITHIC_MAPPING_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "mor_application_variables.h"

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

/**
 * @class MonolithicMappingProcess
 *
 * @ingroup MORApplication
 *
 * @brief This method provides mapping info to conditions
 * @details Based on the provided mapping matrix, the nodes of two model parts are paired
*/
template<typename MatrixType>
class MonolithicMappingProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;

    /// Pointer definition of MonolithicMappingProcess
    KRATOS_CLASS_POINTER_DEFINITION(MonolithicMappingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MonolithicMappingProcess(
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        MatrixType& rMappingMatrix
        ) : mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mrMappingMatrix(rMappingMatrix)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~MonolithicMappingProcess() override = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
        std::cout << "execute nice process!!\n";
        // KRATOS_WATCH(mrMappingMatrix)
        const std::size_t n_nodes_dest = mrMappingMatrix.size1();
        std::map<std::size_t, std::size_t> matrix_id_map;
        for( std::size_t i=0; i<mrModelPartOrigin.Nodes().size(); ++i ) {
            std::size_t id = mrModelPartOrigin.NodesArray()[i]->Id();
            matrix_id_map[id] = i;
            // std::cout << "inserting id=" << id << "i=" << i << std::endl;
            // matrix_id_set.insert(std::make_pair(id, i));
        }

        for( auto& condition : mrModelPartOrigin.Conditions() ) {
            const std::size_t n_nodes = condition.GetGeometry().size();
            NodeVector destination_nodes;
            std::vector<std::size_t> destination_nodes_matrix_id;
            std::vector<std::size_t> origin_nodes_matrix_id;

            for( std::size_t i=0; i<n_nodes_dest; ++i ) {
                double tmp = 0;
                // KRATOS_WATCH(condition.GetGeometry()[i])
                // std::cout << matrix_id_map[condition.GetGeometry()[i].Id()] << std::endl;
                for( std::size_t j=0; j<n_nodes; ++j ){
                    const size_t this_id = matrix_id_map[condition.GetGeometry()[j].Id()];
                    tmp += mrMappingMatrix(i,this_id);
                    origin_nodes_matrix_id.push_back(this_id);
                }

                //append to list of "geometry nodes" for the condition
                if( std::abs(1-tmp) < std::numeric_limits<double>::epsilon()*10000 ) {
                    const std::size_t dest_id = mrModelPartDestination.NodesArray()[i]->Id();
                    destination_nodes.push_back(mrModelPartDestination.pGetNode(dest_id));
                    destination_nodes_matrix_id.push_back(i);
                }

            }
            // KRATOS_WATCH(nodes)
            condition.SetValue(MAPPING_NODES, destination_nodes);

            //now write the mapping factors
            const size_t n_destination_nodes = destination_nodes.size();
            Matrix mapping_factors(n_destination_nodes, n_nodes);
            for( std::size_t i=0; i<n_destination_nodes; ++i ) {
                for( std::size_t j=0; j<n_nodes; ++j ) {
                    mapping_factors(i,j) = mrMappingMatrix(destination_nodes_matrix_id[i], origin_nodes_matrix_id[j]);
                }
            }
            // KRATOS_WATCH(mapping_factors)
            condition.SetValue(MAPPING_FACTOR, mapping_factors);
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
        return "MonolithicMappingProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MonolithicMappingProcess";
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

    ModelPart& mrModelPartOrigin;              // The main model part
    ModelPart& mrModelPartDestination;
    MatrixType& mrMappingMatrix;

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
    ///@na
///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MonolithicMappingProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MonolithicMappingProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MonolithicMappingProcess& operator=(MonolithicMappingProcess const& rOther) = delete;

    /// Copy constructor.
    MonolithicMappingProcess(MonolithicMappingProcess const& rOther) = delete;


    ///@}

}; // Class MonolithicMappingProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   MonolithicMappingProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const MonolithicMappingProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }

}
#endif /* KRATOS_MONOLITHIC_MAPPING_PROCESS_H_INCLUDED defined */
