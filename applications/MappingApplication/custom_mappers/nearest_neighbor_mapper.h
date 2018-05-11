//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

#if !defined(KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"


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

/// Nearest Neighbor Mapper
/** This class implements the Nearest Neighbor Mapping technique.
* Each node on the destination side gets assigned is's closest neighbor on the other side of the interface.
* In the mapping phase every node gets assigned the value of it's neighbor
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/

template<class TSparseSpace, class TDenseSpace>
class NearestNeighborMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:

    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of NearestNeighborMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestNeighborMapper);

    using BaseType = Mapper<TSparseSpace, TDenseSpace>;
    using MapperLocalSystemPointer = typename BaseType::MapperLocalSystemPointer;

    ///@}
    ///@name Life Cycle
    ///@{

    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination)
                          : Mapper<TSparseSpace, TDenseSpace>(rModelPartOrigin,
                                   rModelPartDestination) {}

    NearestNeighborMapper(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters)
                          : Mapper<TSparseSpace, TDenseSpace>(rModelPartOrigin,
                                   rModelPartDestination,
                                   JsonParameters)
    {
        // mpMapperCommunicator->InitializeOrigin(MapperUtilities::Node_Coords);
        // mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
        // mpMapperCommunicator->Initialize();

        // mpInverseMapper.reset(); // explicitly specified to be safe
    }

    /// Destructor.
    virtual ~NearestNeighborMapper() { }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void UpdateInterface(Kratos::Flags MappingOptions, double SearchRadius) override
    {
        // mpMapperCommunicator->UpdateInterface(MappingOptions, SearchRadius);
        // if (mpInverseMapper)
        // {
        //     mpInverseMapper->UpdateInterface(MappingOptions, SearchRadius);
        // }

        // if (MappingOptions.Is(MapperFlags::REMESHED))
        // {
        //     ComputeNumberOfNodesAndConditions();
        // }
    }

    typename Mapper<TSparseSpace, TDenseSpace>::Pointer Clone(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters) override
    {
        return Kratos::make_shared<NearestNeighborMapper<TSparseSpace, TDenseSpace>>(rModelPartOrigin,
                                                          rModelPartDestination,
                                                          JsonParameters);
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
        return "NearestNeighborMapper";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestNeighborMapper";
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

    void InitializeMapperLocalSystem(MapperLocalSystemPointer& pMapperLocalSystem) const override
    {
        // pMapperLocalSystem = Kratos::make_unique<NearestNeighborLocalSystem>();
    }

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // template <typename T>
    // static T GetValueOfNode(InterfaceObject::Pointer pInterfaceObject, //TODO const
    //                         const Variable< T >& rVariable,
    //                         const Kratos::Flags& rOptions,
    //                         const std::vector<double>& rShapeFunctionValues)
    // {
    //     Node<3>* p_base_node = pInterfaceObject->pGetBaseNode();

    //     return p_base_node->FastGetSolutionStepValue(rVariable);
    // }


    // template <typename T>
    // static void SetValueOfNode(InterfaceObject::Pointer pInterfaceObject,
    //                            const T& rValue,
    //                            const Variable< T >& rVariable,
    //                            const Kratos::Flags& rOptions,
    //                            const double Factor)
    // {
    //     Node<3>* p_base_node = pInterfaceObject->pGetBaseNode();

    //     if (rOptions.Is(MapperFlags::ADD_VALUES))
    //     {
    //         p_base_node->FastGetSolutionStepValue(rVariable) += rValue * Factor;
    //     }
    //     else
    //     {
    //         p_base_node->FastGetSolutionStepValue(rVariable) = rValue * Factor;
    //     }
    // }

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
    // NearestNeighborMapper& operator=(NearestNeighborMapper const& rOther) {}

    /// Copy constructor.
    //NearestNeighborMapper(NearestNeighborMapper const& rOther);

    ///@}

}; // Class NearestNeighborMapper

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function


///@}

///@} addtogroup block
}  // namespace Kratos.

#endif // KRATOS_NEAREST_NEIGHBOR_MAPPER_H_INCLUDED  defined