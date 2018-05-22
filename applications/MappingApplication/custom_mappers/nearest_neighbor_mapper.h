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

    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;

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
        // The Initialize function has to be called here bcs it internally calls virtual
        // functions that would not exist yet if it was called from the BaseClass!
        this->Initialize();


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

    MapperUniquePointerType Clone(ModelPart& rModelPartOrigin,
                          ModelPart& rModelPartDestination,
                          Parameters JsonParameters) override
    {
        return Kratos::make_unique<NearestNeighborMapper<TSparseSpace, TDenseSpace>>(rModelPartOrigin,
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

    class NearestNeigborInterfaceData
    {
    public:

        bool SendBack() const
        {
            return mSendingBackRequired;
        }

        void ProcessSearchResult(InterfaceObject::Pointer pInterfaceObject)
        {
            mSendingBackRequired = true; // If this function is called it means that the search was successful

            const double distance = 1.0; // TODO how to get the distance?
            if (distance < mNearestNeighborDistance)
            {
                mNearestNeighborDistance = distance;
                mNearestNeighborId = pInterfaceObject->pGetBaseNode()->GetValue(INTERFACE_EQUATION_ID);
            }
        };

        std::size_t GetNearestNeighborId() const
        {
            return mNearestNeighborId;
        }

        double GetNearestNeighborDistance() const
        {
            return mNearestNeighborDistance;
        }

    private:
        int mNearestNeighborId = -1; // default value, indicates an unsuccessful local search
        double mNearestNeighborDistance = std::numeric_limits<double>::max();
        bool mSendingBackRequired = false; // this is not being serialized since it is not needed after mpi-data-exchange!

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const
        {
            rSerializer.save("NearestNeighborId", mNearestNeighborId);
            rSerializer.save("NearestNeighborDistance", mNearestNeighborDistance);
        }

        virtual void load(Serializer& rSerializer)
        {
            rSerializer.load("NearestNeighborId", mNearestNeighborId);
            rSerializer.load("NearestNeighborDistance", mNearestNeighborDistance);
        }

    };

    template<class TDataHolder>
    class NearestNeighborLocalSystem : public MapperLocalSystem<TDataHolder>
    {
        public:
        using BaseType = MapperLocalSystem<TDataHolder>;
        using NodePointerType = typename BaseType::NodePointerType;

        using MappingWeightsVector = typename BaseType::MappingWeightsVector;
        using EquationIdVectorType = typename BaseType::EquationIdVectorType;

        using SizeType = typename BaseType::IndexType;
        using IndexType = typename BaseType::IndexType;

        NearestNeighborLocalSystem()
        {

        }

        NearestNeighborLocalSystem(NodePointerType pNode) : mpNode(pNode)
        {

        }

        Kratos::unique_ptr<BaseMapperLocalSystem> Create(NodePointerType pNode) const override
        {
            return Kratos::make_unique<NearestNeighborLocalSystem<TDataHolder>>(pNode);
        }

        void CalculateAll(MappingWeightsVector& rMappingWeights,
                          EquationIdVectorType& rOriginIds,
                          EquationIdVectorType& rDestinationIds) const override
        {
            if (rMappingWeights.size() != 1) rMappingWeights.resize(1);
            if (rOriginIds.size() != 1)      rOriginIds.resize(1);
            if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

            if (this->mInterfaceInfos.size() > 0)
            {
                IndexType nearest_neighbor_id = this->mInterfaceInfos[0]->GetInterfaceData().GetNearestNeighborId();
                double nearest_neighbor_distance = this->mInterfaceInfos[0]->GetInterfaceData().GetNearestNeighborDistance();

                for (SizeType i=1; i<this->mInterfaceInfos.size(); ++i)
                {
                    TDataHolder data = this->mInterfaceInfos[i]->GetInterfaceData();
                    double distance = data.GetNearestNeighborDistance();

                    if (distance < nearest_neighbor_distance)
                    {
                        nearest_neighbor_distance = distance;
                        nearest_neighbor_id = data.GetNearestNeighborId();
                    }
                }

                rMappingWeights[0] = 1.0;
                rOriginIds[0] = nearest_neighbor_id;
                // rDestinationIds[0] = this->mrNode.GetValue(INTERFACE_EQUATION_ID); //TODO
            }
            else
            {
                KRATOS_WARNING_IF("NearestNeighborMapper", this->mInterfaceInfos.size() == 0)
                    << "MapperLocalSystem No xxx" << "xxx" << " has not found a neighbor" << std::endl;

                // TODO is this ok? => I guess it would be better to do this in a mor general way in the baseclass...
                // TODO resize to zero, then it wont be assembled! (might be a bit slower though...)
                rMappingWeights[0] = 0.0;
                rOriginIds[0]      = 0;
                rDestinationIds[0] = 0;
            }

        }

        bool UseNodesAsBasis() const override { return true; }

        private:
        NodePointerType mpNode;

    };

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

    MapperLocalSystemPointer GetMapperLocalSystem() const override
    {
        return Kratos::make_unique<NearestNeighborLocalSystem<NearestNeigborInterfaceData>>();
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