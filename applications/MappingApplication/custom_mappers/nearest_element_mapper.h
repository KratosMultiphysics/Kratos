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

#if !defined(KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED )
#define  KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

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

class NearestElementInterfaceInfo : public MapperInterfaceInfo
{
public:

    NearestElementInterfaceInfo(const Point rPoint,
                                const IndexType SourceLocalSystemIndex,
                                const IndexType SourceRank=0)
        : MapperInterfaceInfo(rPoint, SourceLocalSystemIndex, SourceRank)
    {

    }

    MapperInterfaceInfo::Pointer Create(const Point rPoint,
                                        const IndexType SourceLocalSystemIndex,
                                        const IndexType SouceRank) const override
    {
        return Kratos::make_shared<NearestElementInterfaceInfo>(rPoint,
                                                                SourceLocalSystemIndex,
                                                                SouceRank);
    }

    void ProcessSearchResult(InterfaceObject::Pointer pInterfaceObject, const double NeighborDistance) override;

    void GetValue(std::vector<std::size_t>& rValue, const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mNodeIds;
    }

    void GetValue(std::vector<double>& rValue, const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mShapeFunctionValues;
    }

    void GetValue(double& rValue, const InfoType ValueType=MapperInterfaceInfo::InfoType::Dummy) const override
    {
        rValue = mClosestProjectionDistance;
    }

private:

    std::vector<std::size_t> mNodeIds;
    std::vector<double> mShapeFunctionValues;
    double mClosestProjectionDistance = std::numeric_limits<double>::max();

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.save("NodeIds", mNodeIds);
        rSerializer.save("SFValues", mShapeFunctionValues);
        rSerializer.save("ClosestProjectionDistance", mClosestProjectionDistance);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, MapperInterfaceInfo );
        rSerializer.load("NodeIds", mNodeIds);
        rSerializer.load("SFValues", mShapeFunctionValues);
        rSerializer.load("ClosestProjectionDistance", mClosestProjectionDistance);
    }

};

class NearestElementLocalSystem : public MapperLocalSystem
{
public:
    using BaseType = MapperLocalSystem;
    using MapperLocalSystemUniquePointer = typename BaseType::MapperLocalSystemUniquePointer;
    using NodePointerType = typename BaseType::NodePointerType;

    using MappingWeightsVector = typename BaseType::MappingWeightsVector;
    using EquationIdVectorType = typename BaseType::EquationIdVectorType;

    using SizeType = typename BaseType::IndexType;
    using IndexType = typename BaseType::IndexType;

    NearestElementLocalSystem() { }

    NearestElementLocalSystem(NodePointerType pNode) : mpNode(pNode)
    {

    }

    MapperLocalSystemUniquePointer Create(NodePointerType pNode) const override
    {
        return Kratos::make_unique<NearestElementLocalSystem>(pNode);
    }

    void CalculateAll(MappingWeightsVector& rMappingWeights,
                        EquationIdVectorType& rOriginIds,
                        EquationIdVectorType& rDestinationIds) const override;

    bool UseNodesAsBasis() const override { return true; }

private:
    NodePointerType mpNode;

};

/// Interpolative Mapper
/** This class implements the Nearest Element Mapping technique.
* Each node on the destination side gets assigned is's closest condition or element (distance to center)
* on the other side of the interface.
* In the mapping phase every node gets assigned the interpolated value of the condition/element.
* The interpolation is done with the shape funcitons
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
template<class TSparseSpace, class TDenseSpace>
class NearestElementMapper : public Mapper<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestElementMapper
    KRATOS_CLASS_POINTER_DEFINITION(NearestElementMapper);

    using BaseType = Mapper<TSparseSpace, TDenseSpace>;
    using MapperLocalSystemPointer = typename BaseType::MapperLocalSystemPointer;

    using MapperUniquePointerType = typename BaseType::MapperUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination)
                         : Mapper<TSparseSpace, TDenseSpace>(rModelPartOrigin,
                                  rModelPartDestination) {}

    NearestElementMapper(ModelPart& rModelPartOrigin,
                         ModelPart& rModelPartDestination,
                         Parameters JsonParameters)
                         : Mapper<TSparseSpace, TDenseSpace>(rModelPartOrigin,
                                  rModelPartDestination,
                                  JsonParameters)
    {
        // The Initialize function has to be called here bcs it internally calls virtual
        // functions that would not exist yet if it was called from the BaseClass!
        this->Initialize();


        // mpMapperCommunicator->InitializeOrigin(MapperUtilities::Condition_Center);
        // mpMapperCommunicator->InitializeDestination(MapperUtilities::Node_Coords);
        // mpMapperCommunicator->Initialize();

        // mpInverseMapper.reset(); // explicitly specified to be safe
    }

    /// Destructor.
    virtual ~NearestElementMapper() { }


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
        return Kratos::make_unique<NearestElementMapper<TSparseSpace, TDenseSpace>>(rModelPartOrigin,
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
        std::stringstream buffer;
        buffer << "NearestElementMapper" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NearestElementMapper";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


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

    MapperLocalSystemPointer GetMapperLocalSystem() const override
    {
        return Kratos::make_unique<NearestElementLocalSystem>();
    }

    InterfaceObject::ConstructionType GetInterfaceObjectConstructionTypeOrigin() const override
    {
        return InterfaceObject::Geometry_Center;
    }

    InterfaceObject::ConstructionType GetInterfaceObjectConstructionTypeDestination() const override
    {
        return InterfaceObject::Node_Coords;
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
    // NearestElementMapper& operator=(NearestElementMapper const& rOther) {}

    //   /// Copy constructor.
    //   NearestElementMapper(NearestElementMapper const& rOther){}


    ///@}

}; // Class NearestElementMapper

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEAREST_ELEMENT_MAPPER_H_INCLUDED  defined
