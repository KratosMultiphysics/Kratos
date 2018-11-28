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

#if !defined(KRATOS_MAPPER_INTERFACE_INFO_H_INCLUDED)
#define  KRATOS_MAPPER_INTERFACE_INFO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_searching/interface_object.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Object for storing data that is needed to construct the local-mapping-system
/** When constructing th local-mapper-system, some data from the "other" side of the
 * Interface is needed which is different for every mapper.
 * This class stores the data needed and provides it later when the local-system is
 * constructed.
 * For the data-exchange in MPI the data is saved using serialization
*/
class MapperInterfaceInfo
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperInterfaceInfo
    KRATOS_CLASS_POINTER_DEFINITION(MapperInterfaceInfo);

    typedef std::size_t IndexType;

    typedef typename InterfaceObject::CoordinatesArrayType CoordinatesArrayType;

    typedef InterfaceObject::NodeType NodeType;
    typedef InterfaceObject::GeometryType GeometryType;

    ///@}
    ///@name  Enum's
    ///@{

    enum InfoType
    {
        Dummy
    };

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperInterfaceInfo() = default;

    MapperInterfaceInfo(const CoordinatesArrayType& rCoordinates,
                        const IndexType SourceLocalSystemIndex,
                        const IndexType SourceRank)
        : mSourceLocalSystemIndex(SourceLocalSystemIndex),
          mCoordinates(rCoordinates),
          mSourceRank(SourceRank)
    {}

    /// Destructor.
    virtual ~MapperInterfaceInfo() = default;

    ///@}
    ///@name Operations
    ///@{

    /**
    * @brief Processing the result of the search
    * This function processes the results of the search, i.e. it extracts
    * the information needed from the neighbor and saves it such that it can be
    * later accesse by the MapperLocalSystem for assembling it's local system.
    * This happens in the remote partition.
    * @param rpInterfaceObject The InterfaceObject found by the search
    * @param rpInterfaceObject The distance to the InterfaceObject found by the search
    * @author Philipp Bucher
    */
    virtual void ProcessSearchResult(const InterfaceObject::Pointer& rpInterfaceObject,
                                     const double NeighborDistance) = 0;

    /**
    * @brief Processing the result of the search for computing an approximation
    * This function processes the results of the search for computing an approximation.
    * This can be necessary if e.g. a projection fails.
    * In case an approximation is found the function "SetIsApproximation" has to be
    * called to set internal flags for data-exchange and the computation of the
    * local system by the MapperLocalSystem
    * This happens in the remote partition.
    * It's implementation is optional
    * @param rpInterfaceObject The InterfaceObject found by the search
    * @param rpInterfaceObject The distance to the InterfaceObject found by the search
    * @see ProcessSearchResult
    * @see SetIsApproximation
    * @author Philipp Bucher
    */
    virtual void ProcessSearchResultForApproximation(
        const InterfaceObject::Pointer& rpInterfaceObject,
        const double NeighborDistance) {}

    virtual MapperInterfaceInfo::Pointer Create(const CoordinatesArrayType& rCoordinates,
                                                const IndexType SourceLocalSystemIndex,
                                                const IndexType SourceRank=0) const = 0;

    // needed for serialization
    virtual MapperInterfaceInfo::Pointer Create() const = 0;

    IndexType GetLocalSystemIndex() const { return mSourceLocalSystemIndex; }

    IndexType GetSourceRank() const { return mSourceRank; }

    bool GetLocalSearchWasSuccessful() const { return mLocalSearchWasSuccessful; }

    bool GetIsApproximation() const { return mIsApproximation; }

    CoordinatesArrayType& Coordinates()
    {
        return mCoordinates;
    }

    ///@}
    ///@name Access
    ///@{

    virtual void GetValue(int& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::size_t& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(double& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(bool& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(GeometryType& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }

    virtual void GetValue(std::vector<int>& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<std::size_t>& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<double>& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<bool>& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }
    virtual void GetValue(std::vector<GeometryType>& rValue, const InfoType ValueType=InfoType::Dummy) const { KRATOS_ERROR << "Base class function called!" << std::endl; }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "MapperInterfaceInfo";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    // These variables need serialization
    IndexType mSourceLocalSystemIndex;

    // These variables are NOT being serialized bcs they are not needed after searching!
    CoordinatesArrayType mCoordinates;
    IndexType mSourceRank = 0;

    ///@}
    ///@name Protected Operations
    ///@{

    void SetLocalSearchWasSuccessful() { mLocalSearchWasSuccessful = true; }

    void SetIsApproximation()
    {
        // If an approximation is found also means that the local search has been successful!
        // this is needed otherwise it won't be properly processes by the search
        // the MapperLocalSystem has to take care of this!
        mLocalSearchWasSuccessful = true;

        mIsApproximation = true;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    bool mIsApproximation = false;

    bool mLocalSearchWasSuccessful = false; // this is not being serialized since it is not needed after mpi-data-exchange!

    ///@}
    ///@name Private  Access
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        rSerializer.save("LocalSysIdx", mSourceLocalSystemIndex);
        rSerializer.save("IsApproximation", mIsApproximation);
    }

    virtual void load(Serializer& rSerializer)
    {
        rSerializer.load("LocalSysIdx", mSourceLocalSystemIndex);
        rSerializer.load("IsApproximation", mIsApproximation);
    }

    ///@}

}; // Class MapperInterfaceInfo

///@}

///@name Type Definitions
///@{

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_INTERFACE_INFO_H_INCLUDED  defined
