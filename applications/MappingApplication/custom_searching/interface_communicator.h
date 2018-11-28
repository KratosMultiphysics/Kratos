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

#if !defined(KRATOS_INTERFACE_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_INTERFACE_COMMUNICATOR_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/communicator.h"
#include "includes/kratos_parameters.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "custom_searching/custom_configures/interface_object_configure.h"
#include "custom_utilities/mapper_local_system.h"
#include "custom_utilities/mapper_flags.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Object for exchanging data on the Interface
/** Mapping requires knowledge about the "other" side of the Interface. This class communicates the
 * data required by the mappers, hence it also includes the (local) searching
*/
class InterfaceCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceCommunicator);

    typedef Kratos::unique_ptr<MapperInterfaceInfo> MapperInterfaceInfoUniquePointerType;

    typedef Kratos::shared_ptr<MapperInterfaceInfo> MapperInterfaceInfoPointerType;
    typedef std::vector<std::vector<MapperInterfaceInfoPointerType>> MapperInterfaceInfoPointerVectorType;
    typedef Kratos::unique_ptr<MapperInterfaceInfoPointerVectorType> MapperInterfaceInfoPointerVectorPointerType;

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;
    typedef Kratos::shared_ptr<MapperLocalSystemPointerVector> MapperLocalSystemPointerVectorPointer;

    typedef Kratos::unique_ptr<BinsObjectDynamic<InterfaceObjectConfigure>> BinsUniquePointerType;

    typedef InterfaceObjectConfigure::ContainerType InterfaceObjectContainerType;
    typedef Kratos::unique_ptr<InterfaceObjectContainerType> InterfaceObjectContainerUniquePointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceCommunicator(ModelPart& rModelPartOrigin,
                             MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                             Parameters SearchSettings)
        : mrModelPartOrigin(rModelPartOrigin),
          mpMapperLocalSystems(pMapperLocalSystems),
          mSearchSettings(SearchSettings)
    {
        // mEchoLevel = mSearchSettings["echo_level"].GetInt();
        mpMapperInterfaceInfosContainer = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>();
        mpMapperInterfaceInfosContainer->resize(1);
    }

    /// Destructor.
    virtual ~InterfaceCommunicator() = default;

    ///@}
    ///@name Operations
    ///@{

    void ExchangeInterfaceData(const Communicator& rComm,
                               const Kratos::Flags& rOptions,
                               const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                               InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "InterfaceCommunicator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceCommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPartOrigin;
    MapperLocalSystemPointerVectorPointer mpMapperLocalSystems;
    MapperInterfaceInfoPointerVectorPointerType mpMapperInterfaceInfosContainer; // this contains the InterfaceInfos for all ranks! => needed to do the async communication

    BinsUniquePointerType mpLocalBinStructure;

    InterfaceObjectContainerUniquePointerType mpInterfaceObjectsOrigin;

    Parameters mSearchSettings;
    double mSearchRadius;

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operations
    ///@{

    virtual void InitializeSearch(const Kratos::Flags& rOptions,
                                        const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                        InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    virtual void FinalizeSearch();

    virtual void InitializeSearchIteration(const Kratos::Flags& rOptions,
                                           const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    virtual void FinalizeSearchIteration(const MapperInterfaceInfoUniquePointerType& rpInterfaceInfo);

    void FilterInterfaceInfosSuccessfulSearch();

    void AssignInterfaceInfos();

    ///@}

private:
    ///@name Private Operations
    ///@{

    void ConductLocalSearch();

    void CreateInterfaceObjectsOrigin(InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    void UpdateInterfaceObjectsOrigin();

    void InitializeBinsSearchStructure();

    // this function performs the search and the exchange of the data on the interface
    void ConductSearchIteration(const Kratos::Flags& rOptions,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo);

    bool AllNeighborsFound(const Communicator& rComm) const;

    ///@}

}; // Class InterfaceCommunicator

///@}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_COMMUNICATOR_H_INCLUDED  defined
