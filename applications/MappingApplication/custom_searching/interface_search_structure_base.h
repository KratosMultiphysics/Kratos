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

#if !defined(KRATOS_INTERFACE_SEARCH_STRUCTURE_BASE_H_INCLUDED )
#define  KRATOS_INTERFACE_SEARCH_STRUCTURE_BASE_H_INCLUDED

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

/// Rank local searching
/** This class provides features for rank local searching. It is basically a wrapper for the Bin Search Function
* It computes local neighbors and selects the best neighbor out of the search results. How these are selected is
* defined in the "EvaluateResult" function of the InterfaceObject
* If no neighbors are found, the search radius is increased by a factor ("increase_factor" in "Search")
* Look into the class description of the MapperCommunicator to see how this Object is used in the application
*/
class InterfaceSearchStructureBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceSearchStructureBase
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceSearchStructureBase);

    using MapperInterfaceInfoUniquePointerType = Kratos::unique_ptr<MapperInterfaceInfo>;

    using MapperInterfaceInfoPointerType = Kratos::shared_ptr<MapperInterfaceInfo>;
    using MapperInterfaceInfoPointerVectorType = std::vector<std::vector<MapperInterfaceInfoPointerType>>;
    using MapperInterfaceInfoPointerVectorPointerType = Kratos::unique_ptr<MapperInterfaceInfoPointerVectorType>;

    using MapperLocalSystemPointer = Kratos::unique_ptr<MapperLocalSystem>;
    using MapperLocalSystemPointerVector = std::vector<MapperLocalSystemPointer>;
    using MapperLocalSystemPointerVectorPointer = Kratos::shared_ptr<MapperLocalSystemPointerVector>;

    using BinsUniquePointerType = Kratos::unique_ptr<BinsObjectDynamic<InterfaceObjectConfigure>>;

    using InterfaceObjectContainerType = InterfaceObjectConfigure::ContainerType;
    using InterfaceObjectContainerUniquePointerType = Kratos::unique_ptr<InterfaceObjectContainerType>;

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceSearchStructureBase(ModelPart& rModelPartOrigin,
                             MapperLocalSystemPointerVectorPointer pMapperLocalSystems,
                             Parameters SearchSettings)
        : mrModelPartOrigin(rModelPartOrigin),
          mpMapperLocalSystems(pMapperLocalSystems),
          mSearchSettings(SearchSettings)
    {
        mEchoLevel = mSearchSettings["echo_level"].GetInt();
        mpMapperInterfaceInfosContainer = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>();
    }


    /// Destructor.
    virtual ~InterfaceSearchStructureBase() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void ExchangeInterfaceData(const Communicator& rComm,
                               const Kratos::Flags& rOptions,
                               const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                               InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "InterfaceSearchStructureBase" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceSearchStructureBase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}


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

    ModelPart& mrModelPartOrigin;
    MapperLocalSystemPointerVectorPointer mpMapperLocalSystems;
    MapperInterfaceInfoPointerVectorPointerType mpMapperInterfaceInfosContainer; // this contains the InterfaceInfos for all ranks! => needed to do the async communication

    BinsUniquePointerType mpLocalBinStructure;

    InterfaceObjectContainerUniquePointerType mpInterfaceObjectsOrigin;

    Parameters mSearchSettings;
    double mSearchRadius;

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void ConductLocalSearch();

    void CreateInterfaceObjectsOrigin(InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    void UpdateInterfaceObjectsOrigin();

    void InitializeBinsSearchStructure();

    virtual void PrepareSearch(const Kratos::Flags& rOptions,
                                        const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                        InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    virtual void FinalizeSearch();

    virtual void PrepareSearchIteration(const Kratos::Flags& rOptions,
                                        const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                        InterfaceObject::ConstructionType InterfaceObjectTypeOrigin) = 0;

    virtual void FinalizeSearchIteration() = 0;

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

    // int m_omp_threshold_num_nodes = 1000; // TODO constexpr???

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    // this function performs the search and the exchange of the data on the interface
    void ConductSearchIteration(const Kratos::Flags& rOptions,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    bool AllNeighborsFound(const Communicator& rComm) const;

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
    // InterfaceSearchStructureBase& operator=(InterfaceSearchStructureBase const& rOther){}

    //   /// Copy constructor.
    //   InterfaceSearchStructureBase(InterfaceSearchStructureBase const& rOther){}


    ///@}

}; // Class InterfaceSearchStructureBase

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_SEARCH_STRUCTURE_BASE_H_INCLUDED  defined
