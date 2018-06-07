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
    using MapperInterfaceInfoPointerVectorType = std::vector<MapperInterfaceInfoPointerType>;
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
        mpMapperInterfaceInfos = Kratos::make_unique<MapperInterfaceInfoPointerVectorType>();
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



    // void Search(const double SearchRadius, const int MaxSearchIterations)
    // {
    //     // mSearchRadius = SearchRadius;
    //     // mMaxSearchIterations = MaxSearchIterations;
    //     // const int increase_factor = 4;
    //     // int num_iteration = 1;
    //     // bool last_iteration = false;

    //     // if (mMaxSearchIterations == 1)   // in case only one search iteration is conducted
    //     // {
    //     //     last_iteration = true;
    //     // }

    //     // // First Iteration is done outside the search loop bcs it has
    //     // // to be done in any case
    //     // // one search iteration should be enough in most cases (if the search
    //     // // radius was either computed or specified properly)
    //     // // only if some points did not find a neighbor or dont have a valid
    //     // // projection, more search iterations are necessary
    //     // ConductSearchIteration(last_iteration);

    //     // while (num_iteration < mMaxSearchIterations && !mpInterfaceObjectManager->AllNeighborsFound())
    //     // {
    //     //     mSearchRadius *= increase_factor;
    //     //     ++num_iteration;

    //     //     if (num_iteration == mMaxSearchIterations)
    //     //     {
    //     //         last_iteration = true;
    //     //     }

    //     //     if (mEchoLevel >= 2 && mCommRank == 0)
    //     //     {
    //     //         std::cout << "MAPPER WARNING, search radius was increased, "
    //     //                   << "another search iteration is conducted, "
    //     //                   << "search iteration " << num_iteration << " / "
    //     //                   << mMaxSearchIterations << ", search radius "
    //     //                   << mSearchRadius << std::endl;
    //     //     }

    //     //     ConductSearchIteration(last_iteration);
    //     // }
    //     // if (mEchoLevel >= 2)
    //     // {
    //     //     mpInterfaceObjectManager->CheckResults();
    //     // }
    // }


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
    MapperInterfaceInfoPointerVectorPointerType mpMapperInterfaceInfos;

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

    // virtual void Initialize(InterfaceObject::ConstructionType InterfaceObjectTypeOrigin)
    // {
    //     CreateInterfaceObjectsOrigin(InterfaceObjectTypeOrigin);
    //     InitializeBinsSearchStructure();
    //     mInitializeIsPerformed = true;
    // }

    // This function constructs the InterfaceObjects on the Destination
    // In serial it only does it once, whereas in MPI this involves Data-Exchange!
    // Imagine a sliding interface, there the partitions might change!
    virtual void PrepareSearch(const Kratos::Flags& rOptions,
                                        const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                        InterfaceObject::ConstructionType InterfaceObjectTypeOrigin) = 0;

    virtual void FinalizeSearch() = 0;

    // This function constructs the InterfaceObjects on the Destination
    // In serial it only does it once, whereas in MPI this involves Data-Exchange!
    // Imagine a sliding interface, there the partitions might change!
    virtual void PrepareSearchIteration(const Kratos::Flags& rOptions,
                                        const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                        InterfaceObject::ConstructionType InterfaceObjectTypeOrigin) = 0;

    virtual void FinalizeSearchIteration() = 0;


    // void FindLocalNeighbors(InterfaceObjectConfigure::ContainerType& rInterfaceObjects,
    //                         const int InterfaceObjectsSize, std::vector<InterfaceObject::Pointer>& rInterfaceObjectResults,
    //                         std::vector<double>& rMinDistances, std::vector<std::vector<double>>& rShapeFunctionValues,
    //                         std::vector<int>& rPairingIndices)
    // {
    //     // This function finds neighbors of the InterfaceObjects in rInterfaceObjects in bin_structure
    //     // It must be executable by serial and parallel version!
    //     // InterfaceObjectsSize must be passed bcs rInterfaceObjects might contain old entries (it has
    //     // the max receive buffer size as size)!

    //     std::size_t num_inteface_obj_bin = mpInterfaceObjectsOrigin->size();

    //     if (num_inteface_obj_bin > 0)   // this partition has a bin structure
    //     {
    //         InterfaceObjectConfigure::ResultContainerType neighbor_results(num_inteface_obj_bin);
    //         std::vector<double> neighbor_distances(num_inteface_obj_bin);

    //         InterfaceObjectConfigure::IteratorType interface_object_itr;
    //         InterfaceObjectConfigure::ResultIteratorType results_itr;
    //         std::vector<double>::iterator distance_itr;

    //         //   Searching the neighbors
    //         for (int i = 0; i < InterfaceObjectsSize; ++i)
    //         {
    //             interface_object_itr = rInterfaceObjects.begin() + i;
    //             double search_radius = mSearchRadius; // reset search radius

    //             results_itr = neighbor_results.begin();
    //             distance_itr = neighbor_distances.begin();

    //             std::size_t number_of_results = mpLocalBinStructure->SearchObjectsInRadius(
    //                                                 *interface_object_itr, search_radius, results_itr,
    //                                                 distance_itr, num_inteface_obj_bin);

    //             if (number_of_results > 0)   // neighbors were found
    //             {
    //                 SelectBestResult(interface_object_itr, neighbor_results,
    //                                  neighbor_distances, number_of_results,
    //                                  rInterfaceObjectResults[i], rMinDistances[i],
    //                                  rShapeFunctionValues[i], rPairingIndices[i]);
    //             }
    //             else
    //             {
    //                 rMinDistances[i] = -1.0f; // indicates that the search was not succesful
    //                 rInterfaceObjectResults[i].reset(); // Release an old pointer, that is probably existing from a previous search
    //             }
    //         }
    //     }
    //     else     // this partition has no part of the point receiving interface, i.e. the origin of the mapped values
    //     {
    //         for (int i = 0; i < InterfaceObjectsSize; ++i)   // no results in this partition
    //         {
    //             rMinDistances[i] = -1.0f; // indicates that the search was not succesful
    //             rInterfaceObjectResults[i].reset(); // Release an old pointer, that is probably existing from a previous search
    //         }
    //     }
    // }

    // void SelectBestResult(const InterfaceObjectConfigure::IteratorType& rPoint,
    //                       const InterfaceObjectConfigure::ResultContainerType& rResultList,
    //                       const std::vector<double>& rDistances, const std::size_t NumResults,
    //                       InterfaceObject::Pointer& rVecClosestResults, double& rClosestDistance,
    //                       std::vector<double>& rShapeFunctionsValues, int& rPairingStatus)
    // {

    //     double min_distance = std::numeric_limits<double>::max();
    //     rClosestDistance = -1.0f; // indicate a failed search in case no result is good
    //     rPairingStatus = InterfaceObject::PairingStatus::NoNeighbor;

    //     for (int i = 0; i < static_cast<int>(NumResults); ++i)   // find index of best result
    //     {
    //         if (rResultList[i]->EvaluateResult((*rPoint)->Coordinates(), min_distance,
    //                                            rDistances[i], rShapeFunctionsValues))
    //         {
    //             rClosestDistance = min_distance;
    //             rVecClosestResults = rResultList[i];
    //             rPairingStatus = InterfaceObject::PairingStatus::NeighborFound;
    //         }
    //     }

    //     if (rPairingStatus != InterfaceObject::PairingStatus::NeighborFound)
    //     {
    //         for (int i = 0; i < static_cast<int>(NumResults); ++i)   // find index of best result
    //         {
    //             if (rResultList[i]->ComputeApproximation((*rPoint)->Coordinates(), min_distance,
    //                     rShapeFunctionsValues))
    //             {
    //                 rClosestDistance = min_distance;
    //                 rVecClosestResults = rResultList[i];
    //                 rPairingStatus = InterfaceObject::PairingStatus::Approximation;
    //             }
    //         }
    //     }
    // }


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

    // virtual void ConductSearchIteration(const bool LastIteration)
    // {
    //     // InterfaceObjectConfigure::ContainerType interface_objects;
    //     // mpInterfaceObjectManager->GetInterfaceObjectsSerialSearch(interface_objects);

    //     // int num_objects = interface_objects.size();

    //     // std::vector<InterfaceObject::Pointer> interface_object_results(num_objects);
    //     // std::vector<double> min_distances(num_objects);
    //     // std::vector<std::vector<double>> shape_function_values(num_objects);
    //     // std::vector<int> pairing_indices(num_objects);

    //     // FindLocalNeighbors(interface_objects, num_objects, interface_object_results,
    //     //                    min_distances, shape_function_values, pairing_indices);

    //     // mpInterfaceObjectManagerBins->StoreSearchResults(min_distances, interface_object_results, shape_function_values);
    //     // mpInterfaceObjectManager->PostProcessReceivedResults(interface_objects, min_distances,
    //     //         pairing_indices);
    // }



    // this function performs the search and the exchange of the data on the interface
    void ConductSearchIteration(const Kratos::Flags& rOptions,
                                const MapperInterfaceInfoUniquePointerType& rpRefInterfaceInfo,
                                InterfaceObject::ConstructionType InterfaceObjectTypeOrigin);

    bool AllNeighborsFound(const Communicator& rComm) const;


    // void InitializeInterfaceNodeManager(ModelPart& rModelPart)
    // {
    //     mInterfaceObjects.resize(rModelPart.GetCommunicator().LocalMesh().NumberOfNodes());

    //     int i = 0;
    //     for (auto &local_node : rModelPart.GetCommunicator().LocalMesh().Nodes())
    //     {
    //         mInterfaceObjects[i] = InterfaceObject::Pointer( new InterfaceNode(local_node, mEchoLevel) );
    //         ++i;
    //     }
    // }

    // void InitializeInterfaceGeometryObjectManager(ModelPart& rModelPart,
    //         GeometryData::IntegrationMethod IntegrationMethod,
    //         const double ApproximationTolerance)
    // {
    //     bool construct_with_center;
    //     int size_factor = 1;
    //     if (IntegrationMethod == GeometryData::NumberOfIntegrationMethods)
    //     {
    //         construct_with_center = true;
    //         size_factor = 1;
    //     }
    //     else
    //     {
    //         construct_with_center = false;
    //         if (rModelPart.GetCommunicator().LocalMesh().NumberOfConditions() > 0)
    //         {
    //             size_factor = rModelPart.GetCommunicator().LocalMesh().ConditionsBegin()->GetGeometry().IntegrationPointsNumber(IntegrationMethod);
    //         }
    //         else if (rModelPart.GetCommunicator().LocalMesh().NumberOfElements() > 0)
    //         {
    //             size_factor = rModelPart.GetCommunicator().LocalMesh().ElementsBegin()->GetGeometry().IntegrationPointsNumber(IntegrationMethod);
    //         }
    //     }

    //     mInterfaceObjects.reserve(size_factor * rModelPart.GetCommunicator().LocalMesh().NumberOfConditions());

    //     if (construct_with_center)   // construct with condition center point
    //     {
    //         for (auto& condition : rModelPart.GetCommunicator().LocalMesh().Conditions())
    //         {
    //             mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(condition.GetGeometry(),
    //                                         ApproximationTolerance,
    //                                         mEchoLevel,
    //                                         0) ));
    //         }
    //         for (auto& element : rModelPart.GetCommunicator().LocalMesh().Elements())
    //         {
    //             mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(element.GetGeometry(),
    //                                         ApproximationTolerance,
    //                                         mEchoLevel,
    //                                         0) ));
    //         }
    //     }
    //     else     // construct with condition gauss points
    //     {
    //         KRATOS_ERROR << "This is not implemented at the moment" << std::endl;
    //         for (auto& condition : rModelPart.GetCommunicator().LocalMesh().Conditions())
    //         {
    //             for (int g = 0; g < 111111; ++g) // TODO fix this, should be number of GPs
    //             {
    //                 mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(condition.GetGeometry(),
    //                                             ApproximationTolerance,
    //                                             mEchoLevel,
    //                                             g, IntegrationMethod) ));
    //             }
    //         }
    //         for (auto& element : rModelPart.GetCommunicator().LocalMesh().Elements())
    //         {
    //             for (int g = 0; g < 111111; ++g) // TODO fix this, should be number of GPs
    //             {
    //                 mInterfaceObjects.push_back(InterfaceObject::Pointer( new InterfaceGeometryObject(element.GetGeometry(),
    //                                             ApproximationTolerance,
    //                                             mEchoLevel,
    //                                             g, IntegrationMethod) ));
    //             }
    //         }

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
