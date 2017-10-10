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

#if !defined(KRATOS_MAPPER_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_MAPPER_COMMUNICATOR_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "interface_object_manager_serial.h"
#include "interface_search_structure.h"
#include "mapper_utilities.h"
#include "mapper_flags.h"


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

/// Interface btw the Mappers and the Core of the MappingApplication
/** This class is the top instance of the Core of the MappingApplication
* It handles the searching and the exchange of Data btw the Interfaces
* It also checks and validates the JSON input (default-parameters: 
* memeber variable "mDefaultParameters").
*
* Available Echo Levels:
* 0 : Mute every output
* == 1 : Print Timing Information (Mapper Construction and the three basic functions)
* >= 2 : Warnings are printed
* >= 3 : Basic Information, recommended for standard debugging
* >= 4 : Very detailed output for every object on the interface!
*       (Only recommended for debugging of small example, otherwise it gets very messy!
*       Should be only needed for developing)
*
* Structure of the how the different classes work with each other:
* (Illustrated in Thesis mentioned in the header)
* It is described in the serial case, if executed mpi-parallel, the MPI-versions
* of the classes are used: 
* MapperCommunicator
*   Constructs the InterfaceObjectManagers both for the Destination and the Origin
*       Which InterfaceObjects should be used is declared in the input of the "Initialize*" functions
*   It then constructs the SearchStructure and initializes the Search
*   For the exchange of values it implements the corresponding functions
*/
class MapperCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(MapperCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    MapperCommunicator(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
                       Parameters JsonParameters, const int CommRank = 0) :
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mJsonParameters(JsonParameters)
    {
        CheckAndValidateJson(CommRank);

        mInitialSearchRadius = mJsonParameters["search_radius"].GetDouble();
        mMaxSearchIterations = mJsonParameters["search_iterations"].GetInt();
        mApproximationTolerance = mJsonParameters["approximation_tolerance"].GetDouble();

        mEchoLevel = mJsonParameters["echo_level"].GetInt();
        
        CheckInterfaceModelParts(CommRank);
    }

    /// Destructor.
    virtual ~MapperCommunicator() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void InitializeOrigin(MapperUtilities::InterfaceObjectConstructionType InterfaceObjectTypeOrigin,
                                  GeometryData::IntegrationMethod IntegrationMethodOrigin = GeometryData::NumberOfIntegrationMethods)
    {

        mpInterfaceObjectManagerOrigin = InterfaceObjectManagerBase::Pointer(
                                             new InterfaceObjectManagerSerial(mrModelPartOrigin, MyPID(),
                                                     TotalProcesses(),
                                                     InterfaceObjectTypeOrigin,
                                                     IntegrationMethodOrigin,
                                                     mEchoLevel,
                                                     mApproximationTolerance) );

        if (mEchoLevel >= 4)
        {
            mpInterfaceObjectManagerOrigin->PrintInterfaceObjects("Origin");
        }

        // Save for updating the interface
        mInterfaceObjectTypeOrigin = InterfaceObjectTypeOrigin;
        mIntegrationMethodOrigin = IntegrationMethodOrigin;
    }

    virtual void InitializeDestination(MapperUtilities::InterfaceObjectConstructionType InterfaceObjectTypeDestination,
                                       GeometryData::IntegrationMethod IntegrationMethodDestination = GeometryData::NumberOfIntegrationMethods)
    {

        mpInterfaceObjectManagerDestination = InterfaceObjectManagerBase::Pointer(
                new InterfaceObjectManagerSerial(mrModelPartDestination, MyPID(),
                        TotalProcesses(),
                        InterfaceObjectTypeDestination,
                        IntegrationMethodDestination,
                        mEchoLevel,
                        mApproximationTolerance) );

        if (mEchoLevel >= 4)
        {
            mpInterfaceObjectManagerDestination->PrintInterfaceObjects("Destination");
        }

        // Save for updating the interface
        mInterfaceObjectTypeDestination = InterfaceObjectTypeDestination;
        mIntegrationMethodDestination = IntegrationMethodDestination;
    }

    void Initialize()
    {
        InitializeSearchStructure();
        InvokeSearch(mInitialSearchRadius, mMaxSearchIterations);
    }

    void UpdateInterface(Kratos::Flags& rOptions, double InitialSearchRadius)
    {
        if (rOptions.Is(MapperFlags::REMESHED))   // recompute the managers and the search structure
        {
            InitializeOrigin(mInterfaceObjectTypeOrigin, mIntegrationMethodOrigin);
            InitializeDestination(mInterfaceObjectTypeDestination, mIntegrationMethodDestination);
            InitializeSearchStructure();
        }
        else     // clear the managers
        {
            // TODO does the same for now, since the InterfaceObjects do not use the refs to their
            // original entities, so their position is not updated!
            InitializeOrigin(mInterfaceObjectTypeOrigin, mIntegrationMethodOrigin);
            InitializeDestination(mInterfaceObjectTypeDestination, mIntegrationMethodDestination);
            InitializeSearchStructure();
            // mpInterfaceObjectManagerOrigin->Clear();
            // mpInterfaceObjectManagerDestination->Clear();
            // InitializeSearchStructure();
        }

        if (InitialSearchRadius < 0.0f)
        {
            InitialSearchRadius = MapperUtilities::ComputeSearchRadius(mrModelPartOrigin,
                                  mrModelPartDestination);
        }
        mInitialSearchRadius = InitialSearchRadius; // update the search radius

        InvokeSearch(mInitialSearchRadius, mMaxSearchIterations);
    }

    virtual void TransferVariableData(std::function<double(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                                      std::function<void(InterfaceObject*, double)> FunctionPointerDestination,
                                      const Variable<double>& rOriginVariable)
    {
        ExchangeDataLocal(FunctionPointerOrigin, FunctionPointerDestination);
    }

    virtual void TransferVariableData(std::function<array_1d<double, 3>(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                                      std::function<void(InterfaceObject*, array_1d<double, 3>)> FunctionPointerDestination,
                                      const Variable< array_1d<double, 3> >& rOriginVariable)
    {
        ExchangeDataLocal(FunctionPointerOrigin, FunctionPointerDestination);
    }

    virtual void TransferVariableData(std::function<void*(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                                      std::function<void(InterfaceObject*, void*)> FunctionPointerDestination )
    {
        ExchangeDataLocal(FunctionPointerOrigin, FunctionPointerDestination);
    }
   

    virtual int MyPID() // Copy from "kratos/includes/communicator.h"
    {
        return 0;
    }

    virtual int TotalProcesses() // Copy from "kratos/includes/communicator.h"
    {
        return 1;
    }

    void PrintTime(const std::string& rMapperName,
                   const std::string& rFunctionName,
                   const double& rElapsedTime)
    {
        if (mEchoLevel == 1 && MyPID() == 0)
        {
            std::cout  << "MAPPER TIMER: \"" << rMapperName << "\", \"" << rFunctionName
                       << "\" took " <<  rElapsedTime << " seconds" << std::endl;
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
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "MapperCommunicator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperCommunicator";
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
    ModelPart& mrModelPartDestination;

    InterfaceObjectManagerBase::Pointer mpInterfaceObjectManagerOrigin;
    InterfaceObjectManagerBase::Pointer mpInterfaceObjectManagerDestination;

    MapperUtilities::InterfaceObjectConstructionType mInterfaceObjectTypeOrigin;
    GeometryData::IntegrationMethod mIntegrationMethodOrigin;
    MapperUtilities::InterfaceObjectConstructionType mInterfaceObjectTypeDestination;
    GeometryData::IntegrationMethod mIntegrationMethodDestination;

    InterfaceSearchStructure::Pointer mpSearchStructure;

    double mInitialSearchRadius;
    int mMaxSearchIterations;
    double mApproximationTolerance;

    int mEchoLevel = 0;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    template <typename T>
    void ExchangeDataLocal(std::function<T(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                           std::function<void(InterfaceObject*, T)> FunctionPointerDestination)
    {
        std::vector< T > values;

        mpInterfaceObjectManagerOrigin->FillBufferWithValues(FunctionPointerOrigin, values);
        mpInterfaceObjectManagerDestination->ProcessValues(FunctionPointerDestination, values);
    }

    // Function used for Debugging
    void PrintPairs()
    {
        KRATOS_ERROR << "Needs to be re-implemented!" << std::endl;
        // mpInterfaceObjectManagerOrigin->WriteNeighborRankAndCoordinates();
        // Kratos::Flags options = Kratos::Flags();
        // options.Set(MapperFlags::NON_HISTORICAL_DATA);
        // TransferNodalData(NEIGHBOR_RANK, NEIGHBOR_RANK, options);
        // TransferNodalData(NEIGHBOR_COORDINATES, NEIGHBOR_COORDINATES, options);
        // mpInterfaceObjectManagerDestination->PrintNeighbors();
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

    Parameters mJsonParameters;
    Parameters mDefaultParameters = Parameters( R"(
      {
             "mapper_type"                           : "",
             "interface_submodel_part_origin"        : "",
             "interface_submodel_part_destination"   : "",
             "search_radius"                         : -1.0,
             "search_iterations"                     : 3,
             "approximation_tolerance"               : -1.0,
             "echo_level"                            : 0
       }  )" );

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    virtual void InitializeSearchStructure()
    {
        mpSearchStructure = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructure(
                                mpInterfaceObjectManagerDestination, mpInterfaceObjectManagerOrigin, mEchoLevel) );
    }

    virtual void InvokeSearch(const double InitialSearchRadius,
                              const int MaxSearchIterations)
    {
        mpSearchStructure->Search(InitialSearchRadius,
                                  MaxSearchIterations);
        if (mEchoLevel >= 4)
        {
            // PrintPairs(); // TODO reimplement!
        }
    }

    // CommRank is used as input bcs the MyPID function of the non-MPI MapperCommunicator is used
    // since this function is called before the MapperMPICommunicato is initialized
    void CheckAndValidateJson(const int CommRank)
    {
        // Check if there is a valid input for the search parameters
        bool compute_search_radius = true;
        if (mJsonParameters.Has("search_radius"))
        {
            compute_search_radius = false;

            if (mJsonParameters["search_radius"].GetDouble() < 0.0f)
            {
                KRATOS_ERROR << "Invalid Search Radius specified" << std::endl;
            }
        }

        if (mJsonParameters.Has("search_iterations"))
        {
            if (mJsonParameters["search_iterations"].GetInt() < 1)
            {
                KRATOS_ERROR << "Number of specified Search Iterations too small" << std::endl;
            }
        }

        if (mJsonParameters.Has("approximation_tolerance"))
        {
            if (mJsonParameters["approximation_tolerance"].GetDouble() < 0.0f)
            {
                KRATOS_ERROR << "Invalid Tolerance for Approximations specified" << std::endl;
            }
        }

        if (mJsonParameters.Has("echo_level"))
        {
            if (mJsonParameters["echo_level"].GetInt() < 0)
            {
                KRATOS_ERROR << "Echo Level cannot be smaller than 0" << std::endl;
            }

            if (mJsonParameters["echo_level"].GetInt() >= 3 && CommRank == 0)
            {
                std::cout << "Mapper JSON Parameters BEFORE validation:\n "
                          << mJsonParameters.PrettyPrintJsonString() << std::endl;
            }
        }

        mJsonParameters.RecursivelyValidateAndAssignDefaults(mDefaultParameters);

        if (mJsonParameters["echo_level"].GetInt() >= 3 && CommRank == 0)
        {
            std::cout << "Mapper JSON Parameters AFTER validation:\n "
                      << mJsonParameters.PrettyPrintJsonString() << std::endl;
        }

        // Compute the search radius in case it was not specified
        if (compute_search_radius)
        {
            double search_radius = MapperUtilities::ComputeSearchRadius(mrModelPartOrigin,
                                   mrModelPartDestination,
                                   mJsonParameters["echo_level"].GetInt());
            mJsonParameters["search_radius"].SetDouble(search_radius);

            if (mJsonParameters["echo_level"].GetInt() >= 3 && CommRank == 0)
            {
                std::cout << "SearchRadius computed by MapperCommunicator = " << search_radius << std::endl;
            }
        }

        if (mJsonParameters["approximation_tolerance"].GetDouble() < 0.0f)   // nothing specified, set to max
        {
            mJsonParameters["approximation_tolerance"].SetDouble(std::numeric_limits<double>::max());
        }
    }

    // CommRank is used as input bcs the MyPID function of the non-MPI MapperCommunicator is used
    // since this function is called before the MapperMPICommunicato is initialized
    void CheckInterfaceModelParts(const int CommRank)
    {
        const int num_nodes_origin = MapperUtilities::ComputeNumberOfNodes(mrModelPartOrigin);
        const int num_conditions_origin = MapperUtilities::ComputeNumberOfConditions(mrModelPartOrigin);
        const int num_elements_origin = MapperUtilities::ComputeNumberOfElements(mrModelPartOrigin);

        const int num_nodes_destination = MapperUtilities::ComputeNumberOfNodes(mrModelPartDestination);
        const int num_conditions_destination = MapperUtilities::ComputeNumberOfConditions(mrModelPartDestination);
        const int num_elements_destination = MapperUtilities::ComputeNumberOfElements(mrModelPartDestination);

        // Check if the ModelPart contains entities
        if (num_nodes_origin + num_conditions_origin + num_elements_origin < 1)
        {
            KRATOS_ERROR << "Neither Nodes nor Conditions nor Elements found "
                         << "in the Origin ModelPart" << std::endl;
        }

        if (num_nodes_destination + num_conditions_destination + num_elements_destination < 1)
        {
            KRATOS_ERROR << "Neither Nodes nor Conditions nor Elements found "
                         << "in the Destination ModelPart" << std::endl;
        }

        // Check if the inpt ModelParts contain both Elements and Conditions
        // This is NOT possible, bcs the InterfaceObjects are constructed
        // with whatever exists in the Modelpart (see the InterfaceObjectManagerBase, 
        // function "InitializeInterfaceGeometryObjectManager")
        if (num_conditions_origin > 0 && num_elements_origin > 0)
        {
            KRATOS_ERROR << "Origin ModelPart contains both Conditions and Elements "
                         << "which is not permitted" << std::endl;
        }

        if (num_conditions_destination > 0 && num_elements_destination > 0)
        {
            KRATOS_ERROR << "Destination ModelPart contains both Conditions and Elements "
                         << "which is not permitted" << std::endl;
        }

        if (mEchoLevel >= 2) {
            std::vector<double> model_part_origin_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartOrigin);
            std::vector<double> model_part_destination_bbox = MapperUtilities::ComputeModelPartBoundingBox(mrModelPartDestination);

            bool bbox_overlapping = MapperUtilities::ComputeBoundingBoxIntersection(
                                                        model_part_origin_bbox,
                                                        model_part_destination_bbox);
            if(CommRank == 0)
            {
                if (!bbox_overlapping) {
                    std::cout << "MAPPER WARNING, the bounding boxes of the "
                              << "Modelparts do not overlap! "
                              << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
                                                                              model_part_destination_bbox)
                              << std::endl;
                } else if (mEchoLevel >= 3)
                {
                    std::cout << MapperUtilities::PrintModelPartBoundingBoxes(model_part_origin_bbox,
                                                                              model_part_destination_bbox)
                              << std::endl;
                }
            }
        }
    }

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
    MapperCommunicator& operator=(MapperCommunicator const& rOther);

    //   /// Copy constructor.
    //   MapperCommunicator(MapperCommunicator const& rOther){}


    ///@}

}; // Class MapperCommunicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MapperCommunicator& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MapperCommunicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_COMMUNICATOR_H_INCLUDED  defined
