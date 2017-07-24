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

/// Short class definition.
/** Detail class definition.
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
                       Parameters& rJsonParameters) :
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mrJsonParameters(rJsonParameters)
    {

        mInitialSearchRadius = mrJsonParameters["search_radius"].GetDouble();
        mMaxSearchIterations = mrJsonParameters["search_iterations"].GetInt();
        mApproximationTolerance = mrJsonParameters["approximation_tolerance"].GetDouble();

        mEchoLevel = mrJsonParameters["echo_level"].GetInt();
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

    Parameters& mrJsonParameters;
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

    bool mComputeSearchRadius;

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
