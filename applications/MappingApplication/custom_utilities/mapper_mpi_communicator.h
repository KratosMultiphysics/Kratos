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

#if !defined(KRATOS_MAPPER_MPI_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_MAPPER_MPI_COMMUNICATOR_H_INCLUDED

// System includes

// External includes

// Project includes
#include "mapper_communicator.h"
#include "interface_search_structure_mpi.h"
#include "interface_object_manager_parallel.h"
#include "mapper_utilities_mpi.h"


namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

typedef matrix<int> GraphType; // GraphColoringProcess

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
class MapperMPICommunicator : public MapperCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperMPICommunicator
    KRATOS_CLASS_POINTER_DEFINITION(MapperMPICommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    MapperMPICommunicator(ModelPart& rModelPartOrigin, ModelPart& rModelPartDestination,
                          Parameters& rJsonParameters) :
        MapperCommunicator(rModelPartOrigin, rModelPartDestination, rJsonParameters) { }

    /// Destructor.
    virtual ~MapperMPICommunicator() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void InitializeOrigin(MapperUtilities::InterfaceObjectConstructionType InterfaceObjectTypeOrigin,
                          GeometryData::IntegrationMethod IntegrationMethodOrigin = GeometryData::NumberOfIntegrationMethods) override
    {

        mpInterfaceObjectManagerOrigin = InterfaceObjectManagerParallel::Pointer(
                                             new InterfaceObjectManagerParallel(mrModelPartOrigin, MyPID(),
                                                     TotalProcesses(),
                                                     InterfaceObjectTypeOrigin,
                                                     IntegrationMethodOrigin,
                                                     mEchoLevel,
                                                     mApproximationTolerance) );

        if (mEchoLevel > 3)
        {
            mpInterfaceObjectManagerOrigin->PrintInterfaceObjects("Origin");
        }

        mInterfaceObjectTypeOrigin = InterfaceObjectTypeOrigin;
        mIntegrationMethodOrigin = IntegrationMethodOrigin;
    }

    void InitializeDestination(MapperUtilities::InterfaceObjectConstructionType InterfaceObjectTypeDestination,
                               GeometryData::IntegrationMethod IntegrationMethodDestination = GeometryData::NumberOfIntegrationMethods) override
    {

        mpInterfaceObjectManagerDestination = InterfaceObjectManagerParallel::Pointer(
                new InterfaceObjectManagerParallel(mrModelPartDestination, MyPID(),
                        TotalProcesses(),
                        InterfaceObjectTypeDestination,
                        IntegrationMethodDestination,
                        mEchoLevel,
                        mApproximationTolerance) );

        if (mEchoLevel > 3)
        {
            mpInterfaceObjectManagerDestination->PrintInterfaceObjects("Destination");
        }

        mInterfaceObjectTypeDestination = InterfaceObjectTypeDestination;
        mIntegrationMethodDestination = IntegrationMethodDestination;
    }

    void TransferVariableData(std::function<double(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                              std::function<void(InterfaceObject*, double)> FunctionPointerDestination,
                              const Variable<double>& rOriginVariable) override
    {
        mrModelPartOrigin.GetCommunicator().SynchronizeVariable(rOriginVariable); // required bcs
        // data interpolation can also involve ghost nodes
        ExchangeDataLocal(FunctionPointerOrigin, FunctionPointerDestination);
        ExchangeDataRemote(FunctionPointerOrigin, FunctionPointerDestination);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    void TransferVariableData(std::function<array_1d<double, 3>(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                              std::function<void(InterfaceObject*, array_1d<double, 3>)> FunctionPointerDestination,
                              const Variable< array_1d<double, 3> >& rOriginVariable) override
    {
        mrModelPartOrigin.GetCommunicator().SynchronizeVariable(rOriginVariable); // required bcs
        // data interpolation can also involve ghost nodes
        ExchangeDataLocal(FunctionPointerOrigin, FunctionPointerDestination);
        ExchangeDataRemote(FunctionPointerOrigin, FunctionPointerDestination);

        MPI_Barrier(MPI_COMM_WORLD);
    }

    int MyPID () override   // Copy from "kratos/includes/mpi_communicator.h"
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return rank;
    }

    int TotalProcesses() override   // Copy from "kratos/includes/mpi_communicator.h"
    {
        int nproc;
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        return nproc;
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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "MapperMPICommunicator" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MapperMPICommunicator";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}


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

    int mMaxSendBufferSize;
    int mMaxReceiveBufferSize;

    GraphType mColoredGraph;
    int mMaxColors;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void InitializeSearchStructure() override
    {
        mpSearchStructure = InterfaceSearchStructure::Pointer ( new InterfaceSearchStructureMPI(
                                mpInterfaceObjectManagerDestination, mpInterfaceObjectManagerOrigin, mrModelPartOrigin,
                                MyPID(), TotalProcesses(), mEchoLevel) );
    }

    void InvokeSearch(const double InitialSearchRadius,
                      const int MaxSearchIterations) override
    {
        mpSearchStructure->Search(InitialSearchRadius,
                                  MaxSearchIterations);

        mpInterfaceObjectManagerDestination->ComputeBufferSizesAndCommunicationGraph(
            mMaxSendBufferSize,
            mMaxReceiveBufferSize,
            mColoredGraph,
            mMaxColors);
        if (mEchoLevel > 3)
        {
            PrintPairs();
        }
    }

    template <typename T>
    void ExchangeDataRemote(std::function<T(InterfaceObject*, const std::vector<double>&)> FunctionPointerOrigin,
                            std::function<void(InterfaceObject*, T)> FunctionPointerDestination)
    {
        int send_buffer_size = 0;
        int receive_buffer_size = 0;

        int max_send_buffer_size = mMaxSendBufferSize;
        int max_receive_buffer_size = mMaxReceiveBufferSize;

        int buffer_size_factor = MapperUtilitiesMPI::SizeOfVariable(T());

        max_send_buffer_size *= buffer_size_factor;
        max_receive_buffer_size *= buffer_size_factor;

        double* send_buffer = new double[max_send_buffer_size];
        double* receive_buffer = new double[max_receive_buffer_size];

        for (int i = 0; i < mMaxColors; ++i)   // loop over communication steps (aka. colour)
        {
            int comm_partner = mColoredGraph(MyPID(), i); // get the partner rank
            if (comm_partner != -1)   // check if rank is communicating in this communication step (aka. colour)
            {
                mpInterfaceObjectManagerOrigin->FillBufferWithValues(send_buffer, send_buffer_size,
                        comm_partner, FunctionPointerOrigin);

                MapperUtilitiesMPI::MpiSendRecv(send_buffer, receive_buffer, send_buffer_size, receive_buffer_size,
                                                max_send_buffer_size, max_receive_buffer_size, comm_partner);

                mpInterfaceObjectManagerDestination->ProcessValues(receive_buffer, receive_buffer_size,
                        comm_partner, FunctionPointerDestination);

            } // if I am communicating in this loop (comm_partner != -1)
        } // loop colors

        delete [] send_buffer;
        delete [] receive_buffer;
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
    MapperMPICommunicator& operator=(MapperMPICommunicator const& rOther);

    //   /// Copy constructor.
    //   MapperMPICommunicator(MapperMPICommunicator const& rOther){}


    ///@}

}; // Class MapperMPICommunicator

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MapperMPICommunicator& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MapperMPICommunicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_MPI_COMMUNICATOR_H_INCLUDED  defined