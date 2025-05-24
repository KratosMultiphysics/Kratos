//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Jordi Cotela
//                   Carlos Roig
//

#pragma once

// Project includes
#include "includes/model_part_io.h"
#include "custom_processes/metis_divide_heterogeneous_input_process.h"


namespace Kratos
{
///@addtogroup MetisApplication
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

/// Call Metis to divide an heterogeneous mesh, by partitioning its nodal graph.
class KRATOS_API(METIS_APPLICATION) MetisDivideHeterogeneousInputInMemoryProcess : public MetisDivideHeterogeneousInputProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MetisDivideHeterogeneousInputInMemoryProcess
    KRATOS_CLASS_POINTER_DEFINITION(MetisDivideHeterogeneousInputInMemoryProcess);

    typedef MetisDivideHeterogeneousInputProcess BaseType;

    using BaseType::SizeType;
    using BaseType::GraphType;
    using BaseType::idxtype;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MetisDivideHeterogeneousInputInMemoryProcess(IO& rIO, ModelPartIO& rSerialIO, const DataCommunicator& rDataComm, int Dimension = 3, int Verbosity = 0, bool SynchronizeConditions = false):
        BaseType(rIO,rDataComm.Size(),Dimension,Verbosity,SynchronizeConditions), mrSerialIO(rSerialIO), mrDataComm(rDataComm)
    {
        KRATOS_ERROR_IF_NOT(mrDataComm.IsDistributed()) << "DataCommunicator must be distributed!" << std::endl;
    }

    /// Destructor.
    virtual ~MetisDivideHeterogeneousInputInMemoryProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        this->Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    /// Generate a partition using Metis.
    /** Partitioned input is written as <problem name>_<mpi rank>.mdpa
     */
    void Execute() override;

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
        return "MetisDivideHeterogeneousInputInMemoryProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "MetisDivideHeterogeneousInputInMemoryProcess";
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

    ModelPartIO& mrSerialIO;

    const DataCommunicator& mrDataComm;

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

    // Copy constructor.
    //MetisDivideHeterogeneousInputInMemoryProcess(MetisDivideHeterogeneousInputInMemoryProcess const& rOther);


    ///@}

}; // Class MetisDivideHeterogeneousInputInMemoryProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} // addtogroup block

}
