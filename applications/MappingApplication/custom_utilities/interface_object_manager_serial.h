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

#if !defined(KRATOS_INTERFACE_OBJECT_MANAGER_SERIAL_H_INCLUDED )
#define  KRATOS_INTERFACE_OBJECT_MANAGER_SERIAL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "interface_object_manager_base.h"


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
class InterfaceObjectManagerSerial : public InterfaceObjectManagerBase
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfaceObjectManagerSerial
    KRATOS_CLASS_POINTER_DEFINITION(InterfaceObjectManagerSerial);

    ///@}
    ///@name Life Cycle
    ///@{

    InterfaceObjectManagerSerial(ModelPart& rModelPart, int CommRank, int CommSize,
                                 MapperUtilities::InterfaceObjectConstructionType InterfaceObjectType,
                                 GeometryData::IntegrationMethod IntegrationMethod, const int EchoLevel,
                                 const double ApproximationTolerance) :
        InterfaceObjectManagerBase(
            rModelPart, CommRank, CommSize, InterfaceObjectType,
            IntegrationMethod, EchoLevel, ApproximationTolerance) { }

    /// Destructor.
    virtual ~InterfaceObjectManagerSerial() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // **********************************************************************
    // Side we want to find neighbors for aka destination *******************
    // **********************************************************************
    void GetInterfaceObjectsSerialSearch(InterfaceObjectConfigure::ContainerType& rCandidateSendObjects) override
    {
        InitializeSizes();
        for (auto interface_obj : mInterfaceObjects)
        {
            if (!interface_obj->NeighborOrApproximationFound())   // check if the interface object already found a neighbor
            {
                rCandidateSendObjects.push_back(interface_obj);
            }
        }
    }

    void PostProcessReceivedResults(const InterfaceObjectConfigure::ContainerType& rCandidateSendObjects,
                                    const std::vector<double>& rDistances,
                                    const std::vector<int>& rPairingIndices) override
    {
        int i = 0;
        for (auto interface_obj : rCandidateSendObjects)
        {
            if (rDistances[i] > -0.5f)   // failed search has value "-1"
            {
                interface_obj->ProcessSearchResult(rDistances[i], rPairingIndices[i], mCommRank);
                mSendObjects[mCommRank].push_back(interface_obj);
            }
            ++i;
        }
    }

    // **********************************************************************
    // Side where we search neighbors aka origin ****************************
    // **********************************************************************
    void StoreSearchResults(const std::vector<double>& rDistances,
                            const std::vector<InterfaceObject::Pointer> TempClosestResults,
                            const std::vector<std::vector<double>> TempShapeFunctionValues) override
    {
        for (std::size_t i = 0; i < rDistances.size(); ++i)
        {
            if (rDistances[i] > -0.5f)   // failed search has value "-1"
            {
                mReceiveObjects[mCommRank].push_back(TempClosestResults[i]);
                mShapeFunctionValues[mCommRank].push_back(TempShapeFunctionValues[i]);
            }
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
        buffer << "InterfaceObjectManagerSerial" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "InterfaceObjectManagerSerial";
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void InitializeSizes()
    {
        int size = mInterfaceObjects.size();
        mSendObjects[mCommRank].reserve(size);
        mReceiveObjects[mCommRank].reserve(size);
        mShapeFunctionValues[mCommRank].reserve(size);
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
    InterfaceObjectManagerSerial& operator=(InterfaceObjectManagerSerial const& rOther);

    //   /// Copy constructor.
    //   InterfaceObjectManagerSerial(InterfaceObjectManagerSerial const& rOther){}


    ///@}

}; // Class InterfaceObjectManagerSerial

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InterfaceObjectManagerSerial& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InterfaceObjectManagerSerial& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_OBJECT_MANAGER_SERIAL_H_INCLUDED  defined
