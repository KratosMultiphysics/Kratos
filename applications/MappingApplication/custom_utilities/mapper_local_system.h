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

#if !defined(KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED )
#define  KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "mapper_interface_info.h"


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

/// This is the "Condition" of the Mappers.
/** Detail class definition.
*/
class MapperLocalSystem
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperLocalSystem
    KRATOS_CLASS_POINTER_DEFINITION(MapperLocalSystem);

    using MapperLocalSystemUniquePointer = std::unique_ptr<MapperLocalSystem>;
    using MapperInterfaceInfoPointer = MapperInterfaceInfo::Pointer;

    using MappingWeightsVector = std::vector<double>;
    using OriginIdVector       = std::vector<int>;
    using DestinationIdVector  = std::vector<int>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperLocalSystem();

    /// Destructor.
    virtual ~MapperLocalSystem() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual MapperLocalSystemUniquePointer Create() = 0;

    virtual void CalculateAll(MappingWeightsVector& rMappingWeights,
                              OriginIdVector&       rOriginIds,
                              DestinationIdVector&  rDestinationIds) = 0;

    void AddInterfaceInfo(MapperInterfaceInfoPointer pInterfaceInfo)
    {
        mInterfaceInfos.push_back(pInterfaceInfo);
    }

    virtual void Clear()
    {
        mInterfaceInfos.clear();
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
    virtual std::string Info() const {return "MapperLocalSystem";}

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

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

    std::vector<MapperInterfaceInfoPointer> mInterfaceInfos;


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


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    // /// Assignment operator.
    // MapperLocalSystem& operator=(MapperLocalSystem const& rOther);

    // /// Copy constructor.
    // MapperLocalSystem(MapperLocalSystem const& rOther);


    ///@}

}; // Class MapperLocalSystem

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                MapperLocalSystem& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const MapperLocalSystem& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_LOCAL_SYSTEM_H_INCLUDED  defined
