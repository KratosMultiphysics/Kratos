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
class MapperInterfaceInfo
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperInterfaceInfo
    KRATOS_CLASS_POINTER_DEFINITION(MapperInterfaceInfo);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MapperInterfaceInfo();

    /// Destructor.
    virtual ~MapperInterfaceInfo() {
        std::cout << "Destructor of MapperInterfaceInfo called" << std::endl;
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // void ProcessSearchResult(InterfaceObject::Pointer pInterfaceObject)
    // {

    // }


    std::vector<int> GetNeighborIds() const
    {
        std::vector<int> neighbor_ids(1);
        // neighbor_ids[0] = mNeighborId;
        return neighbor_ids;
    }

    std::vector<double> GetNeighborDistances() const
    {
        std::vector<double> neighbor_distances(1);
        // neighbor_distances[0] = mNeighborDistance;
        return neighbor_distances;
    }

    MapperInterfaceInfo::Pointer Create() const
    {
        return Kratos::make_shared<MapperInterfaceInfo>();
    }

    void Clear() {}


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
        return "MapperInterfaceInfo";
    }

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

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const {
        KRATOS_ERROR << "This Object cannot be serialized!" << std::endl;
    }

    virtual void load(Serializer& rSerializer) {
        KRATOS_ERROR << "This Object cannot be serialized!" << std::endl;
    }


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    // MapperInterfaceInfo& operator=(MapperInterfaceInfo const& rOther) {}

    /// Copy constructor.
    // MapperInterfaceInfo(MapperInterfaceInfo const& rOther) {}

    ///@}

}; // Class MapperInterfaceInfo

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

inline std::istream & operator >> (std::istream& rIStream, MapperInterfaceInfo& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const MapperInterfaceInfo& rThis) {
//   rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
//   rThis.PrintData(rOStream);
  return rOStream;
}

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const std::vector<MapperInterfaceInfo::Pointer>& rThis) {
//   rThis.PrintInfo(rOStream);
  rOStream << " : " << std::endl;
//   rThis.PrintData(rOStream);
  return rOStream;
}


///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_INTERFACE_INFO_H_INCLUDED  defined
