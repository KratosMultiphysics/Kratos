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

#if !defined(KRATOS_MAPPER_FACTORY_H_INCLUDED )
#define  KRATOS_MAPPER_FACTORY_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"

#include "custom_mappers/mapper.h"


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

/// Python Interface of the MappingApplication
/** This class constructs the mappers and exposes them to Python
* Some checks are performed to see if the Input (ModelParts and JSON-Parameters) are valid
* For information abt the available echo_levels and the JSON default-parameters
* look into the class description of the MapperCommunicator
*/
class MapperFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MapperFactory
    KRATOS_CLASS_POINTER_DEFINITION(MapperFactory);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~MapperFactory() { }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    

    static Mapper::Pointer CreateMapper(ModelPart& rModelPartOrigin,
                                        ModelPart& rModelPartDestination,
                                        Parameters JsonParameters);

    static void RegisterMapper(const std::string MapperName,
                               Mapper::Pointer pMapperPrototype);


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
        buffer << "MapperFactory" ;
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "MapperFactory";
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

    /// Default constructor.
    MapperFactory() {}

    static ModelPart& ReadInterfaceModelPart(ModelPart& rModelPart,
                                             Parameters InterfaceParameters,
                                             const std::string& InterfaceSide);

    static std::unordered_map<std::string, Mapper::Pointer>& GetRegisteredMappersList();

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{s


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    MapperFactory& operator=(MapperFactory const& rOther);

    //   /// Copy constructor.
    //   MapperFactory(MapperFactory const& rOther){}


    ///@}

}; // Class MapperFactory

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  MapperFactory& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MapperFactory& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MAPPER_FACTORY_H_INCLUDED  defined