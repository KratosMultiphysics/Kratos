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

#if !defined(KRATOS_INTERFACE_PREPROCESSOR_H)
#define  KRATOS_INTERFACE_PREPROCESSOR_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

#include "mapper_local_system.h"


namespace Kratos
{
///@addtogroup MappingApplication
///@{

///@name Kratos Classes
///@{

/// Preparing the interface by construting the MapperLocalSystems
/** The MapperLocalSystems are used for the mapping-system, this class
 * constructs them based on the available geometry
*/
class InterfacePreprocessor
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of InterfacePreprocessor
    KRATOS_CLASS_POINTER_DEFINITION(InterfacePreprocessor);

    typedef Kratos::unique_ptr<MapperLocalSystem> MapperLocalSystemPointer;
    typedef std::vector<MapperLocalSystemPointer> MapperLocalSystemPointerVector;
    typedef Kratos::shared_ptr<MapperLocalSystemPointerVector> MapperLocalSystemPointerVectorPointer;

    ///@}
    ///@name Life Cycle
    ///@{

    InterfacePreprocessor(ModelPart& rModelPartDestination,
                          MapperLocalSystemPointerVectorPointer pMapperLocalSystems)
        : mrModelPartDestination(rModelPartDestination),
          mpMapperLocalSystems(pMapperLocalSystems)
    {
    }

    /// Destructor.
    virtual ~InterfacePreprocessor() = default;

    ///@}
    ///@name Operations
    ///@{

    void CreateMapperLocalSystems(const MapperLocalSystemPointer& rpLocalSystem);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "InterfacePreprocessor";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const {}

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {}

    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPartDestination;
    MapperLocalSystemPointerVectorPointer mpMapperLocalSystems;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CreateMapperLocalSystemsFromNodes(const MapperLocalSystemPointer& rpLocalSystem);

    void CreateMapperLocalSystemsFromGeometries(const MapperLocalSystemPointer& rpLocalSystem);

    ///@}

}; // Class InterfacePreprocessor

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERFACE_PREPROCESSOR_H  defined
