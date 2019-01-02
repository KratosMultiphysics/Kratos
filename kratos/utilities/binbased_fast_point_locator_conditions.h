//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                       license: license.txt
//
//  License:          BSD License
//  Main authors:     Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_BINBASED_FAST_POINT_LOCATOR_CONDITIONS_INCLUDED )
#define  KRATOS_BINBASED_FAST_POINT_LOCATOR_CONDITIONS_INCLUDED

// System includes

// External includes


// Project includes
#include "utilities/binbased_fast_point_locator.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size definition
    typedef std::size_t SizeType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class BinBasedFastPointLocatorConditions
 * @ingroup KratosCore
 * @brief This class is designed to allow the fast location of MANY points on the top of a 3D mesh. (for conditions)
 * @details The utility relies on the creation of a Bin of objects that allows finding quikly a reduced number of condition candidates for the location of a point.
 * The basic idea is to allow finding the condition in which a given spatial position sits
 * The user should call the function "UpdateSearchDatabase" to mount the bin and subsequently locate the points as needed
 * @author  Vicente Mataix Ferrandiz
 * @note The location function is threadsafe, and can be used in OpenMP loops
 * @tparam TDim If we work in a 2D or 3D space
 */
template< SizeType TDim>
class BinBasedFastPointLocatorConditions
    : public BinBasedFastPointLocator<TDim, SpatialContainersConfigure<TDim, Condition>>
{
public:
    ///@name Type Definitions
    ///@{

    /// The base type definition
    typedef BinBasedFastPointLocator<TDim, SpatialContainersConfigure<TDim, Condition>> BaseType;

    /// Pointer definition of BinBasedFastPointLocatorConditions
    KRATOS_CLASS_POINTER_DEFINITION(BinBasedFastPointLocatorConditions);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief This is the default constructor
     * @param rModelPart The model part of the mesh used in the search
     */
    BinBasedFastPointLocatorConditions(ModelPart& rModelPart)
        : BaseType(rModelPart)
    {
    }

    /// Destructor.
    ~BinBasedFastPointLocatorConditions() override = default;
    
    ///@}
};
    
} // namespace Kratos.

#endif // KRATOS_BINBASED_FAST_POINT_LOCATOR_CONDITIONS_INCLUDED  defined


