//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//

#if !defined(CHIMERA_HOLE_CUTTING_UTILITY_H_INCLUDED)
#define CHIMERA_HOLE_CUTTING_UTILITY_H_INCLUDED

// System includes
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"                     // Skin face geometry template
#include "geometries/line_2d_2.h"
// Application includes
#include "chimera_application_variables.h"

namespace Kratos
{

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
class KRATOS_API(CHIMERA_APPLICATION) ChimeraHoleCuttingUtility
{
public:

    typedef std::size_t IndexType;
    ///@name Type Definitions
    ///@{

    enum SideToExtract
    {
        INSIDE=0,
        OUTSIDE=1
    };

    enum Domain
    {
        MAIN_BACKGROUND=1,
        OTHER=-1
    };

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ChimeraHoleCuttingUtility
    KRATOS_CLASS_POINTER_DEFINITION(ChimeraHoleCuttingUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    ChimeraHoleCuttingUtility() = default;

    /// Destructor.
    ~ChimeraHoleCuttingUtility() = default;


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a hole on the given rModelPart at a Distance (-ve) from the zero distance layer.
     * @param rModelPart The modelpart where the hole is to be cut.
     * @param rHoleModelPart The modelpart containing the nodes and corresponding elements from the cut hole. (deactivated elements)
     * @param rHoleBoundaryModelPart Boundary modelpart of the rHoleModelPart
     * @param Distance is the the distance (magnitude) at which hole is to be cut from the zero distance layer.
     */
    template<int TDim>
    void CreateHoleAfterDistance(ModelPart &rModelPart,
                                 ModelPart &rHoleModelPart,
                                 ModelPart &rHoleBoundaryModelPart,
                                 const double Distance);


    /**
     * @brief Removes the elements which are out of the domain.
     *          An element is removed even if one of its nodes is out of the domain (-ve or +ve) as indicated by GetInside and MainDomainOrNot
     * @param rModelPart The modelpart From where the elements are to be removed.
     * @param rModifiedModelPart The modified modelpart without the the elements which are out side.
     * @param DomainType says which sign (-ve or +ve) is inside
     * @param OverLapDistance is the the distance (magnitude) at which hole is to be cut from the zero distance layer.
     * @param GetInside works in combination with MainDomainOrNot to get the feeling of what is inside or what is outside.
     */
    template<int TDim>
    void RemoveOutOfDomainElements(ModelPart &rModelPart,
                                   ModelPart &rModifiedModelPart,
                                   const ChimeraHoleCuttingUtility::Domain DomainType,
                                   const double OverLapDistance=0.0,
                                   const ChimeraHoleCuttingUtility::SideToExtract Side = ChimeraHoleCuttingUtility::SideToExtract::OUTSIDE);

    /**
     * @brief Extracts the outside surface/edges of a modelpart.This uses the flag CHIMERA_INTERNAL_BOUNDARY
     *                  to check if there is an internal boundary in the given ModelPart. The flag GetInternal
     *                  specifies weather to get the internal boundary marked by CHIMERA_INTERNAL_BOUNDARY or the outside one.
     * @param rVolumeModelPart The modelpart on which the boundary is to be found.
     * @param rExtractedBoundaryModelPart The extracted surface/edge modelpart.
     * @param GetInternal A bool specifying which surface/edge extracted. The one marked by CHIMERA_INTERNAL_BOUNDARY or the outside one.
     */
    template<int TDim>
    void ExtractBoundaryMesh( ModelPart &rVolumeModelPart,
                              ModelPart &rExtractedBoundaryModelPart,
                              const ChimeraHoleCuttingUtility::SideToExtract GetInternal = ChimeraHoleCuttingUtility::SideToExtract::OUTSIDE);

    /// Assignment operator.
    ChimeraHoleCuttingUtility &operator=(ChimeraHoleCuttingUtility const &rOther) = delete;

    /// Copy constructor.
    ChimeraHoleCuttingUtility(ChimeraHoleCuttingUtility const& rOther) = delete;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

}; // Class ChimeraHoleCuttingUtility

} // namespace Kratos.

#endif // CHIMERA_HOLE_CUTTING_UTILITY_H_INCLUDED  defined
