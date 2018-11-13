//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//	                 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_GENERATE_EMBEDDED_SKIN_UTILITY_H_INCLUDED )
#define  KRATOS_GENERATE_EMBEDDED_SKIN_UTILITY_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/geometry_data.h"
#include "utilities/divide_geometry.h"


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

/// Utility to compute the skin representation from a distance function.
/** Provided either a continuous or discontinuous distance function, this 
 *  utility reconstructs the skin representation coming from such distance
 *  function. This is done by computing the element intersections and saving
 *  them in an empty provided model part. Note that such skin representation
 *  is discontinuous even for a provided continuous distance field.
 */
class KRATOS_API(KRATOS_CORE) GenerateEmbeddedSkinUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GenerateEmbeddedSkinUtility
    KRATOS_CLASS_POINTER_DEFINITION(GenerateEmbeddedSkinUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GenerateEmbeddedSkinUtility(
        ModelPart &rModelPart,
        ModelPart& rSkinModelPart,
        const std::string LevelSetType = "continuous") :
        mrModelPart(rModelPart),
        mrSkinModelPart(rSkinModelPart),
        mLevelSetType(LevelSetType) {};

    /// Destructor.
    virtual ~GenerateEmbeddedSkinUtility() {};

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute();

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
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;
    ModelPart &mrSkinModelPart;
    const std::string mLevelSetType;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void Clear();

    void ComputeElementSkin(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances,
        unsigned int &rTempNodeId,
        unsigned int &rTempCondId,
        Properties::Pointer pCondProp,
        ModelPart::NodesContainerType &rNewNodesVect,
        ModelPart::ConditionsContainerType &rNewCondsVect);

    const bool inline ElementIsSplit(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

    void RenumberAndAddSkinEntities(
        const ModelPart::NodesContainerType &rNewNodesVect,
        const ModelPart::ConditionsContainerType &rNewCondsVect);

    const Vector SetDistancesVector(ModelPart::ElementIterator ItElem);

    /**
     * Sets the the divide geometry utility according to the geometry type.
     * @param pGeometry Pointer to the element geometry
     * @param rNodalDistances Vector containing the distance values
     * @return A pointer to the divide geometry utility
     */
    DivideGeometry::Pointer SetDivideGeometryUtility(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

    /**
     * Sets the new interface condition geometry
     * @param rOriginGeometryType Interface subgeometry type
     * @param rNewNodesArray Nodes that conform the new interface geometry
     * @return A pointer to the new geometry
     */
    Geometry< Node<3> >::Pointer SetNewConditionGeometry(
        const GeometryData::KratosGeometryType &rOriginGeometryType,
        const Condition::NodesArrayType &rNewNodesArray);

    Properties::Pointer SetSkinEntitiesProperties();

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
    GenerateEmbeddedSkinUtility& operator=(GenerateEmbeddedSkinUtility const& rOther) = delete;

    /// Copy constructor.
    GenerateEmbeddedSkinUtility(GenerateEmbeddedSkinUtility const& rOther) = delete;

    ///@}

}; // Class GenerateEmbeddedSkinUtility

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos.

#endif // KRATOS_GENERATE_EMBEDDED_SKIN_UTILITY_H_INCLUDED  defined
