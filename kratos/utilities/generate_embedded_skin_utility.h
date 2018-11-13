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
        ModelPart &rSkinModelPart,
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

    /**
     * @brief Call to generate the embedded skin model part
     * This method collects all the operations required to generate
     * the embedded skin model part. The new geometries will be stored 
     * inside the skin model part provided in the constructor. 
     */
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

    /**
     * @brief Geometry clear operation
     * This method is called before any construction of the skin.
     * It removes the existent nodes, elements and conditions in
     * the provided skin model part.
     */
    void Clear();

    /**
     * @brief Computes the skin entities of one element
     * For an intersected element, this method computes the skin
     * intersection geometries.
     * @param rGeometry geometry of the element of interest
     * @param rNodalDistances vector containing the element node distances
     * @param rTempNodeId temporal id for the new nodes
     * @param rTempCondId temporal id for the new conditions
     * @param pCondProp pointer to the properties for the new skin conditions
     * @param rNewNodesVect vector to temporary store the new skin nodes
     * @param rNewCondsVect vector to temporary store the new skin conditions
     */
    void ComputeElementSkin(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances,
        unsigned int &rTempNodeId,
        unsigned int &rTempCondId,
        Properties::Pointer pCondProp,
        ModelPart::NodesContainerType &rNewNodesVect,
        ModelPart::ConditionsContainerType &rNewCondsVect);

    /**
     * @brief Checks if an element is split
     * This method checks if an element geometry is split
     * @param rGeometry geometry of the element of interest
     * @param rNodalDistances vector containing the element node distances
     * @return true if the element is split
     * @return false if the element is not split
     */
    const bool inline ElementIsSplit(
        const Geometry<Node<3>> &rGeometry,
        const Vector &rNodalDistances);

    /**
     * @brief Renumber and saves the new skin entities
     * This method renumbers the new skin geometrical entities (MPI compatible)
     * and add them to the provided skin model part.
     * @param rNewNodesVect vector that stores the new skin nodes
     * @param rNewCondsVect vector that stores the new skin conditions
     */
    void RenumberAndAddSkinEntities(
        const ModelPart::NodesContainerType &rNewNodesVect,
        const ModelPart::ConditionsContainerType &rNewCondsVect);

    /**
     * @brief Set the Distances Vector object
     * For a given element, this method sets the vector containing the
     * element node distances.
     * @param ItElem iterator to the element of interest
     * @return const Vector vector containing the element node distances
     */
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

    /**
     * @brief Set the Skin Entities Properties 
     * This method checks which is the last properties id. and 
     * sets a new one accordingly to be used as skin conditions property
     * @return Properties::Pointer pointer to the new skin entities property
     */
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
