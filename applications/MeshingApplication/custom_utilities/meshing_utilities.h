// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_MESHING_UTILITIES_H)
#define KRATOS_MESHING_UTILITIES_H

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

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
/**
 * @namespace MeshingUtilities
 * @ingroup MeshingApplication
 * @brief This namespace includes several utilities necessaries for the computation of the meshing techniques
 * @author Vicente Mataix Ferrandiz
 */
namespace MeshingUtilities
{
    /// The size type definition
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// Definition of the node and geometry
    typedef Node<3> NodeType;
    typedef Geometry<NodeType> GeometryType;

    /// The arrays of elements and nodes
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// Definition of the iterators
    typedef NodesArrayType::iterator NodeItType;
    typedef ElementsArrayType::iterator ElementItType;
    typedef ConditionsArrayType::iterator ConditionItType;

    /**
     * @brief This computes the element size depending of a whole model part and it assigns to the ELEMENT_H variable
     * @param ModelPart The model part where compute the  and block them
     * @param ThisParameters The parameters
     */
    void KRATOS_API(MESHING_APPLICATION) BlockThresholdSizeElements(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief This computes the element size depending of a whole model part and it assigns to the ELEMENT_H variable
     * @param ModelPart The model part where compute the ELEMENT_H
     */
    void KRATOS_API(MESHING_APPLICATION) ComputeElementsSize(ModelPart& rModelPart);

    /**
     * @brief This computes the element size depending of the geometry and it assigns to the ELEMENT_H variable
     * @param itElement The element iterator
     */
    void KRATOS_API(MESHING_APPLICATION) ComputeElementSize(ElementItType itElement);

    /**
     * @brief This process splits a hexahedra mesh into a tetrahedra mesh
     * @details One hexahedra can be splitted into 174 different combinations of tetrahedra (https://arxiv.org/pdf/1801.01288)
     * The node ordering for a hexahedron corresponds with:
     *             v
     *      3----------2
     *      |\     ^   |\
     *      | \    |   | \
     *      |  \   |   |  \
     *      |   7------+---6
     *      |   |  +-- |-- | -> u
     *      0---+---\--1   |
     *       \  |    \  \  |
     *        \ |     \  \ |
     *         \|      w  \|
     *          4----------5
     * The 6 tetrahedra generated on this method corresponds with:
     *      - 0-4-5-7
     *      - 0-5-1-2
     *      - 0-7-5-2
     *      - 0-7-2-3
     *      - 5-7-6-2
     * @param ModelPart The model part to be splitted
     * @param ThisParameters The parameters (additional configurations)
     */
    void KRATOS_API(MESHING_APPLICATION) HexahedraMeshToTetrahedraMesh(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

}; // namespace MeshingUtilities
}  // namespace Kratos
#endif /* KRATOS_MESHING_UTILITIES_H defined */
