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

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /// Definition of the iterators
    typedef NodesArrayType::iterator NodeItType;
    typedef ElementsArrayType::iterator ElementItType;

    /**
     * @brief This method ensured that the properties of elements and conditions are on the model part (it does recursively in all model parts)
     * @param ModelPart The model part where ensure properties
     * @param RemovePreviousProperties If we clear previous properties and ensure only the properties existing in the elements and conditions (true by default)
     */
    void KRATOS_API(MESHING_APPLICATION) RecursiveEnsureModelPartOwnsProperties(
        ModelPart& rModelPart,
        const bool RemovePreviousProperties = true
        );

    /**
     * @brief This method ensured that the properties of elements and conditions are on the model part
     * @param ModelPart The model part where ensure properties
     * @param RemovePreviousProperties If we clear previous properties and ensure only the properties existing in the elements and conditions (true by default)
     */
    void KRATOS_API(MESHING_APPLICATION) EnsureModelPartOwnsProperties(
        ModelPart& rModelPart,
        const bool RemovePreviousProperties = true
        );

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

}; // namespace MeshingUtilities
}  // namespace Kratos
#endif /* KRATOS_MESHING_UTILITIES_H defined */
