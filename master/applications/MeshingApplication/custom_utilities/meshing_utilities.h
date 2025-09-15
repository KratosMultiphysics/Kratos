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
     * @param rElement The element reference
     */
    void KRATOS_API(MESHING_APPLICATION) ComputeElementSize(Element& rElement);

}; // namespace MeshingUtilities
}  // namespace Kratos
#endif /* KRATOS_MESHING_UTILITIES_H defined */
