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

// System includes

// External includes

// Project includes
#include "utilities/auxiliar_model_part_utilities.h"
#include "custom_utilities/meshing_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
namespace MeshingUtilities
{

void BlockThresholdSizeElements(
    ModelPart& rModelPart,
    Parameters ThisParameters
    )
{
    Parameters default_parameters = Parameters(R"(
    {
        "minimal_size" : 0.1,
        "maximal_size" : 10.0
    })" );

    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Getting values
    const double minimal_size = ThisParameters["minimal_size"].GetDouble();
    const double maximal_size = ThisParameters["maximal_size"].GetDouble();

    // Compute elements size
    ComputeElementsSize(rModelPart);

    // Iterating over the elements
    block_for_each(rModelPart.Elements(),
        [&minimal_size,&maximal_size](Element& rElement) {

        if (rElement.IsNot(BLOCKED)) {
            // Getting ELEMENT_H
            const double element_h = rElement.GetValue(ELEMENT_H);

            // Blocking elements in the threshold sizes
            if (element_h <= minimal_size || element_h >= maximal_size) {
                rElement.Set(BLOCKED, true);
            }
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeElementsSize(ModelPart& rModelPart)
{
    // Iterating over the elements
    block_for_each(rModelPart.Elements(),
        [&](Element& rElement) {
        ComputeElementSize(rElement);
    });
}

/***********************************************************************************/
/***********************************************************************************/

void ComputeElementSize(Element& rElement)
{
    auto& this_geometry = rElement.GetGeometry();

    // Here we compute the element size. This process is designed for triangles and tetrahedra, so we only specify for this geometries. Otherwise we take the length (and we throw a warning)
    if (this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Triangle2D3){ // Triangular elements
        rElement.SetValue(ELEMENT_H, 2.0 * this_geometry.Circumradius());
    } else if(this_geometry.GetGeometryType() == GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4){ // Tetrahedral elements
        rElement.SetValue(ELEMENT_H,std::pow(12.0 * this_geometry.Volume()/std::sqrt(2.0), 1.0/3.0));
    } else { // In any othe case just considers the length of the element
        KRATOS_WARNING("MetricErrorProcess") << "This process is designed for tetrahedra (3D) and triangles (2D). Error expected" << std::endl;
        rElement.SetValue(ELEMENT_H, this_geometry.Length());
    }
}

} // namespace MeshingUtilities
} // namespace Kratos
