//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "voxel_mesher_entity_generation_factory.h"
#include "generate_elements_with_cell_color.h"
#include "generate_tetrahedral_elements_with_cell_color.h"
#include "generate_quadratic_elements_with_cell_color.h"
#include "generate_conditions_with_face_color.h"
#include "generate_triangular_conditions_with_face_color.h"
#include "generate_stl_intersection_with_cells.h"

namespace Kratos {

template class KratosComponents<VoxelMesherEntityGenerationFactory::RegisteredType>;

void RegisterVoxelMesherEntityGeneration()
{
    VoxelMesherEntityGenerationFactory::Register<GenerateElementsWithCellColor>("elements_with_cell_color");
    VoxelMesherEntityGenerationFactory::Register<GenerateTetrahedralElementsWithCellColor>("tetrahedral_elements_with_cell_color");
    VoxelMesherEntityGenerationFactory::Register<GenerateQuadraticElementsWithCellColor>("quadratic_elements_with_cell_color");
    VoxelMesherEntityGenerationFactory::Register<GenerateConditionsWithFaceColor>("conditions_with_face_color");
    VoxelMesherEntityGenerationFactory::Register<GenerateTriangularConditionsWithFaceColor>("triangular_conditions_with_face_color");
    VoxelMesherEntityGenerationFactory::Register<GenerateStlIntersectionWithCells>("stl_intersection_with_cells");
}

}