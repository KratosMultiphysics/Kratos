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
#include "voxel_mesher_coloring_factory.h"

#include "modeler/coloring/voxel_mesher_coloring.h"
#include "modeler/coloring/color_cells_with_inside_center.h"
#include "modeler/coloring/color_cells_in_touch.h"
#include "modeler/coloring/color_and_intersect_with_cells_in_touch.h"
#include "modeler/coloring/color_cell_faces.h"
#include "modeler/coloring/color_cell_faces_between_colors.h"
#include "modeler/coloring/color_outer_faces_of_cells_with_colors.h"
#include "modeler/coloring/color_cells_with_center_in_sphere_around_nodes.h"
#include "modeler/coloring/color_connected_cells_in_touch.h"

namespace Kratos {

template class KratosComponents<VoxelMesherColoringFactory::RegisteredType>;

void RegisterVoxelMesherColoring()
{
    VoxelMesherColoringFactory::Register<ColorCellsWithInsideCenter>("cells_with_inside_center");
    VoxelMesherColoringFactory::Register<ColorCellsInTouch>("cells_in_touch");
    VoxelMesherColoringFactory::Register<ColorAndIntersectWithCellsInTouch>("intersect_with_cells_in_touch");
    VoxelMesherColoringFactory::Register<ColorCellFaces>("cells_faces");
    VoxelMesherColoringFactory::Register<ColorCellFacesBetweenColors>("cells_faces_between_colors");
    VoxelMesherColoringFactory::Register<ColorOuterFacesOfCellsWithColors>("outer_faces_of_cells_with_color");
    VoxelMesherColoringFactory::Register<ColorCellsWithCenterInSphereArounNodes>("cells_with_center_in_sphere_arround_nodes");
    VoxelMesherColoringFactory::Register<ColorConnectedCellsInTouch>("connected_cells_in_touch");
}

}