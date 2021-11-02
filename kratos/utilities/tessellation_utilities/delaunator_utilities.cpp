//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#define REAL double
#ifdef USE_TRIANGLE_NONFREE_TPL
#include "triangle.h"
#else 
#include <delaunator.hpp>
#endif

// Project includes
#include "includes/model_part.h"
#include "utilities/tessellation_utilities/delaunator_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
namespace DelaunatorUtilities
{
#ifdef USE_TRIANGLE_NONFREE_TPL
extern "C" {
    void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
}

/**
* @brief This method initializes the triangulate IO
* @param rTriangles The triangles to be initialized
*/
void InitializeTriangulateIO( triangulateio& rTriangles )
{
    rTriangles.pointlist                  = (REAL*) nullptr;
    rTriangles.pointattributelist         = (REAL*) nullptr;
    rTriangles.pointmarkerlist            = (int*) nullptr;
    rTriangles.numberofpoints             = 0;
    rTriangles.numberofpointattributes    = 0;
    rTriangles.trianglelist               = (int*) nullptr;
    rTriangles.triangleattributelist      = (REAL*) nullptr;
    rTriangles.trianglearealist           = (REAL*) nullptr;
    rTriangles.neighborlist               = (int*) nullptr;
    rTriangles.numberoftriangles          = 0;
    rTriangles.numberofcorners            = 3;
    rTriangles.numberoftriangleattributes = 0;
    rTriangles.segmentlist                = (int*) nullptr;
    rTriangles.segmentmarkerlist          = (int*) nullptr;
    rTriangles.numberofsegments           = 0;
    rTriangles.holelist                   = (REAL*) nullptr;
    rTriangles.numberofholes              = 0;
    rTriangles.regionlist                 = (REAL*) nullptr;
    rTriangles.numberofregions            = 0;
    rTriangles.edgelist                   = (int*) nullptr;
    rTriangles.edgemarkerlist             = (int*) nullptr;
    rTriangles.normlist                   = (REAL*) nullptr;
    rTriangles.numberofedges              = 0;
};

/**
* @brief This method cleans the triangulate IO
* @param rTriangles The triangles to be cleaned
*/
void CleanTriangulateIO( triangulateio& rTriangles )
{
    if(rTriangles.pointlist != nullptr) free(rTriangles.pointlist );
    if(rTriangles.pointattributelist != nullptr) free(rTriangles.pointattributelist );
    if(rTriangles.pointmarkerlist != nullptr) free(rTriangles.pointmarkerlist   );
    if(rTriangles.trianglelist != nullptr) free(rTriangles.trianglelist  );
    if(rTriangles.triangleattributelist != nullptr) free(rTriangles.triangleattributelist );
    if(rTriangles.trianglearealist != nullptr) free(rTriangles.trianglearealist );
    if(rTriangles.neighborlist != nullptr) free(rTriangles.neighborlist   );
    if(rTriangles.segmentlist != nullptr) free(rTriangles.segmentlist    );
    if(rTriangles.segmentmarkerlist != nullptr) free(rTriangles.segmentmarkerlist   );
    if(rTriangles.holelist != nullptr) free(rTriangles.holelist      );
    if(rTriangles.regionlist != nullptr) free(rTriangles.regionlist  );
    if(rTriangles.edgelist != nullptr) free(rTriangles.edgelist   );
    if(rTriangles.edgemarkerlist != nullptr) free(rTriangles.edgemarkerlist   );
    if(rTriangles.normlist != nullptr) free(rTriangles.normlist  );
};
#endif
/***********************************************************************************/
/***********************************************************************************/

void CreateTriangleMeshFromNodes(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Ensure node order
    auto& r_nodes_root_array = rModelPart.GetRootModelPart().Nodes();
    const auto it_node_root_begin = r_nodes_root_array.begin();
    IndexPartition<std::size_t>(r_nodes_root_array.size()).for_each(
        [&it_node_root_begin](std::size_t i_node)
        { (it_node_root_begin + i_node)->SetId(i_node + 1); }
    );

    // Getting nodes array
    const auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // The vector to fill
    const std::size_t number_of_nodes = r_nodes_array.size();
    std::vector<double> coordinates(2*number_of_nodes, 0.0);

    // NOTE: 2D asssumed
    for(std::size_t i=0; i<number_of_nodes; ++i) {
        auto it_node = it_node_begin + i;

        // Filling coordinates buffer
        coordinates[2*i]   = it_node->X();
        coordinates[2*i+1] = it_node->Y();
    }

    // Creating the triangles
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    const auto& r_triangles = ComputeTrianglesConnectivity(coordinates);
    std::size_t counter = rModelPart.GetRootModelPart().Elements().size() + 1;
    for (std::size_t i = 0; i < r_triangles.size(); i += 3) {
        rModelPart.CreateNewElement("Element2D3N", counter, {{r_triangles[i] + 1,r_triangles[i + 1] + 1, r_triangles[i + 2] + 1}}, p_elem_prop);
        ++counter;
    }

    // Check orientation
    for (auto& r_elem : rModelPart.Elements()) {
        if(r_elem.GetGeometry().Area() < 0.0) {
            r_elem.GetGeometry()(0).swap(r_elem.GetGeometry()(1));
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::size_t> ComputeTrianglesConnectivity(const std::vector<double>& rCoordinates)
{
#ifdef USE_TRIANGLE_NONFREE_TPL
    // Creating the containers for the input and output
    struct triangulateio in_mid, out_mid, vorout_mid;

    InitializeTriangulateIO(in_mid);
    InitializeTriangulateIO(out_mid);
    InitializeTriangulateIO(vorout_mid);

    in_mid.numberofpoints = rCoordinates.size()/2;
    in_mid.pointlist = (REAL *) malloc(in_mid.numberofpoints * 2 * sizeof(REAL));

    for(std::size_t i = 0; i < rCoordinates.size(); ++i) {
        in_mid.pointlist[i] = rCoordinates[i];
    }

    // "P" suppresses the output .poly file. Saves disk space, but you
    // lose the ability to maintain constraining segments  on later refinements of the mesh.
    // "B" Suppresses boundary markers in the output .node, .poly, and .edge output files
    // "e" outputs edge list (i.e. all the "connectivities")
    // "Q" Quiet:  No terminal output except errors.
    // "z" Numbers all items starting from zero (rather than one)
    // "c" Encloses the convex hull with segments
    // "D" Conforming Delaunay:  all triangles are truly Delaunay
    char options1[] = "QPez";
    triangulate(options1, &in_mid, &out_mid, &vorout_mid);

    const std::size_t number_of_triangles = out_mid.numberoftriangles;

    std::vector<std::size_t> connectivities(3 * number_of_triangles);

    const auto& r_triangles_list = out_mid.trianglelist;

    // Must be copied into a std::vector
    for (std::size_t i = 0; i < number_of_triangles * 3; i += 3) {
        connectivities[i] = r_triangles_list[i];
        connectivities[i + 1] = r_triangles_list[i + 1];
        connectivities[i + 2] = r_triangles_list[i + 2];
    }

    CleanTriangulateIO(in_mid);
    CleanTriangulateIO(out_mid);
    CleanTriangulateIO(vorout_mid);

    return connectivities;
#else
    // Calling the library Delaunator
    delaunator::Delaunator delaunator(rCoordinates);
    const auto& r_triangles = delaunator.triangles;
    return r_triangles;
#endif
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::size_t> ComputeTrianglesConnectivity(const std::vector<Point>& rPoints)
{
    const std::size_t number_of_nodes = rPoints.size();
    std::vector<double> coordinates(2*number_of_nodes, 0.0);

    // NOTE: 2D asssumed
    // Filling coordinates buffer
    for(std::size_t i=0; i<number_of_nodes; ++i) {
        const auto& r_point = rPoints[i];
        coordinates[2*i]   = r_point.X();
        coordinates[2*i+1] = r_point.Y();
    }
    return ComputeTrianglesConnectivity(coordinates);
}

/***********************************************************************************/
/***********************************************************************************/

std::pair<std::vector<std::size_t>, std::vector<double>> ComputeTrianglesConnectivity(
    const std::vector<double>& rCoordinates,
    const std::vector<std::array<double,2>>& rSegments,
    const double AreaConstraint
    )
{
#ifdef USE_TRIANGLE_NONFREE_TPL
    // Creating the containers for the input and output
    struct triangulateio in_mid, out_mid, vorout_mid;

    InitializeTriangulateIO(in_mid);
    InitializeTriangulateIO(out_mid);
    InitializeTriangulateIO(vorout_mid);

    // Initialize the boundary points coordinates list
    // Note 1: that InitializeTriangulateIO allocates nothing for this by default
    // Note 2: this will be deallocated within the ClearTriangulateIO call below
    in_mid.numberofpoints = rCoordinates.size()/2;
    in_mid.pointlist = (REAL *) malloc(in_mid.numberofpoints * 2 * sizeof(REAL));

    for(std::size_t i = 0; i < rCoordinates.size(); ++i) {
        in_mid.pointlist[i] = rCoordinates[i];
    }

    // Initilize the segment list (note that default zero markers are assumed)
    in_mid.numberofsegments = rSegments.size();
    in_mid.segmentlist = (int*) malloc(in_mid.numberofsegments * 2 * sizeof(int));
    in_mid.segmentmarkerlist = (int*) malloc(in_mid.numberofsegments * sizeof(int));
    for (std::size_t i = 0; i < rSegments.size(); ++i) {
        const auto& r_segment = rSegments[i];
        in_mid.segmentlist[2*i] = r_segment[0];
        in_mid.segmentlist[2*i + 1] = r_segment[1];
    }

    // Check https://www.cs.cmu.edu/~quake/triangle.switch.html for a detailed Triangle switches description
    // "Q" quiet (no terminal output except errors)
    // "q" quality mesh generation with no angles smaller than 20 degrees
    // "p" triangulates a Planar Straight Line Graph
    // "z" numbers all items starting from zero (rather than one)
    // "a" imposes a maximum triangle area constrain
    std::string meshing_options = AreaConstraint > 0.0 ? "Qqpza" + std::to_string(AreaConstraint) : "Qqpz";
    triangulate(&meshing_options[0], &in_mid, &out_mid, &vorout_mid);

    // Save the obtained connectivities in an output std::vector
    const auto& r_triangles_list = out_mid.trianglelist;
    const std::size_t number_of_triangles = out_mid.numberoftriangles;
    std::vector<std::size_t> connectivities(3 * number_of_triangles);
    IndexPartition<std::size_t>(number_of_triangles).for_each([&](std::size_t iTriangle){
        connectivities[3 * iTriangle] = r_triangles_list[3 * iTriangle];
        connectivities[3 * iTriangle + 1] = r_triangles_list[3 * iTriangle + 1];
        connectivities[3 * iTriangle + 2] = r_triangles_list[3 * iTriangle + 2];
    });

    // Save the obtained coordinates in an output std::vector
    const auto& r_out_points_list = out_mid.pointlist;
    const std::size_t number_of_output_points = out_mid.numberofpoints;
    std::vector<double> output_coordinates(2 * number_of_output_points);
    IndexPartition<std::size_t>(number_of_output_points).for_each([&](std::size_t iPoint){
        output_coordinates[2 * iPoint] = r_out_points_list[2 * iPoint];
        output_coordinates[2 * iPoint + 1] = r_out_points_list[2 * iPoint + 1];
    });

    // Clean the triangle database
    CleanTriangulateIO(in_mid);
    CleanTriangulateIO(out_mid);
    CleanTriangulateIO(vorout_mid);

    return std::make_pair(connectivities, output_coordinates);
#else
    KRATOS_ERROR << "The current implementation requires Triangle. Please avoid this utility" << std::endl;
    std::vector<std::size_t> connectivities;
    std::vector<double> output_coordinates;
    return std::make_pair(connectivities, output_coordinates);
#endif
}

} // namespace DelaunatorUtilities
} // namespace Kratos

#undef REAL
