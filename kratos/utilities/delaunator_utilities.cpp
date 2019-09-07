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

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

#ifdef INCLUDE_DELAUNATORCPP
#include "delaunator.hpp"
#endif

// Project includes
#include "utilities/delaunator_utilities.h"

namespace Kratos
{
namespace DelaunatorUtilities
{
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

/***********************************************************************************/
/***********************************************************************************/

void CreateTriangleMeshFromNodes(ModelPart& rModelPart)
{
    KRATOS_TRY

    // Ensure node order
    auto& r_nodes_root_array = rModelPart.GetRootModelPart().Nodes();
    const auto it_node_root_begin = r_nodes_root_array.begin();
    #pragma omp parallel for
    for(int i=0; i<static_cast<int>(r_nodes_root_array.size()); ++i) {
        auto it_node = it_node_root_begin + i;
        it_node->SetId(i + 1);
    }

    // Getting nodes array
    const auto& r_nodes_array = rModelPart.Nodes();
    const auto it_node_begin = r_nodes_array.begin();

    // The vector to fill
    std::vector<double> coordinates;

    // NOTE: 2D asssumed

    #pragma omp parallel
    {
        std::vector<double> coordinates_buffer;

        #pragma omp for
        for(int i=0; i<static_cast<int>(r_nodes_array.size()); ++i) {
            auto it_node = it_node_begin + i;

            // Filling coordinates buffer
            coordinates_buffer.push_back(it_node->X());
            coordinates_buffer.push_back(it_node->Y());
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(coordinates_buffer.begin(),coordinates_buffer.end(),back_inserter(coordinates));
        }
    }

    // Creating the triangles
    Properties::Pointer p_elem_prop = rModelPart.CreateNewProperties(0);
    const auto& r_triangles = ComputeTrianglesConnectivity(coordinates);
    std::size_t counter = rModelPart.GetRootModelPart().Elements().size() + 1;
    for (std::size_t i = 0; i < r_triangles.size(); i += 3) {
        rModelPart.CreateNewElement("Element2D3N", counter, {{r_triangles[i] + 1,r_triangles[i + 1] + 1, r_triangles[i + 2] + 1}}, p_elem_prop);
        ++counter;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::size_t> ComputeTrianglesConnectivity(const std::vector<double>& rCoordinates)
{
#ifdef INCLUDE_DELAUNATORCPP
    // Calling the delaunator library
    delaunator::Delaunator delaunator(rCoordinates);

    // Creating the triangles
    return delaunator.triangles;
#else
    // Use triangle library
    return ComputeTrianglesConnectivityWithTriangle(rCoordinates);
#endif
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<std::size_t> ComputeTrianglesConnectivityWithTriangle(const std::vector<double>& rCoordinates)
{
    // Creating the containers for the input and output
    struct triangulateio in_mid;
    struct triangulateio out_mid;
    struct triangulateio vorout_mid;

    InitializeTriangulateIO(in_mid);
    InitializeTriangulateIO(out_mid);
    InitializeTriangulateIO(vorout_mid);

    in_mid.numberofpoints = rCoordinates.size()/2;
    in_mid.pointlist = (REAL *) malloc(in_mid.numberofpoints * 2 * sizeof(REAL));

    for(std::size_t i = 0; i < rCoordinates.size(); i++) {
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
}

} // namespace DelaunatorUtilities
} // namespace Kratos
