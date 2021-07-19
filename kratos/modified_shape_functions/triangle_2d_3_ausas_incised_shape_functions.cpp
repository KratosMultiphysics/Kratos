//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes

// External includes

// Project includes
#include "modified_shape_functions/triangle_2d_3_ausas_incised_shape_functions.h"

namespace Kratos
{

/// Triangle2D3AusasIncisedShapeFunctions implementation
/// Default constructor
Triangle2D3AusasIncisedShapeFunctions::Triangle2D3AusasIncisedShapeFunctions(
    const GeometryPointerType pInputGeometry,
    const Vector& rNodalDistancesWithExtrapolated,
    const Vector& rExtrapolatedEdgeRatios)
    : Triangle2D3AusasModifiedShapeFunctions(pInputGeometry, rNodalDistancesWithExtrapolated)
    , mExtraEdgeRatios(rExtrapolatedEdgeRatios)
{};

/// Destructor
Triangle2D3AusasIncisedShapeFunctions::~Triangle2D3AusasIncisedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle2D3AusasIncisedShapeFunctions::Info() const
{
    return "Triangle2D3N Ausas incised shape functions computation class.";
}

/// Print information about this object.
void Triangle2D3AusasIncisedShapeFunctions::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "Triangle2D3N Ausas incised shape functions computation class.";
}

/// Print object's data.
void Triangle2D3AusasIncisedShapeFunctions::PrintData(std::ostream& rOStream) const
{
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    const Vector extra_edge_ratios = this->GetExtrapolatedEdgeRatios();

    rOStream << "Triangle2D3N Ausas incised shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::ostringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tNodal distance values including extrapolated intersections: " << distances_buffer.str() << "\n";
    std::stringstream ratios_buffer;
    std::ostringstream stm2;
    for (unsigned int i = 0; i < extra_edge_ratios.size(); ++i) {
        stm2 << extra_edge_ratios(i);
        ratios_buffer << stm2.str() << " ";
    }
    rOStream << "\tEdge ratios of extrapolated intersections: " << ratios_buffer.str();
}

// Returns the nodal distances vector.
const Vector& Triangle2D3AusasIncisedShapeFunctions::GetExtrapolatedEdgeRatios() const
{
    return mExtraEdgeRatios;
}

void Triangle2D3AusasIncisedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    Triangle2D3AusasIncisedShapeFunctions::SetPositiveSideCondensationMatrix(
        rPosSideCondMatrix,
        Triangle2D3AusasModifiedShapeFunctions::mpTriangleSplitter->mEdgeNodeI,
        Triangle2D3AusasModifiedShapeFunctions::mpTriangleSplitter->mEdgeNodeJ,
        Triangle2D3AusasModifiedShapeFunctions::mpTriangleSplitter->mSplitEdges);
}

void Triangle2D3AusasIncisedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    Triangle2D3AusasIncisedShapeFunctions::SetNegativeSideCondensationMatrix(
        rNegSideCondMatrix,
        Triangle2D3AusasModifiedShapeFunctions::mpTriangleSplitter->mEdgeNodeI,
        Triangle2D3AusasModifiedShapeFunctions::mpTriangleSplitter->mEdgeNodeJ,
        Triangle2D3AusasModifiedShapeFunctions::mpTriangleSplitter->mSplitEdges);
}

// Sets the condensation matrix to transform the subdivsion positive side values to entire element ones.
void Triangle2D3AusasIncisedShapeFunctions::SetPositiveSideCondensationMatrix(
    Matrix& rPosSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges)
{
    const std::size_t n_nodes = 3;
    const std::size_t n_edges = 3;

    // Initialize intersection points condensation matrix
    rPosSideCondMatrix = ZeroMatrix(n_nodes + n_edges, n_nodes);

    // Get the nodal distances vector
    const Vector& nodal_distances = this->GetNodalDistances();

    // Fill the original geometry points main diagonal
    for (std::size_t i = 0; i < n_nodes; ++i) {
        rPosSideCondMatrix(i,i) = (nodal_distances(i) > 0.0) ? 1.0 : 0.0;
    }

    // Compute the intersection points contributions
    for (std::size_t id_edge = 0; id_edge < n_edges; ++id_edge) {
        // Check if the edge has an intersection point
        if (rSplitEdges[n_nodes+id_edge] != -1) {
            // Get the nodes that compose the edge
            const std::size_t edge_node_i = rEdgeNodeI[id_edge];
            const std::size_t edge_node_j = rEdgeNodeJ[id_edge];

            // Transform definition of edge ID from divide_triangle_2d_3.h to definition of triangle_2d_3.h (geometry)
            unsigned int ratio_edge_id = edge_id_for_geometry[id_edge];
            // Check if edge is intersected by extrapolated skin geometry (or actual skin geometry)
            if (mExtraEdgeRatios[ratio_edge_id] > 0.0) {
                // Set shape function value according to the edge ratio of the extrapolated intersection.
                rPosSideCondMatrix(n_nodes+id_edge, node_ids_for_geometry[id_edge][0]) = 1.0 - mExtraEdgeRatios[ratio_edge_id];
                rPosSideCondMatrix(n_nodes+id_edge, node_ids_for_geometry[id_edge][1]) = mExtraEdgeRatios[ratio_edge_id];
            } else {
                // Set to one the shape function value along the positive side of the edge.
                rPosSideCondMatrix(n_nodes+id_edge, edge_node_i) = (nodal_distances(edge_node_i) > 0.0) ? 1.0 : 0.0;
                rPosSideCondMatrix(n_nodes+id_edge, edge_node_j) = (nodal_distances(edge_node_j) > 0.0) ? 1.0 : 0.0;
            }
        }
    }
}

// Sets the condensation matrix to transform the subdivsion negative side values to entire element ones.
void Triangle2D3AusasIncisedShapeFunctions::SetNegativeSideCondensationMatrix(
    Matrix& rNegSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges)
{
    const std::size_t n_nodes = 3;
    const std::size_t n_edges = 3;

    // Initialize intersection points condensation matrix
    rNegSideCondMatrix = ZeroMatrix(n_nodes + n_edges, n_nodes);

    // Get the nodal distances vector
    const Vector& nodal_distances = this->GetNodalDistances();

    // Fill the original geometry points main diagonal
    for (std::size_t i = 0; i < n_nodes; ++i) {
        rNegSideCondMatrix(i,i) = (nodal_distances(i) < 0.0) ? 1.0 : 0.0;
    }

    // Compute the intersection points contributions
    for (std::size_t id_edge = 0; id_edge < n_edges; ++id_edge) {
        // Check if the edge has an intersection point
        if (rSplitEdges[n_nodes+id_edge] != -1) {
            // Get the nodes that compose the edge
            const std::size_t edge_node_i = rEdgeNodeI[id_edge];
            const std::size_t edge_node_j = rEdgeNodeJ[id_edge];

            // Transform definition of edge ID from divide_triangle_2d_3.h to definition of triangle_2d_3.h (geometry)
            unsigned int ratio_edge_id = edge_id_for_geometry[id_edge];
            // Check if edge is intersected by extrapolated skin geometry (or actual skin geometry)
            if (mExtraEdgeRatios[ratio_edge_id] > 0.0) {
                // Set shape function value according to the edge ratio of the extrapolated intersection.
                rNegSideCondMatrix(n_nodes+id_edge, node_ids_for_geometry[id_edge][0]) = 1.0 - mExtraEdgeRatios[ratio_edge_id];
                rNegSideCondMatrix(n_nodes+id_edge, node_ids_for_geometry[id_edge][1]) = mExtraEdgeRatios[ratio_edge_id];
            } else {
                // Set to one the shape function value along the negative side of the edge.
                rNegSideCondMatrix(n_nodes+id_edge, edge_node_i) = (nodal_distances(edge_node_i) < 0.0) ? 1.0 : 0.0;
                rNegSideCondMatrix(n_nodes+id_edge, edge_node_j) = (nodal_distances(edge_node_j) < 0.0) ? 1.0 : 0.0;
            }
        }
    }
}

}; //namespace Kratos
