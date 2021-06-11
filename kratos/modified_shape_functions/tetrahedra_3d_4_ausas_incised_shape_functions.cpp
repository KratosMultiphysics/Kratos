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
#include "modified_shape_functions/tetrahedra_3d_4_ausas_incised_shape_functions.h"

namespace Kratos
{

/// Tetrahedra3D4AusasIncisedShapeFunctions implementation
/// Default constructor
Tetrahedra3D4AusasIncisedShapeFunctions::Tetrahedra3D4AusasIncisedShapeFunctions(
    const GeometryPointerType pInputGeometry,
    const Vector& rNodalDistancesWithExrapolated,
    const Vector& rExtrapolatedEdgeRatios)
    : Tetrahedra3D4AusasModifiedShapeFunctions(pInputGeometry, rNodalDistancesWithExrapolated)
    , mExtraEdgeRatios(rExtrapolatedEdgeRatios)
{};

/// Destructor
Tetrahedra3D4AusasIncisedShapeFunctions::~Tetrahedra3D4AusasIncisedShapeFunctions() {};

/// Turn back information as a string.
std::string Tetrahedra3D4AusasIncisedShapeFunctions::Info() const
{
    return "Tetrahedra3D4N Ausas incised shape functions computation class.";
}

/// Print information about this object.
void Tetrahedra3D4AusasIncisedShapeFunctions::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "Tetrahedra3D4N Ausas incised shape functions computation class.";
}

/// Print object's data.
void Tetrahedra3D4AusasIncisedShapeFunctions::PrintData(std::ostream& rOStream) const
{
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    const Vector extra_edge_ratios = this->GetExtrapolatedEdgeRatios();

    rOStream << "Tetrahedra3D4N Ausas incised shape functions computation class:\n";
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
const Vector& Tetrahedra3D4AusasIncisedShapeFunctions::GetExtrapolatedEdgeRatios() const
{
    return mExtraEdgeRatios;
}

void Tetrahedra3D4AusasIncisedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    Tetrahedra3D4AusasIncisedShapeFunctions::SetPositiveSideCondensationMatrix(
        rPosSideCondMatrix,
        Tetrahedra3D4AusasModifiedShapeFunctions::mpTetrahedraSplitter->mEdgeNodeI,
        Tetrahedra3D4AusasModifiedShapeFunctions::mpTetrahedraSplitter->mEdgeNodeJ,
        Tetrahedra3D4AusasModifiedShapeFunctions::mpTetrahedraSplitter->mSplitEdges);
}

void Tetrahedra3D4AusasIncisedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    Tetrahedra3D4AusasIncisedShapeFunctions::SetNegativeSideCondensationMatrix(
        rNegSideCondMatrix,
        Tetrahedra3D4AusasModifiedShapeFunctions::mpTetrahedraSplitter->mEdgeNodeI,
        Tetrahedra3D4AusasModifiedShapeFunctions::mpTetrahedraSplitter->mEdgeNodeJ,
        Tetrahedra3D4AusasModifiedShapeFunctions::mpTetrahedraSplitter->mSplitEdges);
}

// Sets the condensation matrix to transform the subdivsion positive side values to entire element ones.
void Tetrahedra3D4AusasIncisedShapeFunctions::SetPositiveSideCondensationMatrix(
    Matrix& rPosSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges)
{
    const std::size_t n_nodes = 4;
    const std::size_t n_edges = 6;

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

            // Transform definition of edge ID from divide_tetrahedra_3d_4.h to definition of tetrahedra_3d_4.h (geometry)
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
void Tetrahedra3D4AusasIncisedShapeFunctions::SetNegativeSideCondensationMatrix(
    Matrix& rNegSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges)
{
    const std::size_t n_nodes = 4;
    const std::size_t n_edges = 6;

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

            // Transform definition of edge ID from divide_tetrahedra_3d_4.h to definition of tetrahedra_3d_4.h (geometry)
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
