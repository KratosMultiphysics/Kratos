//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "modified_shape_functions/ausas_modified_shape_functions.h"

namespace Kratos
{

/// AusasModifiedShapeFunctions implementation
/// Default constructor
AusasModifiedShapeFunctions::AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances) {
};

/// Destructor
AusasModifiedShapeFunctions::~AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string AusasModifiedShapeFunctions::Info() const {
    return "Ausas modified shape functions computation base class.";
};

/// Print information about this object.
void AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Ausas modified shape functions computation base class.";
};

/// Print object's data.
void AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Ausas modified shape functions computation base class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::ostringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void AusasModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod)
{
    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side);

        // Compute the positive side values
        const auto& r_splitting_util = *(pGetSplittingUtil());
        this->ComputeValuesOnOneSide(
            rPositiveSideShapeFunctionsValues,
            rPositiveSideShapeFunctionsGradientsValues,
            rPositiveSideWeightsValues,
            r_splitting_util.GetPositiveSubdivisions(),
            p_matrix_pos_side,
            IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the \'ComputePositiveSideShapeFunctionsAndGradientsValues\' method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void AusasModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod)
{
    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side);

        // Compute the negative side values
        const auto& r_splitting_util = *(pGetSplittingUtil());
        this->ComputeValuesOnOneSide(
            rNegativeSideShapeFunctionsValues,
            rNegativeSideShapeFunctionsGradientsValues,
            rNegativeSideWeightsValues,
            r_splitting_util.GetNegativeSubdivisions(),
            p_matrix_neg_side,
            IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the \'ComputeNegativeSideShapeFunctionsAndGradientsValues\' method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void AusasModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod)
{
    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side);

        // Compute the positive side interface values
        const auto& r_splitting_util = *(pGetSplittingUtil());
        this->ComputeFaceValuesOnOneSide(
            rInterfacePositiveSideShapeFunctionsValues,
            rInterfacePositiveSideShapeFunctionsGradientsValues,
            rInterfacePositiveSideWeightsValues,
            r_splitting_util.GetPositiveInterfaces(),
            r_splitting_util.GetPositiveSubdivisions(),
            r_splitting_util.GetPositiveInterfacesParentIds(),
            p_matrix_pos_side,
            IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the \'ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues\' method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void AusasModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod)
{
    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side);

        // Compute the positive side interface values
        const auto& r_splitting_util = *(pGetSplittingUtil());
        this->ComputeFaceValuesOnOneSide(
            rInterfaceNegativeSideShapeFunctionsValues,
            rInterfaceNegativeSideShapeFunctionsGradientsValues,
            rInterfaceNegativeSideWeightsValues,
            r_splitting_util.GetNegativeInterfaces(),
            r_splitting_util.GetNegativeSubdivisions(),
            r_splitting_util.GetNegativeInterfacesParentIds(),
            p_matrix_neg_side,
            IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the \'ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues\' method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void AusasModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rPositiveExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
    Vector &rPositiveExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod)
{
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side);

        // Get the external faces
        auto& r_splitting_util = *(pGetSplittingUtil());
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        r_splitting_util.GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            r_splitting_util.GetPositiveSubdivisions(),
            FaceId);

        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rPositiveExteriorFaceShapeFunctionsValues,
            rPositiveExteriorFaceShapeFunctionsGradientsValues,
            rPositiveExteriorFaceWeightsValues,
            exterior_faces_vector,
            r_splitting_util.GetPositiveSubdivisions(),
            exterior_faces_parent_ids_vector,
            p_matrix_pos_side,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the \'ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues\' method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void AusasModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rNegativeExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
    Vector &rNegativeExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod)
{
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side);

        // Get the external faces
        auto& r_splitting_util = *(pGetSplittingUtil());
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        r_splitting_util.GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            r_splitting_util.GetNegativeSubdivisions(),
            FaceId);

        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rNegativeExteriorFaceShapeFunctionsValues,
            rNegativeExteriorFaceShapeFunctionsGradientsValues,
            rNegativeExteriorFaceWeightsValues,
            exterior_faces_vector,
            r_splitting_util.GetNegativeSubdivisions(),
            exterior_faces_parent_ids_vector,
            p_matrix_neg_side,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the \'ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues\' method for a non divided geometry.";
    }
};

// Sets the condensation matrix to transform the subdivsion positive side values to entire element ones.
void AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
    Matrix& rPosSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const unsigned int nedges = (this->GetInputGeometry())->EdgesNumber();
    const unsigned int nnodes = (this->GetInputGeometry())->PointsNumber();

    // Initialize intersection points condensation matrix
    rPosSideCondMatrix = ZeroMatrix(nnodes + nedges, nnodes);

    // Get the nodal distances vector
    const Vector& nodal_distances = this->GetNodalDistances();

    // Fill the original geometry points main diagonal
    for (unsigned int i = 0; i < nnodes; ++i) {
        rPosSideCondMatrix(i,i) = (nodal_distances(i) > 0.0) ? 1.0 : 0.0;
    }

    // Compute the intersection points contributions
    unsigned int row = nnodes;
    for (unsigned int idedge = 0; idedge < nedges; ++idedge) {
        // Check if the edge has an intersection point
        if (rSplitEdges[nnodes+idedge] != -1) {
            // Get the nodes that compose the edge
            const unsigned int edge_node_i = rEdgeNodeI[idedge];
            const unsigned int edge_node_j = rEdgeNodeJ[idedge];

            // Set to one the shape function value along the positive side of the edge.
            rPosSideCondMatrix(row, edge_node_i) = (nodal_distances(edge_node_i) > 0.0) ? 1.0 : 0.0;
            rPosSideCondMatrix(row, edge_node_j) = (nodal_distances(edge_node_j) > 0.0) ? 1.0 : 0.0;
        }
        row++;
    }
}

// Sets the condensation matrix to transform the subdivsion negative side values to entire element ones.
void AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
    Matrix& rNegSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const unsigned int nedges = (this->GetInputGeometry())->EdgesNumber();
    const unsigned int nnodes = (this->GetInputGeometry())->PointsNumber();

    // Initialize intersection points condensation matrix
    rNegSideCondMatrix = ZeroMatrix(nnodes + nedges, nnodes);

    // Get the nodal distances vector
    const Vector& nodal_distances = this->GetNodalDistances();

    // Fill the original geometry points main diagonal
    for (unsigned int i = 0; i < nnodes; ++i) {
        rNegSideCondMatrix(i,i) = (nodal_distances(i) < 0.0) ? 1.0 : 0.0;
    }

    // Compute the intersection points contributions
    unsigned int row = nnodes;
    for (unsigned int idedge = 0; idedge < nedges; ++idedge) {
        // Check if the edge has an intersection point
        if (rSplitEdges[nnodes+idedge] != -1) {
            // Get the nodes that compose the edge
            const unsigned int edge_node_i = rEdgeNodeI[idedge];
            const unsigned int edge_node_j = rEdgeNodeJ[idedge];

            // Set to one the shape function value along the negative side of the edge.
            rNegSideCondMatrix(row, edge_node_i) = (nodal_distances(edge_node_i) < 0.0) ? 1.0 : 0.0;
            rNegSideCondMatrix(row, edge_node_j) = (nodal_distances(edge_node_j) < 0.0) ? 1.0 : 0.0;
        }
        row++;
    }
}

}; //namespace Kratos
