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
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"

namespace Kratos
{

/// Triangle2D3AusasModifiedShapeFunctions implementation
/// Default constructor
Triangle2D3AusasModifiedShapeFunctions::Triangle2D3AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTriangleSplitter(boost::make_shared<DivideTriangle2D3>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTriangleSplitter->GenerateDivision();
    mpTriangleSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Triangle2D3AusasModifiedShapeFunctions::~Triangle2D3AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle2D3AusasModifiedShapeFunctions::Info() const {
    return "Triangle2D3N modified shape functions computation class.";
};

/// Print information about this object.
void Triangle2D3AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Triangle2D3N modified shape functions computation class.";
};

/// Print object's data.
void Triangle2D3AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Triangle2D3N modified shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        distances_buffer << std::to_string(nodal_distances(i)) << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Returns true if the element is splitting
bool Triangle2D3AusasModifiedShapeFunctions::IsSplit() {
    return mpTriangleSplitter->mIsSplit;
};

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void Triangle2D3AusasModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_pos_side;
        SetPositiveSideCondensationMatrix(p_matrix_pos_side,
                                          mpTriangleSplitter->mEdgeNodeI,
                                          mpTriangleSplitter->mEdgeNodeJ,
                                          mpTriangleSplitter->mSplitEdges);

        // Compute the positive side values
        this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                     rPositiveSideShapeFunctionsGradientsValues,
                                     rPositiveSideWeightsValues,
                                     mpTriangleSplitter->mPositiveSubdivisions,
                                     p_matrix_pos_side,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void Triangle2D3AusasModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_neg_side;
        SetNegativeSideCondensationMatrix(p_matrix_neg_side,
                                          mpTriangleSplitter->mEdgeNodeI,
                                          mpTriangleSplitter->mEdgeNodeJ,
                                          mpTriangleSplitter->mSplitEdges);

        // Compute the negative side values
        this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                     rNegativeSideShapeFunctionsGradientsValues,
                                     rNegativeSideWeightsValues,
                                     mpTriangleSplitter->mNegativeSubdivisions,
                                     p_matrix_neg_side,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void Triangle2D3AusasModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side,
                                                mpTriangleSplitter->mEdgeNodeI,
                                                mpTriangleSplitter->mEdgeNodeJ,
                                                mpTriangleSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeInterfaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
                                              rInterfacePositiveSideShapeFunctionsGradientsValues,
                                              rInterfacePositiveSideWeightsValues,
                                              mpTriangleSplitter->mPositiveInterfaces,
                                              mpTriangleSplitter->mPositiveSubdivisions,
                                              mpTriangleSplitter->mPositiveInterfacesParentIds,
                                              p_matrix_pos_side,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void Triangle2D3AusasModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side,
                                                mpTriangleSplitter->mEdgeNodeI,
                                                mpTriangleSplitter->mEdgeNodeJ,
                                                mpTriangleSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeInterfaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
                                              rInterfaceNegativeSideShapeFunctionsGradientsValues,
                                              rInterfaceNegativeSideWeightsValues,
                                              mpTriangleSplitter->mNegativeInterfaces,
                                              mpTriangleSplitter->mNegativeSubdivisions,
                                              mpTriangleSplitter->mNegativeInterfacesParentIds,
                                              p_matrix_neg_side,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void Triangle2D3AusasModifiedShapeFunctions::ComputePositiveSideInterfaceUnitNormals(
    std::vector<Vector> &rPositiveSideInterfaceUnitNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeInterfaceNormalOnOneSide(rPositiveSideInterfaceUnitNormal,
                                              mpTriangleSplitter->mPositiveInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceUnitNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void Triangle2D3AusasModifiedShapeFunctions::ComputeNegativeSideInterfaceUnitNormals(
    std::vector<Vector> &rNegativeSideInterfaceUnitNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeInterfaceNormalOnOneSide(rNegativeSideInterfaceUnitNormal,
                                              mpTriangleSplitter->mNegativeInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceUnitNormals method for a non divided geometry.";
    }
};

// Sets the condensation matrix to transform the subdivsion positive side values to entire element ones.
void Triangle2D3AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
    Matrix& rPosSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const unsigned int nedges = (this->GetInputGeometry())->EdgesNumber();
    const unsigned int nnodes = (this->GetInputGeometry())->PointsNumber();
        
    // Initialize intersection points condensation matrix
    rPosSideCondMatrix = ZeroMatrix(nnodes + nedges, nnodes);

    // Get the nodal distances vector
    const array_1d<double, 3> nodal_distances = this->GetNodalDistances();

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
void Triangle2D3AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
    Matrix& rNegSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const unsigned int nedges = (this->GetInputGeometry())->EdgesNumber();
    const unsigned int nnodes = (this->GetInputGeometry())->PointsNumber();
        
    // Initialize intersection points condensation matrix
    rNegSideCondMatrix = ZeroMatrix(nnodes + nedges, nnodes);

    // Get the nodal distances vector
    const array_1d<double, 3> nodal_distances = this->GetNodalDistances();

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
