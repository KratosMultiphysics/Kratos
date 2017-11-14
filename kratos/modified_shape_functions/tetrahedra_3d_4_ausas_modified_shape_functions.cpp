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
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

namespace Kratos
{

/// Tetrahedra3D4AusasModifiedShapeFunctions implementation
/// Default constructor
Tetrahedra3D4AusasModifiedShapeFunctions::Tetrahedra3D4AusasModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(boost::make_shared<DivideTetrahedra3D4>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Tetrahedra3D4AusasModifiedShapeFunctions::~Tetrahedra3D4AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Tetrahedra3D4AusasModifiedShapeFunctions::Info() const {
    return "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print information about this object.
void Tetrahedra3D4AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print object's data.
void Tetrahedra3D4AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Tetrahedra3D4N modified shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        distances_buffer << std::to_string(nodal_distances(i)) << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Returns true if the element is splitting
bool Tetrahedra3D4AusasModifiedShapeFunctions::IsSplit() {
    return mpTetrahedraSplitter->mIsSplit;
};

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side,
                                                mpTetrahedraSplitter->mEdgeNodeI,
                                                mpTetrahedraSplitter->mEdgeNodeJ,
                                                mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side values
        this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                     rPositiveSideShapeFunctionsGradientsValues,
                                     rPositiveSideWeightsValues,
                                     mpTetrahedraSplitter->mPositiveSubdivisions,
                                     p_matrix_pos_side,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side,
                                                mpTetrahedraSplitter->mEdgeNodeI,
                                                mpTetrahedraSplitter->mEdgeNodeJ,
                                                mpTetrahedraSplitter->mSplitEdges);

        // Compute the negative side values
        this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                     rNegativeSideShapeFunctionsGradientsValues,
                                     rNegativeSideWeightsValues,
                                     mpTetrahedraSplitter->mNegativeSubdivisions,
                                     p_matrix_neg_side,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side,
                                                mpTetrahedraSplitter->mEdgeNodeI,
                                                mpTetrahedraSplitter->mEdgeNodeJ,
                                                mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeInterfaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
                                              rInterfacePositiveSideShapeFunctionsGradientsValues,
                                              rInterfacePositiveSideWeightsValues,
                                              mpTetrahedraSplitter->mPositiveInterfaces,
                                              mpTetrahedraSplitter->mPositiveSubdivisions,
                                              mpTetrahedraSplitter->mPositiveInterfacesParentIds,
                                              p_matrix_pos_side,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side,
                                                mpTetrahedraSplitter->mEdgeNodeI,
                                                mpTetrahedraSplitter->mEdgeNodeJ,
                                                mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeInterfaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
                                              rInterfaceNegativeSideShapeFunctionsGradientsValues,
                                              rInterfaceNegativeSideWeightsValues,
                                              mpTetrahedraSplitter->mNegativeInterfaces,
                                              mpTetrahedraSplitter->mNegativeSubdivisions,
                                              mpTetrahedraSplitter->mNegativeInterfacesParentIds,
                                              p_matrix_neg_side,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputePositiveSideInterfaceUnitNormals(
    std::vector<Vector> &rPositiveSideInterfaceUnitNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeInterfaceNormalOnOneSide(rPositiveSideInterfaceUnitNormal,
                                              mpTetrahedraSplitter->mPositiveInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceUnitNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeNegativeSideInterfaceUnitNormals(
    std::vector<Vector> &rNegativeSideInterfaceUnitNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeInterfaceNormalOnOneSide(rNegativeSideInterfaceUnitNormal,
                                              mpTetrahedraSplitter->mNegativeInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceUnitNormals method for a non divided geometry.";
    }
};

// Sets the condensation matrix to transform the subdivsion positive side values to entire element ones.
void Tetrahedra3D4AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
    Matrix& rPosSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const unsigned int nedges = (this->GetInputGeometry())->EdgesNumber();
    const unsigned int nnodes = (this->GetInputGeometry())->PointsNumber();
        
    // Initialize intersection points condensation matrix
    rPosSideCondMatrix = ZeroMatrix(nnodes + nedges, nnodes);

    // Get the nodal distances vector
    const array_1d<double, 4> nodal_distances = this->GetNodalDistances();

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
void Tetrahedra3D4AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
    Matrix& rNegSideCondMatrix,
    const std::vector<int>& rEdgeNodeI,
    const std::vector<int>& rEdgeNodeJ,
    const std::vector<int>& rSplitEdges) {

    const unsigned int nedges = (this->GetInputGeometry())->EdgesNumber();
    const unsigned int nnodes = (this->GetInputGeometry())->PointsNumber();
        
    // Initialize intersection points condensation matrix
    rNegSideCondMatrix = ZeroMatrix(nnodes + nedges, nnodes);

    // Get the nodal distances vector
    const array_1d<double, 4> nodal_distances = this->GetNodalDistances();

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
