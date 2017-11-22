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
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{

/// Triangle2D3ModifiedShapeFunctions implementation
/// Default constructor
Triangle2D3ModifiedShapeFunctions::Triangle2D3ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTriangleSplitter(boost::make_shared<DivideTriangle2D3>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTriangleSplitter->GenerateDivision();
    mpTriangleSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Triangle2D3ModifiedShapeFunctions::~Triangle2D3ModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle2D3ModifiedShapeFunctions::Info() const {
    return "Triangle2D3N modified shape functions computation class.";
};

/// Print information about this object.
void Triangle2D3ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Triangle2D3N modified shape functions computation class.";
};

/// Print object's data.
void Triangle2D3ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
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
bool Triangle2D3ModifiedShapeFunctions::IsSplit() {
    return mpTriangleSplitter->mIsSplit;
};

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void Triangle2D3ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        SetCondensationMatrix(p_matrix,
                              mpTriangleSplitter->mEdgeNodeI,
                              mpTriangleSplitter->mEdgeNodeJ,
                              mpTriangleSplitter->mSplitEdges);

        // Compute the positive side values
        this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                     rPositiveSideShapeFunctionsGradientsValues,
                                     rPositiveSideWeightsValues,
                                     mpTriangleSplitter->mPositiveSubdivisions,
                                     p_matrix,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void Triangle2D3ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        SetCondensationMatrix(p_matrix,
                              mpTriangleSplitter->mEdgeNodeI,
                              mpTriangleSplitter->mEdgeNodeJ,
                              mpTriangleSplitter->mSplitEdges);

        // Compute the negative side values
        this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                     rNegativeSideShapeFunctionsGradientsValues,
                                     rNegativeSideWeightsValues,
                                     mpTriangleSplitter->mNegativeSubdivisions,
                                     p_matrix,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void Triangle2D3ModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
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
                                              p_matrix,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void Triangle2D3ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
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
                                              p_matrix,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Triangle2D3ModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    std::vector<Vector> &rPositiveSideInterfaceAreaNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars area normal values
        this->ComputeInterfaceNormalOnOneSide(rPositiveSideInterfaceAreaNormal,
                                              mpTriangleSplitter->mPositiveInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Triangle2D3ModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    std::vector<Vector> &rNegativeSideInterfaceAreaNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars area normal values
        this->ComputeInterfaceNormalOnOneSide(rNegativeSideInterfaceAreaNormal,
                                              mpTriangleSplitter->mNegativeInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

}; //namespace Kratos
