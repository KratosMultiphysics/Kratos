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
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{

/// Tetrahedra3D4ModifiedShapeFunctions implementation
/// Default constructor
Tetrahedra3D4ModifiedShapeFunctions::Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(boost::make_shared<DivideTetrahedra3D4>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Tetrahedra3D4ModifiedShapeFunctions::~Tetrahedra3D4ModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Tetrahedra3D4ModifiedShapeFunctions::Info() const {
    return "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print information about this object.
void Tetrahedra3D4ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print object's data.
void Tetrahedra3D4ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
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
bool Tetrahedra3D4ModifiedShapeFunctions::IsSplit() {
    return mpTetrahedraSplitter->mIsSplit;
};

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void Tetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side values
        this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                     rPositiveSideShapeFunctionsGradientsValues,
                                     rPositiveSideWeightsValues,
                                     mpTetrahedraSplitter->mPositiveSubdivisions,
                                     p_matrix,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        // Compute the negative side values
        this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                     rNegativeSideShapeFunctionsGradientsValues,
                                     rNegativeSideWeightsValues,
                                     mpTetrahedraSplitter->mNegativeSubdivisions,
                                     p_matrix,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void Tetrahedra3D4ModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
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
                                              p_matrix,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void Tetrahedra3D4ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    std::vector<Matrix> &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
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
                                              p_matrix,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void Tetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    std::vector<Vector> &rPositiveSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeInterfaceNormalOnOneSide(rPositiveSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->mPositiveInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    std::vector<Vector> &rNegativeSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeInterfaceNormalOnOneSide(rNegativeSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->mNegativeInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

}; //namespace Kratos
