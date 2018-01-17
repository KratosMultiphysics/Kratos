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
    AusasModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Tetrahedra3D4AusasModifiedShapeFunctions::~Tetrahedra3D4AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Tetrahedra3D4AusasModifiedShapeFunctions::Info() const {
    return "Tetrahedra3D4N Ausas modified shape functions computation class.";
};

/// Print information about this object.
void Tetrahedra3D4AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N Ausas modified shape functions computation class.";
};

/// Print object's data.
void Tetrahedra3D4AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Tetrahedra3D4N Ausas modified shape functions computation class:\n";
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
    ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
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
    ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
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
    ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side,
                                                mpTetrahedraSplitter->mEdgeNodeI,
                                                mpTetrahedraSplitter->mEdgeNodeJ,
                                                mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
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
    ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
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
        this->ComputeFaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
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

// Given a face id, computes the positive side subdivision shape function values in that face.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rPositiveExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
    Vector &rPositiveExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(
            p_matrix_pos_side,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);
        
        // Get the external faces
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mPositiveSubdivisions,
            FaceId);
        
        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rPositiveExteriorFaceShapeFunctionsValues,
            rPositiveExteriorFaceShapeFunctionsGradientsValues,
            rPositiveExteriorFaceWeightsValues,
            exterior_faces_vector,
            mpTetrahedraSplitter->mPositiveSubdivisions,
            exterior_faces_parent_ids_vector,
            p_matrix_pos_side,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rNegativeExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
    Vector &rNegativeExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(
            p_matrix_neg_side,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);
        
        // Get the external faces
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mNegativeSubdivisions,
            FaceId);
        
        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rNegativeExteriorFaceShapeFunctionsValues,
            rNegativeExteriorFaceShapeFunctionsGradientsValues,
            rNegativeExteriorFaceWeightsValues,
            exterior_faces_vector,
            mpTetrahedraSplitter->mNegativeSubdivisions,
            exterior_faces_parent_ids_vector,
            p_matrix_neg_side,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    std::vector<Vector> &rPositiveSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeFaceNormalOnOneSide(rPositiveSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->mPositiveInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    std::vector<Vector> &rNegativeSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeFaceNormalOnOneSide(rNegativeSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->mNegativeInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side exteriorface outwards area normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputePositiveExteriorFaceAreaNormals(
    std::vector<Vector> &rPositiveExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod)
{

    if (this->IsSplit())
    {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mPositiveSubdivisions,
            FaceId);

        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rPositiveExteriorFaceAreaNormal,
                                              exterior_faces_vector,
                                              IntegrationMethod);
    }
    else
    {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side exterior face outwards area normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeNegativeExteriorFaceAreaNormals(
    std::vector<Vector> &rNegativeExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod)
{

    if (this->IsSplit())
    {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTetrahedraSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTetrahedraSplitter->mNegativeSubdivisions,
            FaceId);

        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rNegativeExteriorFaceAreaNormal,
                                              exterior_faces_vector,
                                              IntegrationMethod);
    }
    else
    {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceAreaNormals method for a non divided geometry.";
    }
};

}; //namespace Kratos
