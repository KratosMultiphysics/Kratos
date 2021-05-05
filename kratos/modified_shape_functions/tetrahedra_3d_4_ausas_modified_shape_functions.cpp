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
    std::stringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};


// Returns a pointer to the splitting utility
const DivideGeometry::Pointer Tetrahedra3D4AusasModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTetrahedraSplitter;
};

// Returns true if the element is splitting
bool Tetrahedra3D4AusasModifiedShapeFunctions::IsSplit() const
{
    return mpTetrahedraSplitter->mIsSplit;
};

void Tetrahedra3D4AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
        rPosSideCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

void Tetrahedra3D4AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
        rNegSideCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

// Compute the positive side interface outwards area normal vector values.
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    std::vector<Vector> &rPositiveSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars unit normal values
        this->ComputeFaceNormalOnOneSide(rPositiveSideInterfaceAreaNormals,
                                              mpTetrahedraSplitter->GetPositiveInterfaces(),
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
                                              mpTetrahedraSplitter->GetNegativeInterfaces(),
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
            mpTetrahedraSplitter->GetPositiveSubdivisions(),
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
            mpTetrahedraSplitter->GetNegativeSubdivisions(),
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

// Computes the positive side shape function values in the edges intersections
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeShapeFunctionsOnPositiveEdgeIntersections(
    Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues){

    if (this->IsSplit()) {
        // Get the positive side condensation matrix
         Matrix p_matrix_pos_side;
        this->SetPositiveSideCondensationMatrix(p_matrix_pos_side);

        // Compute the edge intersections shape function values
        this->ComputeEdgeIntersectionValuesOnOneSide(
            p_matrix_pos_side,
            rPositiveEdgeIntersectionsShapeFunctionsValues);

    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnPositiveEdgeIntersections method for a non divided geometry.";
    }
};

// Computes the negative side shape function values in the edges intersections
void Tetrahedra3D4AusasModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(
    Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues){

    if (this->IsSplit()) {
        // Get the positive side condensation matrix
        Matrix p_matrix_neg_side;
        this->SetNegativeSideCondensationMatrix(p_matrix_neg_side);

        // Compute the edge intersections shape function values
        this->ComputeEdgeIntersectionValuesOnOneSide(
            p_matrix_neg_side,
            rNegativeEdgeIntersectionsShapeFunctionsValues);

    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnNegativeEdgeIntersections method for a non divided geometry.";
    }
};

}; //namespace Kratos
