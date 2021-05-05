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
    AusasModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTriangleSplitter(Kratos::make_shared<DivideTriangle2D3>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTriangleSplitter->GenerateDivision();
    mpTriangleSplitter->GenerateIntersectionsSkin();
};

/// Destructor
Triangle2D3AusasModifiedShapeFunctions::~Triangle2D3AusasModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string Triangle2D3AusasModifiedShapeFunctions::Info() const {
    return "Triangle2D3N Ausas modified shape functions computation class.";
};

/// Print information about this object.
void Triangle2D3AusasModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Triangle2D3N Ausas modified shape functions computation class.";
};

/// Print object's data.
void Triangle2D3AusasModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Triangle2D3N Ausas modified shape functions computation class:\n";
    rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
    std::stringstream distances_buffer;
    std::ostringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Returns a pointer to the splitting utility
const DivideGeometry::Pointer Triangle2D3AusasModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTriangleSplitter;
};

void Triangle2D3AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetPositiveSideCondensationMatrix(
        rPosSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle2D3AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    AusasModifiedShapeFunctions::SetNegativeSideCondensationMatrix(
        rNegSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

// Compute the positive side interface outwards area normal vector values.
void Triangle2D3AusasModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    std::vector<Vector> &rPositiveSideInterfaceAreaNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rPositiveSideInterfaceAreaNormal,
                                              mpTriangleSplitter->GetPositiveInterfaces(),
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Triangle2D3AusasModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    std::vector<Vector> &rNegativeSideInterfaceAreaNormal,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rNegativeSideInterfaceAreaNormal,
                                              mpTriangleSplitter->GetNegativeInterfaces(),
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side exteriorface outwards area normal vector values.
void Triangle2D3AusasModifiedShapeFunctions::ComputePositiveExteriorFaceAreaNormals(
    std::vector<Vector> &rPositiveExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod)
{

    if (this->IsSplit())
    {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTriangleSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTriangleSplitter->GetPositiveSubdivisions(),
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
void Triangle2D3AusasModifiedShapeFunctions::ComputeNegativeExteriorFaceAreaNormals(
    std::vector<Vector> &rNegativeExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod)
{

    if (this->IsSplit())
    {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTriangleSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTriangleSplitter->GetNegativeSubdivisions(),
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
void Triangle2D3AusasModifiedShapeFunctions::ComputeShapeFunctionsOnPositiveEdgeIntersections(
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
void Triangle2D3AusasModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(
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
