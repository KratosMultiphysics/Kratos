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
#include "custom_tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{

/// Tetrahedra3D4ModifiedShapeFunctions implementation
/// Default constructor
CustomTetrahedra3D4ModifiedShapeFunctions::CustomTetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
    CustomModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(Kratos::make_shared<ContactLineDivideTetrahedra3D4<Node>>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

CustomTetrahedra3D4ModifiedShapeFunctions::CustomTetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances, const std::vector<int> rStructureNodes) :
    CustomModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTetrahedraSplitter(Kratos::make_shared<ContactLineDivideTetrahedra3D4<Node>>(*pInputGeometry, rNodalDistances, rStructureNodes)) {

    // Perform the element splitting
    mpTetrahedraSplitter->GenerateDivision();
    mpTetrahedraSplitter->GenerateIntersectionsSkin();
};

/// Destructor
CustomTetrahedra3D4ModifiedShapeFunctions::~CustomTetrahedra3D4ModifiedShapeFunctions() {};

/// Turn back information as a string.
std::string CustomTetrahedra3D4ModifiedShapeFunctions::Info() const {
    return "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print information about this object.
void CustomTetrahedra3D4ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
    rOStream << "Tetrahedra3D4N modified shape functions computation class.";
};

/// Print object's data.
void CustomTetrahedra3D4ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
    const GeometryPointerType p_geometry = this->GetInputGeometry();
    const Vector nodal_distances = this->GetNodalDistances();
    rOStream << "Tetrahedra3D4N modified shape functions computation class:\n";
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
const ContactLineDivideGeometry<Node>::Pointer CustomTetrahedra3D4ModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTetrahedraSplitter;
};

// Returns true if the element is splitting
bool CustomTetrahedra3D4ModifiedShapeFunctions::IsSplit() {
    return mpTetrahedraSplitter->mIsSplit;
};

void CustomTetrahedra3D4ModifiedShapeFunctions::SetCondensationMatrix(Matrix& rIntPointCondMatrix)
{
    CustomModifiedShapeFunctions::SetCondensationMatrix(
        rIntPointCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

void CustomTetrahedra3D4ModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    CustomModifiedShapeFunctions::SetCondensationMatrix(
        rPosSideCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

void CustomTetrahedra3D4ModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    CustomModifiedShapeFunctions::SetCondensationMatrix(
        rNegSideCondMatrix,
        mpTetrahedraSplitter->mEdgeNodeI,
        mpTetrahedraSplitter->mEdgeNodeJ,
        mpTetrahedraSplitter->mSplitEdges);
}

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    // mpTetrahedraSplitter->mEdgeNodeI,
                                    // mpTetrahedraSplitter->mEdgeNodeJ,
                                    // mpTetrahedraSplitter->mSplitEdges);

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
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    // mpTetrahedraSplitter->mEdgeNodeI,
                                    // mpTetrahedraSplitter->mEdgeNodeJ,
                                    // mpTetrahedraSplitter->mSplitEdges);

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
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    // mpTetrahedraSplitter->mEdgeNodeI,
                                    // mpTetrahedraSplitter->mEdgeNodeJ,
                                    // mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
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
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    // mpTetrahedraSplitter->mEdgeNodeI,
                                    // mpTetrahedraSplitter->mEdgeNodeJ,
                                    // mpTetrahedraSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
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

// Returns all the shape function values for the contact line.
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
    std::vector<unsigned int>& ContactLineIndices,
    std::vector<Matrix> &rContactLineNegativeSideShapeFunctionsValues,
    std::vector<ShapeFunctionsGradientsType> &rContactLineNegativeSideShapeFunctionsGradientsValues,
    std::vector<Vector> &rContactLineNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod){

    rContactLineNegativeSideShapeFunctionsValues.clear();
    rContactLineNegativeSideShapeFunctionsGradientsValues.clear();
    rContactLineNegativeSideWeightsValues.clear();

    if (ContactLineIndices.size() > 0){
        std::vector < unsigned int >& contact_interface_ids = mpTetrahedraSplitter->mContactInterface;

        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    // mpTetrahedraSplitter->mEdgeNodeI,
                                    // mpTetrahedraSplitter->mEdgeNodeJ,
                                    // mpTetrahedraSplitter->mSplitEdges);

        for (unsigned int i_cl = 0; i_cl < ContactLineIndices.size(); i_cl++){
            const unsigned int cl_id = ContactLineIndices[i_cl];
            const unsigned int i_parent = 
                    mpTetrahedraSplitter->mNegativeInterfacesParentIds[(contact_interface_ids[cl_id])];
            const auto& r_parent_geom = mpTetrahedraSplitter->mNegativeSubdivisions[i_parent];
            const auto& r_interface_geom = (mpTetrahedraSplitter->mContactLine)[cl_id];

            Matrix contact_line_negative_side_shape_function_values;
            ShapeFunctionsGradientsType contact_line_negative_side_shape_function_gradient_values;
            Vector contact_line_negative_side_weight_values;

            // Compute the positive side interface values
            this->ComputeFaceValuesOnOneSide(contact_line_negative_side_shape_function_values,
                                        contact_line_negative_side_shape_function_gradient_values,
                                        contact_line_negative_side_weight_values,
                                        r_interface_geom,
                                        r_parent_geom,
                                        p_matrix,
                                        IntegrationMethod);

            rContactLineNegativeSideShapeFunctionsValues.push_back(contact_line_negative_side_shape_function_values);
            rContactLineNegativeSideShapeFunctionsGradientsValues.push_back(contact_line_negative_side_shape_function_gradient_values);
            rContactLineNegativeSideWeightsValues.push_back(contact_line_negative_side_weight_values);
        }
        
    } 

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/*     const int i_contact_face = mpTetrahedraSplitter->mContactInterface;

    if (i_contact_face > -1) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix,
                                    mpTetrahedraSplitter->mEdgeNodeI,
                                    mpTetrahedraSplitter->mEdgeNodeJ,
                                    mpTetrahedraSplitter->mSplitEdges);

        const unsigned int i_parent = mpTetrahedraSplitter->mNegativeInterfacesParentIds[mpTetrahedraSplitter->mContactInterface];
        const auto& r_interface_geom = mpTetrahedraSplitter->mContactLine;
        const auto& r_parent_geom = mpTetrahedraSplitter->mNegativeSubdivisions[i_parent];

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rContactLineNegativeSideShapeFunctionsValues,
                                        rContactLineNegativeSideShapeFunctionsGradientsValues,
                                        rContactLineNegativeSideWeightsValues,
                                        r_interface_geom,
                                        r_parent_geom,
                                        p_matrix,
                                        IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues method for a geometry without contact line.";
    } */
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

};

// Given a face id, computes the positive side subdivision shape function values in that face.
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rPositiveExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
    Vector &rPositiveExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
            // mpTetrahedraSplitter->mEdgeNodeI,
            // mpTetrahedraSplitter->mEdgeNodeJ,
            // mpTetrahedraSplitter->mSplitEdges);
        
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
            p_matrix,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rNegativeExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
    Vector &rNegativeExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
            // mpTetrahedraSplitter->mEdgeNodeI,
            // mpTetrahedraSplitter->mEdgeNodeJ,
            // mpTetrahedraSplitter->mSplitEdges);
        
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
            p_matrix,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards unit normal vector values.
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    AreaNormalsContainerType &rPositiveSideInterfaceAreaNormals,
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

// Compute the positive side interface outwards unit normal vector values.
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    AreaNormalsContainerType &rNegativeSideInterfaceAreaNormals,
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
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputePositiveExteriorFaceAreaNormals(
    AreaNormalsContainerType &rPositiveExteriorFaceAreaNormal,
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
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeExteriorFaceAreaNormals(
    AreaNormalsContainerType &rNegativeExteriorFaceAreaNormal,
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

// Computes the positive side shape function values in the edges intersections
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeShapeFunctionsOnPositiveEdgeIntersections(
    Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues){

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
            // mpTetrahedraSplitter->mEdgeNodeI,
            // mpTetrahedraSplitter->mEdgeNodeJ,
            // mpTetrahedraSplitter->mSplitEdges);

        // Compute the edge intersections shape function values
        this->ComputeEdgeIntersectionValuesOnOneSide(
            p_matrix,
            rPositiveEdgeIntersectionsShapeFunctionsValues);

    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnPositiveEdgeIntersections method for a non divided geometry.";
    }
};

// Computes the negative side shape function values in the edges intersections
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(
    Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues){
    
    if (this->IsSplit()) {
        // Note that positive and negative sides values are equal for standard shape functions
        this->ComputeShapeFunctionsOnPositiveEdgeIntersections(
            rNegativeEdgeIntersectionsShapeFunctionsValues);
    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnNegativeEdgeIntersections method for a non divided geometry.";
    }
};

/* bool Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideContactLineVector(
    Vector &rNegativeSideContactLineVector) */
void CustomTetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideContactLineVector(
    std::vector<unsigned int>& FaceIndices,
    std::vector<Vector> &rNegativeSideContactLineVector)
{
    FaceIndices.clear();
    rNegativeSideContactLineVector.clear();

    FaceIndices = mpTetrahedraSplitter->mContactFace;

    //KRATOS_INFO("ComputeNegativeSideContactLineVector") << (mpTetrahedraSplitter->mContactFace).size() << std::endl;

    if (FaceIndices.size() > 0){
        
        std::vector < unsigned int >& contact_interface_ids = mpTetrahedraSplitter->mContactInterface;
        std::vector < unsigned int >& contact_interface_edge_ids = mpTetrahedraSplitter->mContactEdge;

        for (unsigned int i_cl = 0; i_cl < FaceIndices.size(); i_cl ++){
            const unsigned int i_interface = contact_interface_ids[i_cl];
            const unsigned int i_edge = contact_interface_edge_ids[i_cl];

            //KRATOS_INFO("(contact_interface_ids[i_cl])") << i_interface << std::endl;
            //KRATOS_INFO("(contact_interface_edge_ids[i_cl])") << i_edge << std::endl;

            const auto& edges = (mpTetrahedraSplitter->mNegativeInterfaces)[i_interface]->Edges();
            const auto& edgei = edges[i_edge];

            //KRATOS_INFO("edgei[0].Coordinates()") << edgei[0].Coordinates() << std::endl;
            //KRATOS_INFO("edgei[1].Coordinates()") << edgei[1].Coordinates() << std::endl;

            Vector aux_vector = edgei[1].Coordinates() - edgei[0].Coordinates();

            //KRATOS_INFO("aux_vector") << aux_vector << std::endl;

            const double norm = Kratos::norm_2(aux_vector);
            aux_vector = 1.0/norm * aux_vector;

            rNegativeSideContactLineVector.push_back(aux_vector);
        }
    }

    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    /* const int i_contact_face = mpTetrahedraSplitter->mContactInterface;
    const int i_contact_edge = mpTetrahedraSplitter->mContactEdge;

    if (i_contact_face > -1){
        IndexedPointGeometryType edgei = (mpTetrahedraSplitter->mNegativeInterfaces[i_contact_face]->Edges())[i_contact_edge];
        const array_1d<double, 3> aux_vector = 
            edgei[1].Coordinates() - edgei[0].Coordinates();

        rNegativeSideContactLineVector = aux_vector;
        return true;
    }
    else
    {
        rNegativeSideContactLineVector = ZeroVector(3);
        return false;
    } */
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    
};

}; //namespace Kratos
