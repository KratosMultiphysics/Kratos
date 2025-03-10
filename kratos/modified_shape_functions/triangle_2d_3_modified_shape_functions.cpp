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
    mpTriangleSplitter(Kratos::make_shared<DivideTriangle2D3<Node>>(*pInputGeometry, rNodalDistances)) {

    // Perform the element splitting
    mpTriangleSplitter->GenerateDivision();
    mpTriangleSplitter->GenerateIntersectionsSkin();
};

    //
Triangle2D3ModifiedShapeFunctions::Triangle2D3ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector &rNodalDistances, const Vector& rStructureNodes) : 
    ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
    mpTriangleSplitter(Kratos::make_shared<DivideTriangle2D3<Node>>(*pInputGeometry, rNodalDistances, rStructureNodes))
    {

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
    std::ostringstream stm;
    for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
        stm << nodal_distances(i);
        distances_buffer << stm.str() << " ";
    }
    rOStream << "\tDistance values: " << distances_buffer.str();
};

// Returns a pointer to the splitting utility
const DivideGeometry<Node>::Pointer Triangle2D3ModifiedShapeFunctions::pGetSplittingUtil() const {
    return mpTriangleSplitter;
};

// Returns true if the element is splitting
bool Triangle2D3ModifiedShapeFunctions::IsSplit() {
    return mpTriangleSplitter->mIsSplit;
};

void Triangle2D3ModifiedShapeFunctions::SetCondensationMatrix(Matrix& rIntPointCondMatrix)
{
    ModifiedShapeFunctions::SetCondensationMatrix(
        rIntPointCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle2D3ModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix& rPosSideCondMatrix)
{
    ModifiedShapeFunctions::SetCondensationMatrix(
        rPosSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

void Triangle2D3ModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix& rNegSideCondMatrix)
{
    ModifiedShapeFunctions::SetCondensationMatrix(
        rNegSideCondMatrix,
        mpTriangleSplitter->mEdgeNodeI,
        mpTriangleSplitter->mEdgeNodeJ,
        mpTriangleSplitter->mSplitEdges);
}

// Internally computes the splitting pattern and returns all the shape function values for the positive side.
void Triangle2D3ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rPositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
    Vector &rPositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        SetCondensationMatrix(p_matrix);//,
                              //mpTriangleSplitter->mEdgeNodeI,
                              //mpTriangleSplitter->mEdgeNodeJ,
                              //mpTriangleSplitter->mSplitEdges);

        // Compute the positive side values
        this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                     rPositiveSideShapeFunctionsGradientsValues,
                                     rPositiveSideWeightsValues,
                                     mpTriangleSplitter->GetPositiveSubdivisions(),//mPositiveSubdivisions,
                                     p_matrix,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative side.
void Triangle2D3ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
    Vector &rNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the intersection points condensation matrix
        Matrix p_matrix;
        SetCondensationMatrix(p_matrix);//,
                              //mpTriangleSplitter->mEdgeNodeI,
                              //mpTriangleSplitter->mEdgeNodeJ,
                              //mpTriangleSplitter->mSplitEdges);

        // Compute the negative side values
        this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                     rNegativeSideShapeFunctionsGradientsValues,
                                     rNegativeSideWeightsValues,
                                     mpTriangleSplitter->GetNegativeSubdivisions(),//mNegativeSubdivisions,
                                     p_matrix,
                                     IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
void Triangle2D3ModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfacePositiveSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
    Vector &rInterfacePositiveSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    //mpTriangleSplitter->mEdgeNodeI,
                                    //mpTriangleSplitter->mEdgeNodeJ,
                                    //mpTriangleSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
                                              rInterfacePositiveSideShapeFunctionsGradientsValues,
                                              rInterfacePositiveSideWeightsValues,
                                              mpTriangleSplitter->GetPositiveInterfaces(),//mPositiveInterfaces,
                                              mpTriangleSplitter->GetPositiveSubdivisions(),//mPositiveSubdivisions,
                                              mpTriangleSplitter->GetPositiveInterfacesParentIds(),//mPositiveInterfacesParentIds,
                                              p_matrix,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
void Triangle2D3ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
    Matrix &rInterfaceNegativeSideShapeFunctionsValues,
    ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
    Vector &rInterfaceNegativeSideWeightsValues,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(p_matrix);//,
                                    //mpTriangleSplitter->mEdgeNodeI,
                                    //mpTriangleSplitter->mEdgeNodeJ,
                                    //mpTriangleSplitter->mSplitEdges);

        // Compute the positive side interface values
        this->ComputeFaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
                                              rInterfaceNegativeSideShapeFunctionsGradientsValues,
                                              rInterfaceNegativeSideWeightsValues,
                                              mpTriangleSplitter->GetNegativeInterfaces(),//mNegativeInterfaces,
                                              mpTriangleSplitter->GetNegativeSubdivisions(),//mNegativeSubdivisions,
                                              mpTriangleSplitter->GetNegativeInterfacesParentIds(),//mNegativeInterfacesParentIds,
                                              p_matrix,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

    // Returns all the shape function values for the contact line.
    void Triangle2D3ModifiedShapeFunctions::ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
        std::vector<unsigned int> &ContactLineIndices,
        std::vector<Matrix> &rContactLineNegativeSideShapeFunctionsValues,
        std::vector<ShapeFunctionsGradientsType> &rContactLineNegativeSideShapeFunctionsGradientsValues,
        std::vector<Vector> &rContactLineNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod)
    {

        rContactLineNegativeSideShapeFunctionsValues.clear();
        rContactLineNegativeSideShapeFunctionsGradientsValues.clear();
        rContactLineNegativeSideWeightsValues.clear();

        if (ContactLineIndices.size() > 0)
        {
            //std::vector<unsigned int> &contact_interface_ids = mpTetrahedraSplitter->mContactInterface;

            // Get the interface condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            for (unsigned int i_cl = 0; i_cl < ContactLineIndices.size(); i_cl++)
            {
                const unsigned int cl_id = ContactLineIndices[i_cl];
                const unsigned int i_parent =
                    //mpTetrahedraSplitter->GetNegativeInterfacesParentIds()[(contact_interface_ids[cl_id])]; //mNegativeInterfacesParentIds[(contact_interface_ids[cl_id])];
                    mpTriangleSplitter->GetNegativeInterfacesParentIds()[(mpTriangleSplitter->GetContactInterface()[cl_id])]; //mNegativeInterfacesParentIds[(contact_interface_ids[cl_id])];
                auto r_parent_geom = mpTriangleSplitter->GetNegativeSubdivisions()[i_parent];      //mNegativeSubdivisions[i_parent];
                //const auto &r_interface_geom = (mpTetrahedraSplitter->mContactLine)[cl_id];
                auto r_interface_geom = (mpTriangleSplitter->GetContactLine())[cl_id];

                Matrix contact_line_negative_side_shape_function_values;
                ShapeFunctionsGradientsType contact_line_negative_side_shape_function_gradient_values;
                Vector contact_line_negative_side_weight_values;
                ////////
                std::cout << "cl_id: " << cl_id << std::endl;
                std::cout << "Contact interface size: " << mpTriangleSplitter->GetContactInterface().size() << std::endl;
                std::cout << "Contact interface value: " << mpTriangleSplitter->GetContactInterface()[cl_id] << std::endl;
                std::cout << "Negative interfaces ParentIds size: " << mpTriangleSplitter->GetNegativeInterfacesParentIds().size() << std::endl;
                std::cout << "i_parent: " << i_parent << std::endl;
                if (mpTriangleSplitter->GetNegativeSubdivisions().empty()) {
                    std::cerr << "Error: No negative subdivisions found." << std::endl;
                }
                std::cout << "Negative subdivisions size: " << mpTriangleSplitter->GetNegativeSubdivisions().size() << std::endl;
                if (!r_parent_geom) {
                    std::cerr << "Error: r_parent_geom is nullptr at index " << i_parent << std::endl;
                }
                std::cout << "r_parent_geom address: " << r_parent_geom << std::endl;
                std::cout << "Info of parent geometry at index " << i_parent << ": " << r_parent_geom->Info() << std::endl;
                std::cout << "Condensation Matrix (P) size: " << p_matrix.size1() << "x" << p_matrix.size2() << std::endl;
                ////////

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
    };

// Given a face id, computes the positive side subdivision shape function values in that face.
void Triangle2D3ModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rPositiveExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
    Vector &rPositiveExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(
            p_matrix);//,
            //mpTriangleSplitter->mEdgeNodeI,
            //mpTriangleSplitter->mEdgeNodeJ,
            //mpTriangleSplitter->mSplitEdges);

        // Get the external faces
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        mpTriangleSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTriangleSplitter->GetPositiveSubdivisions(),//mPositiveSubdivisions,
            FaceId);
        
        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rPositiveExteriorFaceShapeFunctionsValues,
            rPositiveExteriorFaceShapeFunctionsGradientsValues,
            rPositiveExteriorFaceWeightsValues,
            exterior_faces_vector,
            mpTriangleSplitter->GetPositiveSubdivisions(),//mPositiveSubdivisions,
            exterior_faces_parent_ids_vector,
            p_matrix,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Given a face id, computes the positive side subdivision shape function values in that face.
void Triangle2D3ModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
    Matrix &rNegativeExteriorFaceShapeFunctionsValues,
    ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
    Vector &rNegativeExteriorFaceWeightsValues,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod
) {
    if (this->IsSplit()) {
        // Get the condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(
            p_matrix);//,
            //mpTriangleSplitter->mEdgeNodeI,
            //mpTriangleSplitter->mEdgeNodeJ,
            //mpTriangleSplitter->mSplitEdges);

        // Get the external faces
        std::vector < unsigned int > exterior_faces_parent_ids_vector;
        std::vector < IndexedPointGeometryPointerType > exterior_faces_vector;
        mpTriangleSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTriangleSplitter->GetNegativeSubdivisions(),//mNegativeSubdivisions,
            FaceId);
        
        // Compute the positive side external face values
        this->ComputeFaceValuesOnOneSide(
            rNegativeExteriorFaceShapeFunctionsValues,
            rNegativeExteriorFaceShapeFunctionsGradientsValues,
            rNegativeExteriorFaceWeightsValues,
            exterior_faces_vector,
            mpTriangleSplitter->GetNegativeSubdivisions(),//mNegativeSubdivisions,
            exterior_faces_parent_ids_vector,
            p_matrix,
            IntegrationMethod);

    } else {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Triangle2D3ModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
    AreaNormalsContainerType &rPositiveSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rPositiveSideInterfaceAreaNormals,
                                              mpTriangleSplitter->GetPositiveInterfaces(),//mPositiveInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// Compute the positive side interface outwards area normal vector values.
void Triangle2D3ModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
    AreaNormalsContainerType &rNegativeSideInterfaceAreaNormals,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rNegativeSideInterfaceAreaNormals,
                                              mpTriangleSplitter->GetNegativeInterfaces(),//mNegativeInterfaces,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side face outwards area normal vector values.
void Triangle2D3ModifiedShapeFunctions::ComputePositiveExteriorFaceAreaNormals(
    AreaNormalsContainerType &rPositiveExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTriangleSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTriangleSplitter->GetPositiveSubdivisions(),//mPositiveSubdivisions,
            FaceId);
        
        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rPositiveExteriorFaceAreaNormal,
                                              exterior_faces_vector,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputePositiveExteriorFaceAreaNormals method for a non divided geometry.";
    }
};

// For a given face, computes the positive side face outwards area normal vector values.
void Triangle2D3ModifiedShapeFunctions::ComputeNegativeExteriorFaceAreaNormals(
    AreaNormalsContainerType &rNegativeExteriorFaceAreaNormal,
    const unsigned int FaceId,
    const IntegrationMethodType IntegrationMethod) {

    if (this->IsSplit()) {
        // Get the external faces
        std::vector<unsigned int> exterior_faces_parent_ids_vector;
        std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
        mpTriangleSplitter->GenerateExteriorFaces(
            exterior_faces_vector,
            exterior_faces_parent_ids_vector,
            mpTriangleSplitter->GetNegativeSubdivisions(),//mNegativeSubdivisions,
            FaceId);

        // Compute the positive side interface outwars area normal values
        this->ComputeFaceNormalOnOneSide(rNegativeExteriorFaceAreaNormal,
                                              exterior_faces_vector,
                                              IntegrationMethod);
    } else {
        KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceAreaNormals method for a non divided geometry.";
    }
};

// Computes the positive side shape function values in the edges intersections
void Triangle2D3ModifiedShapeFunctions::ComputeShapeFunctionsOnPositiveEdgeIntersections(
    Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues){

    if (this->IsSplit()) {
        // Get the interface condensation matrix
        Matrix p_matrix;
        this->SetCondensationMatrix(
            p_matrix);//,
            //mpTriangleSplitter->mEdgeNodeI,
            //mpTriangleSplitter->mEdgeNodeJ,
            //mpTriangleSplitter->mSplitEdges);

        // Compute the edge intersections shape function values
        this->ComputeEdgeIntersectionValuesOnOneSide(
            p_matrix,
            rPositiveEdgeIntersectionsShapeFunctionsValues);

    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnPositiveEdgeIntersections method for a non divided geometry.";
    }
};

// Computes the negative side shape function values in the edges intersections
void Triangle2D3ModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(
    Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues){
    
    if (this->IsSplit()) {
        // Note that positive and negative sides values are equal for standard shape functions
        this->ComputeShapeFunctionsOnPositiveEdgeIntersections(
            rNegativeEdgeIntersectionsShapeFunctionsValues);
    } else {
        KRATOS_ERROR << "Using the ComputeShapeFunctionsOnNegativeEdgeIntersections method for a non divided geometry.";
    }
};

void Triangle2D3ModifiedShapeFunctions::ComputeNegativeSideContactLineVector(
    std::vector<unsigned int> &FaceIndices,
    std::vector<Vector> &rNegativeSideContactLineVector)
{
    FaceIndices.clear();
    rNegativeSideContactLineVector.clear();

    FaceIndices = mpTriangleSplitter->GetContactFace();

    if (FaceIndices.size() > 0)
    {

        for (unsigned int i_cl = 0; i_cl < FaceIndices.size(); i_cl++)
        {

            Vector aux_vector(3);
            aux_vector[0] = 0.0;
            aux_vector[1] = 0.0;
            aux_vector[2] = -1.0;

            rNegativeSideContactLineVector.push_back(aux_vector);
        }
    }
};

}; //namespace Kratos
