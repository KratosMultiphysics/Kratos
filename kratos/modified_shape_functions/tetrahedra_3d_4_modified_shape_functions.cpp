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
    Tetrahedra3D4ModifiedShapeFunctions::Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector &rNodalDistances) : ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
                                                                                                                                                        mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4<Node>>(*pInputGeometry, rNodalDistances))
    {

        // Perform the element splitting
        mpTetrahedraSplitter->GenerateDivision();
        mpTetrahedraSplitter->GenerateIntersectionsSkin();
    };

    //
    Tetrahedra3D4ModifiedShapeFunctions::Tetrahedra3D4ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector &rNodalDistances, const Vector& rStructureNodes) : ModifiedShapeFunctions(pInputGeometry, rNodalDistances),
                                                                                                                                                                                                mpTetrahedraSplitter(Kratos::make_shared<DivideTetrahedra3D4<Node>>(*pInputGeometry, rNodalDistances, rStructureNodes))
    {

        // Perform the element splitting
        mpTetrahedraSplitter->GenerateDivision();
        mpTetrahedraSplitter->GenerateIntersectionsSkin();
    };

    /// Destructor
    Tetrahedra3D4ModifiedShapeFunctions::~Tetrahedra3D4ModifiedShapeFunctions() {};

    /// Turn back information as a string.
    std::string Tetrahedra3D4ModifiedShapeFunctions::Info() const
    {
        return "Tetrahedra3D4N modified shape functions computation class.";
    };

    /// Print information about this object.
    void Tetrahedra3D4ModifiedShapeFunctions::PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "Tetrahedra3D4N modified shape functions computation class.";
    };

    /// Print object's data.
    void Tetrahedra3D4ModifiedShapeFunctions::PrintData(std::ostream &rOStream) const
    {
        const GeometryPointerType p_geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Tetrahedra3D4N modified shape functions computation class:\n";
        rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
        std::stringstream distances_buffer;
        std::stringstream stm;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i)
        {
            stm << nodal_distances(i);
            distances_buffer << stm.str() << " ";
        }
        rOStream << "\tDistance values: " << distances_buffer.str();
    };

    // Returns a pointer to the splitting utility
    const DivideGeometry<Node>::Pointer Tetrahedra3D4ModifiedShapeFunctions::pGetSplittingUtil() const
    {
        return mpTetrahedraSplitter;
    };

    // Returns true if the element is splitting
    bool Tetrahedra3D4ModifiedShapeFunctions::IsSplit()
    {
        return mpTetrahedraSplitter->mIsSplit;
    };

    void Tetrahedra3D4ModifiedShapeFunctions::SetCondensationMatrix(Matrix &rIntPointCondMatrix)
    {
        ModifiedShapeFunctions::SetCondensationMatrix(
            rIntPointCondMatrix,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);
    }

    void Tetrahedra3D4ModifiedShapeFunctions::SetPositiveSideCondensationMatrix(Matrix &rPosSideCondMatrix)
    {
        ModifiedShapeFunctions::SetCondensationMatrix(
            rPosSideCondMatrix,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);
    }

    void Tetrahedra3D4ModifiedShapeFunctions::SetNegativeSideCondensationMatrix(Matrix &rNegSideCondMatrix)
    {
        ModifiedShapeFunctions::SetCondensationMatrix(
            rNegSideCondMatrix,
            mpTetrahedraSplitter->mEdgeNodeI,
            mpTetrahedraSplitter->mEdgeNodeJ,
            mpTetrahedraSplitter->mSplitEdges);
    }

    // Internally computes the splitting pattern and returns all the shape function values for the positive side.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rPositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveSideShapeFunctionsGradientsValues,
        Vector &rPositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod)
    {

        if (this->IsSplit())
        {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Compute the positive side values
            this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues,
                                         rPositiveSideShapeFunctionsGradientsValues,
                                         rPositiveSideWeightsValues,
                                         mpTetrahedraSplitter->GetPositiveSubdivisions(), // mPositiveSubdivisions,
                                         p_matrix,
                                         IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the negative side.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeSideShapeFunctionsGradientsValues,
        Vector &rNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod)
    {

        if (this->IsSplit())
        {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Compute the negative side values
            this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues,
                                         rNegativeSideShapeFunctionsGradientsValues,
                                         rNegativeSideWeightsValues,
                                         mpTetrahedraSplitter->GetNegativeSubdivisions(), //->mNegativeSubdivisions,
                                         p_matrix,
                                         IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the positive interface side.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfacePositiveSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfacePositiveSideShapeFunctionsGradientsValues,
        Vector &rInterfacePositiveSideWeightsValues,
        const IntegrationMethodType IntegrationMethod)
    {

        if (this->IsSplit())
        {
            // Get the interface condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Compute the positive side interface values
            this->ComputeFaceValuesOnOneSide(rInterfacePositiveSideShapeFunctionsValues,
                                             rInterfacePositiveSideShapeFunctionsGradientsValues,
                                             rInterfacePositiveSideWeightsValues,
                                             mpTetrahedraSplitter->GetPositiveInterfaces(),          // mPositiveInterfaces,
                                             mpTetrahedraSplitter->GetPositiveSubdivisions(),        // mPositiveSubdivisions,
                                             mpTetrahedraSplitter->GetPositiveInterfacesParentIds(), // mPositiveInterfacesParentIds,
                                             p_matrix,
                                             IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the negative interface side.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        Matrix &rInterfaceNegativeSideShapeFunctionsValues,
        ShapeFunctionsGradientsType &rInterfaceNegativeSideShapeFunctionsGradientsValues,
        Vector &rInterfaceNegativeSideWeightsValues,
        const IntegrationMethodType IntegrationMethod)
    {

        if (this->IsSplit())
        {
            // Get the interface condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Compute the positive side interface values
            this->ComputeFaceValuesOnOneSide(rInterfaceNegativeSideShapeFunctionsValues,
                                             rInterfaceNegativeSideShapeFunctionsGradientsValues,
                                             rInterfaceNegativeSideWeightsValues,
                                             mpTetrahedraSplitter->GetNegativeInterfaces(),          // mNegativeInterfaces,
                                             mpTetrahedraSplitter->GetNegativeSubdivisions(),        // mNegativeSubdivisions,
                                             mpTetrahedraSplitter->GetNegativeInterfacesParentIds(), // mNegativeInterfacesParentIds,
                                             p_matrix,
                                             IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Returns all the shape function values for the contact line.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
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
                    mpTetrahedraSplitter->GetNegativeInterfacesParentIds()[(mpTetrahedraSplitter->GetContactInterface()[cl_id])]; //mNegativeInterfacesParentIds[(contact_interface_ids[cl_id])];
                auto r_parent_geom = mpTetrahedraSplitter->GetNegativeSubdivisions()[i_parent];      //mNegativeSubdivisions[i_parent];
                //const auto &r_interface_geom = (mpTetrahedraSplitter->mContactLine)[cl_id];
                auto r_interface_geom = (mpTetrahedraSplitter->GetContactLine())[cl_id];

                Matrix contact_line_negative_side_shape_function_values;
                ShapeFunctionsGradientsType contact_line_negative_side_shape_function_gradient_values;
                Vector contact_line_negative_side_weight_values;
                ////////
                std::cout << "cl_id: " << cl_id << std::endl;
                std::cout << "Contact interface size: " << mpTetrahedraSplitter->GetContactInterface().size() << std::endl;
                std::cout << "Contact interface value: " << mpTetrahedraSplitter->GetContactInterface()[cl_id] << std::endl;
                std::cout << "Negative interfaces ParentIds size: " << mpTetrahedraSplitter->GetNegativeInterfacesParentIds().size() << std::endl;
                std::cout << "i_parent: " << i_parent << std::endl;
                if (mpTetrahedraSplitter->GetNegativeSubdivisions().empty()) {
                    std::cerr << "Error: No negative subdivisions found." << std::endl;
                }
                std::cout << "Negative subdivisions size: " << mpTetrahedraSplitter->GetNegativeSubdivisions().size() << std::endl;
                if (!r_parent_geom) {
                    std::cerr << "Error: r_parent_geom is nullptr at index " << i_parent << std::endl;
                }
                std::cout << "r_parent_geom address: " << r_parent_geom << std::endl;
                std::cout << "Info of parent geometry at index " << i_parent << ": " << r_parent_geom->Info() << std::endl;
                printf("Matrix size: %d x %d\n", p_matrix.size1(), p_matrix.size2()); // Correct (prints in decimal)
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
    void Tetrahedra3D4ModifiedShapeFunctions::ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rPositiveExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rPositiveExteriorFaceShapeFunctionsGradientsValues,
        Vector &rPositiveExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod)
    {
        if (this->IsSplit())
        {
            // Get the condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Get the external faces
            std::vector<unsigned int> exterior_faces_parent_ids_vector;
            std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
            mpTetrahedraSplitter->GenerateExteriorFaces(
                exterior_faces_vector,
                exterior_faces_parent_ids_vector,
                mpTetrahedraSplitter->GetPositiveSubdivisions(), // mPositiveSubdivisions,
                FaceId);

            // Compute the positive side external face values
            this->ComputeFaceValuesOnOneSide(
                rPositiveExteriorFaceShapeFunctionsValues,
                rPositiveExteriorFaceShapeFunctionsGradientsValues,
                rPositiveExteriorFaceWeightsValues,
                exterior_faces_vector,
                mpTetrahedraSplitter->GetPositiveSubdivisions(), // mPositiveSubdivisions,
                exterior_faces_parent_ids_vector,
                p_matrix,
                IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputePositiveExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Given a face id, computes the positive side subdivision shape function values in that face.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues(
        Matrix &rNegativeExteriorFaceShapeFunctionsValues,
        ShapeFunctionsGradientsType &rNegativeExteriorFaceShapeFunctionsGradientsValues,
        Vector &rNegativeExteriorFaceWeightsValues,
        const unsigned int FaceId,
        const IntegrationMethodType IntegrationMethod)
    {
        if (this->IsSplit())
        {
            // Get the condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Get the external faces
            std::vector<unsigned int> exterior_faces_parent_ids_vector;
            std::vector<IndexedPointGeometryPointerType> exterior_faces_vector;
            mpTetrahedraSplitter->GenerateExteriorFaces(
                exterior_faces_vector,
                exterior_faces_parent_ids_vector,
                mpTetrahedraSplitter->GetNegativeSubdivisions(), // mNegativeSubdivisions,
                FaceId);

            // Compute the positive side external face values
            this->ComputeFaceValuesOnOneSide(
                rNegativeExteriorFaceShapeFunctionsValues,
                rNegativeExteriorFaceShapeFunctionsGradientsValues,
                rNegativeExteriorFaceWeightsValues,
                exterior_faces_vector,
                mpTetrahedraSplitter->GetNegativeSubdivisions(), // mNegativeSubdivisions,
                exterior_faces_parent_ids_vector,
                p_matrix,
                IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeNegativeExteriorFaceShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Compute the positive side interface outwards unit normal vector values.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputePositiveSideInterfaceAreaNormals(
        AreaNormalsContainerType &rPositiveSideInterfaceAreaNormals,
        const IntegrationMethodType IntegrationMethod)
    {

        if (this->IsSplit())
        {
            // Compute the positive side interface outwars unit normal values
            this->ComputeFaceNormalOnOneSide(rPositiveSideInterfaceAreaNormals,
                                             mpTetrahedraSplitter->GetPositiveInterfaces(), // mPositiveInterfaces,
                                             IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputePositiveSideInterfaceAreaNormals method for a non divided geometry.";
        }
    };

    // Compute the positive side interface outwards unit normal vector values.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideInterfaceAreaNormals(
        AreaNormalsContainerType &rNegativeSideInterfaceAreaNormals,
        const IntegrationMethodType IntegrationMethod)
    {

        if (this->IsSplit())
        {
            // Compute the positive side interface outwars unit normal values
            this->ComputeFaceNormalOnOneSide(rNegativeSideInterfaceAreaNormals,
                                             mpTetrahedraSplitter->GetNegativeInterfaces(), // mNegativeInterfaces,
                                             IntegrationMethod);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeNegativeSideInterfaceAreaNormals method for a non divided geometry.";
        }
    };

    // For a given face, computes the positive side exteriorface outwards area normal vector values.
    void Tetrahedra3D4ModifiedShapeFunctions::ComputePositiveExteriorFaceAreaNormals(
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
                mpTetrahedraSplitter->GetPositiveSubdivisions(), // mPositiveSubdivisions,
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
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeExteriorFaceAreaNormals(
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
                mpTetrahedraSplitter->GetNegativeSubdivisions(), // mNegativeSubdivisions,
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
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeShapeFunctionsOnPositiveEdgeIntersections(
        Matrix &rPositiveEdgeIntersectionsShapeFunctionsValues)
    {

        if (this->IsSplit())
        {
            // Get the interface condensation matrix
            Matrix p_matrix;
            this->SetCondensationMatrix(p_matrix); //,
                                                   // mpTetrahedraSplitter->mEdgeNodeI,
                                                   // mpTetrahedraSplitter->mEdgeNodeJ,
                                                   // mpTetrahedraSplitter->mSplitEdges);

            // Compute the edge intersections shape function values
            this->ComputeEdgeIntersectionValuesOnOneSide(
                p_matrix,
                rPositiveEdgeIntersectionsShapeFunctionsValues);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeShapeFunctionsOnPositiveEdgeIntersections method for a non divided geometry.";
        }
    };

    // Computes the negative side shape function values in the edges intersections
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeShapeFunctionsOnNegativeEdgeIntersections(
        Matrix &rNegativeEdgeIntersectionsShapeFunctionsValues)
    {

        if (this->IsSplit())
        {
            // Note that positive and negative sides values are equal for standard shape functions
            this->ComputeShapeFunctionsOnPositiveEdgeIntersections(
                rNegativeEdgeIntersectionsShapeFunctionsValues);
        }
        else
        {
            KRATOS_ERROR << "Using the ComputeShapeFunctionsOnNegativeEdgeIntersections method for a non divided geometry.";
        }
    };

    /* bool Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideContactLineVector(
        Vector &rNegativeSideContactLineVector) */
    void Tetrahedra3D4ModifiedShapeFunctions::ComputeNegativeSideContactLineVector(
        std::vector<unsigned int> &FaceIndices,
        std::vector<Vector> &rNegativeSideContactLineVector)
    {
        FaceIndices.clear();
        rNegativeSideContactLineVector.clear();

        //FaceIndices = mpTetrahedraSplitter->mContactFace;
        FaceIndices = mpTetrahedraSplitter->GetContactFace();

        // KRATOS_INFO("ComputeNegativeSideContactLineVector") << (mpTetrahedraSplitter->mContactFace).size() << std::endl;

        if (FaceIndices.size() > 0)
        {

            //std::vector<unsigned int> &contact_interface_ids = mpTetrahedraSplitter->mContactInterface;
            //std::vector<unsigned int> &contact_interface_edge_ids = mpTetrahedraSplitter->mContactEdge;

            for (unsigned int i_cl = 0; i_cl < FaceIndices.size(); i_cl++)
            {
                //const unsigned int i_interface = contact_interface_ids[i_cl];
                const unsigned int i_interface = mpTetrahedraSplitter->GetContactInterface()[i_cl];
                //const unsigned int i_edge = contact_interface_edge_ids[i_cl];
                const unsigned int i_edge = mpTetrahedraSplitter->GetContactEdge()[i_cl];

                // KRATOS_INFO("(contact_interface_ids[i_cl])") << i_interface << std::endl;
                // KRATOS_INFO("(contact_interface_edge_ids[i_cl])") << i_edge << std::endl;

                const auto &edges = (mpTetrahedraSplitter->GetNegativeInterfaces())[i_interface]->Edges(); // mNegativeInterfaces)[i_interface]->Edges();
                const auto &edgei = edges[i_edge];

                // KRATOS_INFO("edgei[0].Coordinates()") << edgei[0].Coordinates() << std::endl;
                // KRATOS_INFO("edgei[1].Coordinates()") << edgei[1].Coordinates() << std::endl;

                Vector aux_vector = edgei[1].Coordinates() - edgei[0].Coordinates();

                // KRATOS_INFO("aux_vector") << aux_vector << std::endl;

                const double norm = Kratos::norm_2(aux_vector);
                aux_vector = 1.0 / norm * aux_vector;

                rNegativeSideContactLineVector.push_back(aux_vector);
            }
        }
    };

}; // namespace Kratos
