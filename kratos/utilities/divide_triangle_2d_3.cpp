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
#include "utilities/split_triangle.h"
#include "utilities/divide_triangle_2d_3.h"

namespace Kratos
{
    /// Default constructor
    template<class TPointType>
    DivideTriangle2D3<TPointType>::DivideTriangle2D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
        DivideGeometry<TPointType>(rInputGeometry, rNodalDistances), m_structure_node_id(Vector()) {};

    template<class TPointType>
    DivideTriangle2D3<TPointType>::DivideTriangle2D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances, const Vector& rStructureNodes) :
        DivideGeometry<TPointType>(rInputGeometry, rNodalDistances), m_structure_node_id(rStructureNodes) {};

    /// Destructor
    template<class TPointType>
    DivideTriangle2D3<TPointType>::~DivideTriangle2D3() {};

    /// Turn back information as a string.
    template<class TPointType>
    std::string DivideTriangle2D3<TPointType>::Info() const {
        return "Triangle divide operations utility.";
    };

    /// Print information about this object.
    template<class TPointType>
    void DivideTriangle2D3<TPointType>::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Triangle divide operations utility.";
    };

    /// Print object's data.
    template<class TPointType>
    void DivideTriangle2D3<TPointType>::PrintData(std::ostream& rOStream) const {
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Triangle divide operations utility constructed with:\n";
        rOStream << "   Geometry type: " << geometry.Info() << "\n";
        std::stringstream distances_buffer;
        std::stringstream stm;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            stm << nodal_distances(i);
            distances_buffer << stm.str() << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    // Returns the mEdgeNodeI member vector
    template<class TPointType>
    const std::vector<int>& DivideTriangle2D3<TPointType>::GetEdgeIdsI() const {
        return this->mEdgeNodeI;
    }

    // Returns the mEdgeNodeJ member vector
    template<class TPointType>
    const std::vector<int>& DivideTriangle2D3<TPointType>::GetEdgeIdsJ() const {
        return this->mEdgeNodeJ;
    }

    // Returns the mSplitEdges member vector
    template<class TPointType>
    std::vector<int>& DivideTriangle2D3<TPointType>::GetSplitEdges() {
        return this->mSplitEdges;
    }

    // Performs and saves the splitting pattern.
    template<class TPointType>
    void DivideTriangle2D3<TPointType>::GenerateDivision() {

        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();

        // Fill the auxiliar points vector set
        if (this->mIsSplit) {

            // Set the triangle geometry parameters
            const unsigned int n_nodes = 3;
            const unsigned int n_edges = 3;

            // Clear the auxiliar vector points set
            this->mAuxPointsContainer.clear();
            this->mAuxPointsContainer.reserve(6);

            // Store the contact line node (IDs), DEPRECATED: only 1 ID is expected
            std::vector<int> contact_point_node_ids;

            this->mContactFace.clear();

            // Clear the subdivision vectors
            this->mPositiveSubdivisions.clear();
            this->mNegativeSubdivisions.clear();
            this->mPositiveSubdivisions.reserve(2);
            this->mNegativeSubdivisions.reserve(2);

            // Add the original geometry points
            std::vector<int> gl_ids_split_edges(this->mSplitEdges);
            for (unsigned int i = 0; i < n_nodes; ++i) {
                gl_ids_split_edges[i] = geometry[i].Id();
                const array_1d<double, 3> aux_point_coords = geometry[i].Coordinates();
                IndexedPointPointerType paux_point = Kratos::make_shared<IndexedPoint>(aux_point_coords, i);
                this->mAuxPointsContainer.push_back(paux_point);
            }

            // Decide the splitting pattern
            unsigned int aux_node_id = n_nodes;

            for(unsigned int idedge = 0; idedge < n_edges; ++idedge) {
                const unsigned int edge_node_i = this->mEdgeNodeI[idedge];
                const unsigned int edge_node_j = this->mEdgeNodeJ[idedge];

                // Check if the edge is split
                if(nodal_distances(edge_node_i) * nodal_distances(edge_node_j) < 0) {
                    // Set the new node id. in the split edge array corresponding slot
                    this->mSplitEdges[idedge + n_nodes] = aux_node_id;
                    gl_ids_split_edges[idedge + n_nodes] = aux_node_id;

                    // Edge nodes coordinates
                    const array_1d<double, 3> i_node_coords = geometry[edge_node_i].Coordinates();
                    const array_1d<double, 3> j_node_coords = geometry[edge_node_j].Coordinates();

                    // Compute the coordinates of the point on the edge
                    const double aux_node_rel_location = std::abs (nodal_distances(edge_node_i)/(nodal_distances(edge_node_j)-nodal_distances(edge_node_i)));
                    array_1d<double, 3> aux_point_coords;
                    for (unsigned int comp = 0; comp < 3; ++comp) {
                        aux_point_coords(comp) = j_node_coords(comp)*aux_node_rel_location + i_node_coords(comp)*(1.0-aux_node_rel_location);
                    }

                    // Add the intersection point to the auxiliar points array
                    IndexedPointPointerType paux_point = Kratos::make_shared<IndexedPoint>(aux_point_coords, aux_node_id);
                    this->mAuxPointsContainer.push_back(paux_point);

                    if ( m_structure_node_id(edge_node_i) == 1.0 && m_structure_node_id(edge_node_j) == 1.0 ){
                        //KRATOS_INFO("DivideTetrahedra3D4::GenerateDivision()") << idedge + n_nodes << ", " << aux_node_id << std::endl;
                        contact_point_node_ids.push_back(aux_node_id);
                        this->mContactPoint.push_back(paux_point);
                        
                        // Create a PointerVector with the single point
                        Kratos::PointerVector<Kratos::IndexedPoint> point_vector;
                        point_vector.push_back(paux_point);
    
                        // Create the geometry with the points vector
                        auto geometry_point = Kratos::make_shared<Kratos::Geometry<Kratos::IndexedPoint>>(point_vector);

                        this->mContactLine.push_back(geometry_point);

                        this->mContactFace.push_back(idedge);
                        std::cout << "contact_point_node_id: " << aux_node_id << "  contact_line_node_cords: "<< paux_point->Coordinates();
                    }
                }

                aux_node_id++;
            }

            if (!contact_point_node_ids.empty()) {
            std::cout << "contact_line_node_ids: " << contact_point_node_ids ;
                for (const auto& contact_point : this->mContactPoint) {
                   std::cout << " contact_line_node_cords: " << contact_point->Coordinates();
                    }}
            if (!contact_point_node_ids.empty()) {
                std::cout << "contact_line_node_ids: " << contact_point_node_ids;
                for (const auto& geometry_point : this->mContactLine) {
                    const auto& point = geometry_point->Points()[0];
                    std::cout << " contact_line_node_cords2: " << point.Coordinates();
                }
            }

            // Call the splitting mode computation function
            std::vector<int> edge_ids(3);
            TriangleSplit::TriangleSplitMode(gl_ids_split_edges.data(), edge_ids.data());

            // Call the splitting function
            std::vector<int> t(12);     // Ids of the generated subdivisions
            int n_int = 0;              // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            TriangleSplit::Split_Triangle(edge_ids.data(), t.data(), &this->mDivisionsNumber, &this->mSplitEdgesNumber, &n_int);

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < this->mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2;
                TriangleSplit::TriangleGetNewConnectivityGID(idivision, t.data(), this->mSplitEdges.data(), &i0, &i1, &i2);

                // Generate a pointer to an auxiliar triangular geometry made with the subdivision points
                IndexedPointGeometryPointerType p_aux_partition = GenerateAuxiliaryPartitionTriangle(i0, i1, i2);

                // Determine if the subdivision is whether in the negative or the positive side
                // Note that zero distance nodes are also identified and stored in here
                unsigned int neg = 0, pos = 0;
                if(i0 <= 2) {if(nodal_distances(i0) < 0.0) neg++; else if(nodal_distances(i0) > 0.0) pos++; else this->mNodeIsCut.set(i0);};
                if(i1 <= 2) {if(nodal_distances(i1) < 0.0) neg++; else if(nodal_distances(i1) > 0.0) pos++; else this->mNodeIsCut.set(i1);};
                if(i2 <= 2) {if(nodal_distances(i2) < 0.0) neg++; else if(nodal_distances(i2) > 0.0) pos++; else this->mNodeIsCut.set(i2);};

                KRATOS_ERROR_IF(neg > 0 && pos > 0) << "The subgeometry " << i0 << " " << i1 << " " << i2 << " in triangle has nodes in both positive and negative sides." << std::endl;

                bool is_positive = false;
                if(pos > 0) {is_positive = true;}

                // Add the generated triangle to its corresponding partition arrays
                if (is_positive) {
                    this->mPositiveSubdivisions.push_back(p_aux_partition);
                } else {
                    this->mNegativeSubdivisions.push_back(p_aux_partition);
                }
            }

        } else {
            this->mDivisionsNumber = 1;
            this->mSplitEdgesNumber = 0;
        }
    };

    template<class TPointType>
    void DivideTriangle2D3<TPointType>::GenerateIntersectionsSkin() {

        // Set some geometry constant parameters
        const unsigned int n_faces = 3;

        // Clear the interfaces vectors
        this->mPositiveInterfaces.clear();
        this->mNegativeInterfaces.clear();
        this->mPositiveInterfaces.reserve(1);
        this->mNegativeInterfaces.reserve(1);
        this->mPositiveInterfacesParentIds.clear();
        this->mNegativeInterfacesParentIds.clear();
        this->mPositiveInterfacesParentIds.reserve(1);
        this->mNegativeInterfacesParentIds.reserve(1);

        if (this->mIsSplit) {

            const unsigned int n_positive_subdivision = this->mPositiveSubdivisions.size();
            const unsigned int n_negative_subdivision = this->mNegativeSubdivisions.size();

            // Compute the positive side intersection geometries
            for (unsigned int i_subdivision = 0; i_subdivision < n_positive_subdivision; ++i_subdivision) {
                // Get the subdivision geometry
                const IndexedPointGeometryType& r_subdivision_geom = *this->mPositiveSubdivisions[i_subdivision];

                // Faces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    // Get the subdivision face nodal keys
                    int node_i_key = r_subdivision_geom[mEdgeNodeI[i_face]].Id();
                    int node_j_key = r_subdivision_geom[mEdgeNodeJ[i_face]].Id();

                    // Check the nodal keys to state which nodes belong to the interface
                    // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliary interface points
                    // For the zero distance case, the corresponding node is considered as part of the interface
                    if (NodeIsInterface(node_i_key) && NodeIsInterface(node_j_key)) {
                        // Generate an indexed point line geometry pointer with the two interface nodes
                        IndexedPointGeometryPointerType p_intersection_line = this->GenerateIntersectionLine(node_i_key, node_j_key);
                        this->mPositiveInterfaces.push_back(p_intersection_line);
                        this->mPositiveInterfacesParentIds.push_back(i_subdivision);

                        // In triangles, a unique face can belong to the interface
                        break;
                    }
                }
            }

            // Compute the negative side intersection geometries
            for (unsigned int i_subdivision = 0; i_subdivision < n_negative_subdivision; ++i_subdivision) {
                // Get the subdivision geometry
                const IndexedPointGeometryType& r_subdivision_geom = *this->mNegativeSubdivisions[i_subdivision];

                // Faces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    // Get the subdivision face nodal keys
                    int node_i_key = r_subdivision_geom[mEdgeNodeI[i_face]].Id();
                    int node_j_key = r_subdivision_geom[mEdgeNodeJ[i_face]].Id();

                    // Check the nodal keys to state which nodes belong to the interface
                    // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliary interface points
                    // For the zero distance case, the corresponding node is considered as part of the interface
                    if (NodeIsInterface(node_i_key) && NodeIsInterface(node_j_key)) {
                        // Generate an indexed point line geometry pointer with the two interface nodes
                        IndexedPointGeometryPointerType p_intersection_line = this->GenerateIntersectionLine(node_i_key ,node_j_key);
                        this->mNegativeInterfaces.push_back(p_intersection_line);
                        this->mNegativeInterfacesParentIds.push_back(i_subdivision);

                        //////////
                        //KRATOS_INFO("DivideTriangle2D3::GenerateDivision()") <<  "mContactFace.size=" << this->mContactFace.size() << std::endl;
                        //KRATOS_INFO("DivideTriangle2D3::GenerateDivision()") <<  "mContactPoint.size" << this->mContactPoint.size() << std::endl;
                        //////////

                        // if (mContactFace.size() > 0){
                        if (this->mContactFace.size() > 0){

                            for (unsigned int i_contact = 0; i_contact < this->mContactFace.size(); i_contact++){

                                this->mContactInterface.push_back( this->mNegativeInterfaces.size() - 1 );
                                //KRATOS_INFO("this->mNegativeInterfaces.size()") <<  this->mNegativeInterfaces.size() << std::endl;
                            }
                        }
                        //KRATOS_INFO("DivideTetrahedra3D4::GenerateDivision()") <<  "mContactInterface.size()=" << this->mContactInterface.size() << "mContactEdge.size()="<< this->mContactEdge.size()<<std::endl; 

                        // In triangles, a unique face can belong to the interface
                        break;
                    }
                }
            }
        } else {
            KRATOS_ERROR << "Trying to generate the intersection skin in DivideTriangle2D3::GenerateIntersectionsSkin() for a non-split element.";
        }
    };

    template<class TPointType>
    void DivideTriangle2D3<TPointType>::GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer) {

        // Set some geometry constant parameters
        const unsigned int n_faces = 3;

        // Set the exterior faces vector
        rExteriorFacesVector.clear();
        rExteriorFacesVector.reserve(3);
        rExteriorFacesParentSubdivisionsIdsVector.clear();
        rExteriorFacesParentSubdivisionsIdsVector.reserve(3);

        // Iterate the triangle faces
        for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
            std::vector < unsigned int > aux_ext_faces_parent_ids;
            std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > aux_ext_faces;

            DivideTriangle2D3::GenerateExteriorFaces(
                aux_ext_faces,
                aux_ext_faces_parent_ids,
                rSubdivisionsContainer,
                i_face);

            rExteriorFacesVector.insert(rExteriorFacesVector.end(), aux_ext_faces.begin(), aux_ext_faces.end());
            rExteriorFacesParentSubdivisionsIdsVector.insert(rExteriorFacesParentSubdivisionsIdsVector.end(), aux_ext_faces_parent_ids.begin(), aux_ext_faces_parent_ids.end());
        }
    };

    template<class TPointType>
    void DivideTriangle2D3<TPointType>::GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer,
        const unsigned int FatherFaceId) {
        // Set some geometry constant parameters
        const unsigned int n_faces = 3;

        // Set the exterior faces vector
        rExteriorFacesVector.clear();
        rExteriorFacesVector.reserve(2);
        rExteriorFacesParentSubdivisionsIdsVector.clear();
        rExteriorFacesParentSubdivisionsIdsVector.reserve(2);

        if (this->mIsSplit) {
            // Create the face nodes data
            // The position represents the face while the value real and intersection nodes in that face edges
            std::array < std::array < unsigned int, 3 >, 3 > edges_map = {{
                {{1, 2, 4}},     // Face 0
                {{2, 0, 5}},     // Face 1
                {{0, 1, 3}}}};   // Face 2

            // Compute the side exterior faces geometries
            const unsigned int n_subdivision = rSubdivisionsContainer.size();
            for (unsigned int i_subdivision = 0; i_subdivision < n_subdivision; ++i_subdivision) {
                // Get the subdivision geometry
                const IndexedPointGeometryType& r_subdivision_geom = *rSubdivisionsContainer[i_subdivision];

                // Subdivision geometry subfaces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    // Get the subdivision subface nodal keys
                    int node_i_key = r_subdivision_geom[mEdgeNodeI[i_face]].Id();
                    int node_j_key = r_subdivision_geom[mEdgeNodeJ[i_face]].Id();

                    // Get the candidate nodes
                    std::array< unsigned int, 3 > faces_edge_nodes = edges_map[FatherFaceId];

                    // Search the subdivision nodal keys into the parent geometry face key value
                    if (std::find(faces_edge_nodes.begin(), faces_edge_nodes.end(), node_i_key) != faces_edge_nodes.end()) {
                        if (std::find(faces_edge_nodes.begin(), faces_edge_nodes.end(), node_j_key) != faces_edge_nodes.end()) {
                            // If both nodes are in the candidate nodes list, the subface is exterior
                            IndexedPointGeometryPointerType p_subface_line = this->GenerateIntersectionLine(node_i_key, node_j_key);

                            rExteriorFacesVector.push_back(p_subface_line);
                            rExteriorFacesParentSubdivisionsIdsVector.push_back(i_subdivision);
                        }
                    }
                }
            }
        } else {
            KRATOS_ERROR << "Trying to generate the exterior faces in DivideTriangle2D3::GenerateExteriorFaces() for a non-split element.";
        }
    };

    template<class TPointType>
    typename DivideTriangle2D3<TPointType>::IndexedPointGeometryPointerType DivideTriangle2D3<TPointType>::GenerateAuxiliaryPartitionTriangle(
        const int I0,
        const int I1,
        const int I2)
    {
        return Kratos::make_shared<IndexedPointTriangleType>(
            this->mAuxPointsContainer(I0),
            this->mAuxPointsContainer(I1),
            this->mAuxPointsContainer(I2));
    };

    template<class TPointType>
    typename DivideTriangle2D3<TPointType>::IndexedPointGeometryPointerType DivideTriangle2D3<TPointType>::GenerateIntersectionLine(
        const int I0,
        const int I1)
    {
        return Kratos::make_shared<IndexedPointLineType>(
            this->mAuxPointsContainer(I0),
            this->mAuxPointsContainer(I1));
    };

    template<class TPointType>
    bool DivideTriangle2D3<TPointType>::NodeIsInterface(int NodeKey) const
    {
        constexpr int num_nodes = 3;
        return NodeKey >= num_nodes || mNodeIsCut[NodeKey];
    }


    template class DivideTriangle2D3<Node>;
    template class DivideTriangle2D3<IndexedPoint>;

};
