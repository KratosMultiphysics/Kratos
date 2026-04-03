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

#include "utilities/split_tetrahedra.h"
#include "utilities/divide_tetrahedra_3d_4.h"

namespace Kratos
{
    /// Default constructor
    template<class TPointType>
    DivideTetrahedra3D4<TPointType>::DivideTetrahedra3D4(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
        DivideGeometry<TPointType>(rInputGeometry, rNodalDistances) {};

    /// Destructor
    template<class TPointType>
    DivideTetrahedra3D4<TPointType>::~DivideTetrahedra3D4() {};

    /// Turn back information as a string.
    template<class TPointType>
    std::string DivideTetrahedra3D4<TPointType>::Info() const {
        return "Tetrahedra divide operations utility.";
    };

    /// Print information about this object.
    template<class TPointType>
    void DivideTetrahedra3D4<TPointType>::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Tetrahedra divide operations utility.";
    };

    /// Print object's data.
    template<class TPointType>
    void DivideTetrahedra3D4<TPointType>::PrintData(std::ostream& rOStream) const {
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Tetrahedra divide operations utility constructed with:\n";
        rOStream << "   Geometry type: " << geometry.Info() << "\n";
        std::stringstream distances_buffer;
        std::ostringstream stm;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            stm << nodal_distances(i);
            distances_buffer << stm.str() << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    // Returns the mEdgeNodeI member vector
    template<class TPointType>
    const std::vector<int>& DivideTetrahedra3D4<TPointType>::GetEdgeIdsI() const {
        return this->mEdgeNodeI;
    }

    // Returns the mEdgeNodeJ member vector
    template<class TPointType>
    const std::vector<int>& DivideTetrahedra3D4<TPointType>::GetEdgeIdsJ() const {
        return this->mEdgeNodeJ;
    }

    // Returns the mSplitEdges member vector
    template<class TPointType>
    std::vector<int>& DivideTetrahedra3D4<TPointType>::GetSplitEdges() {
        return this->mSplitEdges;
    }

    // Performs and saves the splitting pattern.
    template<class TPointType>
    void DivideTetrahedra3D4<TPointType>::GenerateDivision() {

        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();

        // Fill the auxiliar points vector set
        if (this->mIsSplit) {

            // Set the Tetrahedra geometry parameters
            const unsigned int n_nodes = 4;
            const unsigned int n_edges = 6;

            // Clear the auxiliar vector points set
            this->mAuxPointsContainer.clear();
            this->mAuxPointsContainer.reserve(10);

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
                }

                aux_node_id++;
            }

            // Call the splitting mode computation function
            std::vector<int> edge_ids(6);
            TetrahedraSplit::TetrahedraSplitMode(gl_ids_split_edges.data(), edge_ids.data());

            // Call the splitting function
            std::vector<int> t(56);     // Ids of the generated subdivisions
            int n_int = 0;              // Number of internal nodes (set to 0 since it is not needed for tetrahedra splitting)
            TetrahedraSplit::Split_Tetrahedra(edge_ids.data(), t.data(), &this->mDivisionsNumber, &this->mSplitEdgesNumber, &n_int);

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < this->mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2, i3;
                TetrahedraSplit::TetrahedraGetNewConnectivityGID(idivision, t.data(), this->mSplitEdges.data(), &i0, &i1, &i2, &i3);

                // Generate a pointer to an auxiliar triangular geometry made with the subdivision points
                IndexedPointGeometryPointerType p_aux_partition = Kratos::make_shared<IndexedPointTetrahedraType>(
                    this->mAuxPointsContainer(i0),
                    this->mAuxPointsContainer(i1),
                    this->mAuxPointsContainer(i2),
                    this->mAuxPointsContainer(i3));

                // Determine if the subdivision is whether in the negative or the positive side
                // Note that zero distance nodes are also identified and stored in here 
                unsigned int neg = 0, pos = 0;
                if(i0 <= 3) { if (nodal_distances(i0) < 0.0) neg++; else if (nodal_distances(i0) > 0.0) pos++; else mNodeIsCut.set(i0); }
                if(i1 <= 3) { if (nodal_distances(i1) < 0.0) neg++; else if (nodal_distances(i1) > 0.0) pos++; else mNodeIsCut.set(i1); }
                if(i2 <= 3) { if (nodal_distances(i2) < 0.0) neg++; else if (nodal_distances(i2) > 0.0) pos++; else mNodeIsCut.set(i2); }
                if(i3 <= 3) { if (nodal_distances(i3) < 0.0) neg++; else if (nodal_distances(i3) > 0.0) pos++; else mNodeIsCut.set(i3); }

                if(neg > 0 && pos > 0)
                    KRATOS_ERROR << "The subgeometry " << i0 << " " << i1 << " " << i2 << " " << i3 << " in tetrahedra has nodes in both positive and negative sides." << std::endl;

                bool is_positive = false;
                if(pos > 0) {is_positive = true;}

                // Add the generated tetrahedra to its corresponding partition arrays
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
    void DivideTetrahedra3D4<TPointType>::GenerateIntersectionsSkin() {
        // Set some geometry constant parameters
        const unsigned int n_faces = 4;

        // Clear the interfaces vectors
        this->mPositiveInterfaces.clear();
        this->mNegativeInterfaces.clear();
        this->mPositiveInterfacesParentIds.clear();
        this->mNegativeInterfacesParentIds.clear();

        if (this->mIsSplit) {

            const unsigned int n_positive_subdivision = this->mPositiveSubdivisions.size();
            const unsigned int n_negative_subdivision = this->mNegativeSubdivisions.size();

            // Compute the positive side intersection geometries
            for (unsigned int i_subdivision = 0; i_subdivision < n_positive_subdivision; ++i_subdivision) {
                // Get the subdivision geometry faces
                const IndexedPointGeometryPointerType p_subdivision_geom = this->mPositiveSubdivisions[i_subdivision];
                const IndexedGeometriesArrayType subdivision_faces = p_subdivision_geom->GenerateFaces();

                // Faces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    const IndexedPointGeometryType& r_face = subdivision_faces[i_face];

                    // Get the subdivision face nodal keys
                    int node_i_key = r_face[0].Id();
                    int node_j_key = r_face[1].Id();
                    int node_k_key = r_face[2].Id();

                    // Check the nodal keys to state which nodes belong to the interface
                    // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliary interface points
                    // For the zero distance case, the corresponding node is considered as part of the interface
                    if (NodeIsInterface(node_i_key) && NodeIsInterface(node_j_key) && NodeIsInterface(node_k_key)) {
                        // Generate an indexed point triangle geometry pointer with the two interface nodes
                        IndexedPointGeometryPointerType p_intersection_tri = Kratos::make_shared<IndexedPointTriangleType>(this->mAuxPointsContainer(node_i_key),
                                                                                                                          this->mAuxPointsContainer(node_j_key),
                                                                                                                          this->mAuxPointsContainer(node_k_key));
                        this->mPositiveInterfaces.push_back(p_intersection_tri);
                        this->mPositiveInterfacesParentIds.push_back(i_subdivision);
                    }
                }
            }

            // Compute the negative side intersection geometries
            for (unsigned int i_subdivision = 0; i_subdivision < n_negative_subdivision; ++i_subdivision) {
                // Get the subdivision geometry
                const IndexedPointGeometryPointerType p_subdivision_geom = this->mNegativeSubdivisions[i_subdivision];
                const IndexedGeometriesArrayType subdivision_faces = p_subdivision_geom->GenerateFaces();

                // Faces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    const IndexedPointGeometryType& r_face = subdivision_faces[i_face];

                    // Get the subdivision face nodal keys
                    int node_i_key = r_face[0].Id();
                    int node_j_key = r_face[1].Id();
                    int node_k_key = r_face[2].Id();

                    // Check the nodal keys to state which nodes belong to the interface
                    // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliary interface points
                    // For the zero distance case, the corresponding node is considered as part of the interface
                    if (NodeIsInterface(node_i_key) && NodeIsInterface(node_j_key) && NodeIsInterface(node_k_key)) {
                        // Generate an indexed point triangle geometry pointer with the two interface nodes
                        IndexedPointGeometryPointerType p_intersection_tri = Kratos::make_shared<IndexedPointTriangleType>(this->mAuxPointsContainer(node_i_key),
                                                                                                                          this->mAuxPointsContainer(node_j_key),
                                                                                                                          this->mAuxPointsContainer(node_k_key));
                        this->mNegativeInterfaces.push_back(p_intersection_tri);
                        this->mNegativeInterfacesParentIds.push_back(i_subdivision);
                    }
                }
            }
        } else {
            KRATOS_ERROR << "Trying to generate the intersection skin in DivideTetrahedra3D4::GenerateIntersectionsSkin() for a non-split element.";
        }
    };

    template<class TPointType>
    void DivideTetrahedra3D4<TPointType>::GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer) {

        // Set some geometry constant parameters
        const unsigned int n_faces = 4;

        // Set the exterior faces vectors
        rExteriorFacesVector.clear();
        rExteriorFacesParentSubdivisionsIdsVector.clear();

        // Iterate the triangle faces
        for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
            std::vector < unsigned int > aux_ext_faces_parent_ids;
            std::vector < DivideTetrahedra3D4::IndexedPointGeometryPointerType > aux_ext_faces;

            DivideTetrahedra3D4::GenerateExteriorFaces(
                aux_ext_faces,
                aux_ext_faces_parent_ids,
                rSubdivisionsContainer,
                i_face);

            rExteriorFacesVector.insert(rExteriorFacesVector.end(), aux_ext_faces.begin(), aux_ext_faces.end());
            rExteriorFacesParentSubdivisionsIdsVector.insert(rExteriorFacesParentSubdivisionsIdsVector.end(), aux_ext_faces_parent_ids.begin(), aux_ext_faces_parent_ids.end());
        }
    };

    template<class TPointType>
    void DivideTetrahedra3D4<TPointType>::GenerateExteriorFaces(
        std::vector < IndexedPointGeometryPointerType > &rExteriorFacesVector,
        std::vector < unsigned int > &rExteriorFacesParentSubdivisionsIdsVector,
        const std::vector < IndexedPointGeometryPointerType > &rSubdivisionsContainer,
        const unsigned int FatherFaceId) {
        // Set some geometry constant parameters
        const unsigned int n_faces = 4;

        // Set the exterior faces vector
        rExteriorFacesVector.clear();
        rExteriorFacesParentSubdivisionsIdsVector.clear();

        if (this->mIsSplit) {
            // Create the face nodes data
            // The position represents the face while the value the real and intersection nodes in that face edges
            std::array < std::array< unsigned int, 6 >, 4> edges_map = {{
                {{2, 3, 1, 7, 8, 9}},     // Face 0
                {{0, 2, 3, 5, 6, 9}},     // Face 1
                {{0, 1, 3, 4, 6, 8}},     // Face 2
                {{0, 2, 1, 4, 5, 7}}}};   // Face 3

            // Compute the side exterior faces geometries
            const unsigned int n_subdivision = rSubdivisionsContainer.size();
            for (unsigned int i_subdivision = 0; i_subdivision < n_subdivision; ++i_subdivision) {
                // Get the subdivision faces
                const IndexedPointGeometryPointerType p_subdivision_geom = rSubdivisionsContainer[i_subdivision];
                const IndexedGeometriesArrayType subdivision_faces = p_subdivision_geom->GenerateFaces();

                // Subdivision geometry subfaces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    const IndexedPointGeometryType& r_face = subdivision_faces[i_face];

                    // Get the subdivision face nodal keys
                    int node_i_key = r_face[0].Id();
                    int node_j_key = r_face[1].Id();
                    int node_k_key = r_face[2].Id();

                    // Get the parent geometry face key value (candidate nodes)
                    std::array< unsigned int, 6 > faces_edge_nodes = edges_map[FatherFaceId];

                    // Search the subdivision nodal keys into the parent geometry face key value
                    if (std::find(faces_edge_nodes.begin(), faces_edge_nodes.end(), node_i_key) != faces_edge_nodes.end()) {
                        if (std::find(faces_edge_nodes.begin(), faces_edge_nodes.end(), node_j_key) != faces_edge_nodes.end()) {
                            if (std::find(faces_edge_nodes.begin(), faces_edge_nodes.end(), node_k_key) != faces_edge_nodes.end()) {
                                // If both nodes are in the candidate nodes list, the subface is exterior
                                IndexedPointGeometryPointerType p_subface_triang = Kratos::make_shared<IndexedPointTriangleType>(
                                    this->mAuxPointsContainer(node_i_key),
                                    this->mAuxPointsContainer(node_j_key),
                                    this->mAuxPointsContainer(node_k_key));
                                rExteriorFacesVector.push_back(p_subface_triang);
                                rExteriorFacesParentSubdivisionsIdsVector.push_back(i_subdivision);
                            }
                        }
                    }
                }
            }
        } else {
            KRATOS_ERROR << "Trying to generate the exterior faces in DivideTetrahedra3D4::GenerateExteriorFaces() for a non-split element.";
        }
    };

    template<class TPointType>
    bool DivideTetrahedra3D4<TPointType>::NodeIsInterface(int NodeKey) const
    {
        constexpr int num_nodes = 4;
        return NodeKey >= num_nodes || mNodeIsCut[NodeKey];
    }

    template class KRATOS_API(KRATOS_CORE) DivideTetrahedra3D4<Node>;
    template class KRATOS_API(KRATOS_CORE) DivideTetrahedra3D4<IndexedPoint>;

};
