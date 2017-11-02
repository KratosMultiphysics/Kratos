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

#include "utilities/split_triangle.c"
#include "utilities/divide_triangle_2d_3.h"

namespace Kratos
{
    /// Default constructor
    DivideTriangle2D3::DivideTriangle2D3(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
        DivideGeometry(rInputGeometry, rNodalDistances) {};

    /// Destructor
    DivideTriangle2D3::~DivideTriangle2D3() {};

    /// Turn back information as a string.
    std::string DivideTriangle2D3::Info() const {
        return "Triangle divide operations utility.";
    };

    /// Print information about this object.
    void DivideTriangle2D3::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Triangle divide operations utility.";
    };

    /// Print object's data.
    void DivideTriangle2D3::PrintData(std::ostream& rOStream) const {
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Triangle divide operations utility constructed with:\n";
        rOStream << "   Geometry type: " << geometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(nodal_distances(i)) << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    // Performs and saves the splitting pattern.
    void DivideTriangle2D3::GenerateDivision() {

        const GeometryType geometry = GetInputGeometry();
        const Vector nodal_distances = GetNodalDistances();

        // Fill the auxiliar points vector set
        if (mIsSplit) {

            // Set the triangle geometry parameters
            const unsigned int n_nodes = 3;
            const unsigned int n_edges = 3;

            // Clear the auxiliar vector points set
            mAuxPointsContainer.clear();
            mAuxPointsContainer.reserve(6);

            // Clear the subdivision vectors
            mPositiveSubdivisions.clear();
            mNegativeSubdivisions.clear();
            mPositiveSubdivisions.reserve(2);
            mNegativeSubdivisions.reserve(2);

            // Add the original geometry points
            for (unsigned int i = 0; i < n_nodes; ++i) {
                const array_1d<double, 3> aux_point_coords = geometry[i].Coordinates();
                IndexedPointPointerType paux_point = boost::make_shared<IndexedPoint>(aux_point_coords, i);
                mAuxPointsContainer.push_back(paux_point);
            }

            // Decide the splitting pattern
            unsigned int aux_node_id = n_nodes;

            for(unsigned int idedge = 0; idedge < n_edges; ++idedge) {
                const unsigned int edge_node_i = mEdgeNodeI[idedge];
                const unsigned int edge_node_j = mEdgeNodeJ[idedge];

                // Check if the edge is split
                if(nodal_distances(edge_node_i) * nodal_distances(edge_node_j) < 0) {
                    // Set the new node id. in the split edge array corresponding slot
                    mSplitEdges[idedge + n_nodes] = aux_node_id;

                    // Edge nodes coordinates
                    const array_1d<double, 3> i_node_coords = geometry[edge_node_i].Coordinates();
                    const array_1d<double, 3> j_node_coords = geometry[edge_node_j].Coordinates();

                    // Compute the coordinates of the point on the edge
                    const double aux_node_rel_location = std::abs (nodal_distances(edge_node_i)/(nodal_distances(edge_node_j)-nodal_distances(edge_node_i)));
                    array_1d<double, 3> aux_point_coords;
                    for (unsigned int comp = 0; comp < 3; ++comp) {
                        aux_point_coords(comp) = i_node_coords(comp)*aux_node_rel_location + j_node_coords(comp)*(1.0-aux_node_rel_location);
                    }

                    // Add the intersection point to the auxiliar points array
                    IndexedPointPointerType paux_point = boost::make_shared<IndexedPoint>(aux_point_coords, aux_node_id);
                    mAuxPointsContainer.push_back(paux_point);
                }

                aux_node_id++;
            }

            // Call the splitting mode computation function
            std::vector<int> edge_ids(3);
            TriangleSplitMode(mSplitEdges.data(), edge_ids.data());

            // Call the splitting function
            std::vector<int> t(12);     // Ids of the generated subdivisions
            int n_int = 0;              // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            Split_Triangle(edge_ids.data(), t.data(), &mDivisionsNumber, &mSplitEdgesNumber, &n_int);

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2;
                TriangleGetNewConnectivityGID(idivision, t.data(), mSplitEdges.data(), &i0, &i1, &i2);

                // Generate a pointer to an auxiliar triangular geometry made with the subdivision points
                IndexedPointGeometryPointerType p_aux_partition = boost::make_shared<IndexedPointTriangleType>(mAuxPointsContainer(i0), 
                                                                                                               mAuxPointsContainer(i1), 
                                                                                                               mAuxPointsContainer(i2));

                // Determine if the subdivision is wether in the negative or the positive side
                bool is_positive;
                if ((i0 == 0) || (i0 == 1) || (i0 == 2)) {
                    is_positive = nodal_distances(i0) < 0.0 ? false : true;
                } else if ((i1 == 0) || (i1 == 1) || (i1 == 2)) {
                    is_positive = nodal_distances(i1) < 0.0 ? false : true;
                } else {
                    is_positive = nodal_distances(i2) < 0.0 ? false : true;
                }

                // Add the generated triangle to its corresponding partition arrays
                if (is_positive) {
                    mPositiveSubdivisions.push_back(p_aux_partition);
                } else {
                    mNegativeSubdivisions.push_back(p_aux_partition);
                }
            }

        } else {
            mDivisionsNumber = 1;
            mSplitEdgesNumber = 0;
        }
    };

void DivideTriangle2D3::GenerateIntersectionsSkin() {
        
        // Set some geometry constant parameters
        const int n_nodes = 3;
        const unsigned int n_faces = 3;

        // Clear the interfaces vectors
        mPositiveInterfaces.clear();
        mNegativeInterfaces.clear();
        mPositiveInterfaces.reserve(1);
        mNegativeInterfaces.reserve(1);

        if (mIsSplit) {

            const unsigned int n_positive_subdivision = mPositiveSubdivisions.size();
            const unsigned int n_negative_subdivision = mNegativeSubdivisions.size();
            
            // Compute the positive side intersection geometries
            for (unsigned int i_subdivision = 0; i_subdivision < n_positive_subdivision; ++i_subdivision) {
                // Get the subdivision geometry
                const IndexedPointGeometryType& r_subdivision_geom = *mPositiveSubdivisions[i_subdivision];

                // Faces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    // Get the subdivision face nodal keys
                    int node_i_key = r_subdivision_geom[mEdgeNodeI[i_face]].Id();
                    int node_j_key = r_subdivision_geom[mEdgeNodeJ[i_face]].Id();
                    
                    // Check the nodal keys to state which nodes belong to the interface
                    // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliar interface points
                    if ((node_i_key >= n_nodes) && (node_j_key >= n_nodes)) {
                        // Generate an indexed point line geometry pointer with the two interface nodes
                        IndexedPointGeometryPointerType p_intersection_line = boost::make_shared<IndexedPointLineType>(mAuxPointsContainer(node_i_key), 
                                                                                                                       mAuxPointsContainer(node_j_key));
                        mPositiveInterfaces.push_back(p_intersection_line);

                        // In triangles, a unique face can belong to the interface
                        break; 
                    }
                }
            }

            // Compute the negative side intersection geometries
            for (unsigned int i_subdivision = 0; i_subdivision < n_negative_subdivision; ++i_subdivision) {
                // Get the subdivision geometry
                const IndexedPointGeometryType& r_subdivision_geom = *mNegativeSubdivisions[i_subdivision];

                // Faces iteration
                for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                    // Get the subdivision face nodal keys
                    int node_i_key = r_subdivision_geom[mEdgeNodeI[i_face]].Id();
                    int node_j_key = r_subdivision_geom[mEdgeNodeJ[i_face]].Id();
                    
                    // Check the nodal keys to state which nodes belong to the interface
                    // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliar interface points
                    if ((node_i_key >= n_nodes) && (node_j_key >= n_nodes)) {
                        // Generate an indexed point line geometry pointer with the two interface nodes
                        IndexedPointGeometryPointerType p_intersection_line = boost::make_shared<IndexedPointLineType>(mAuxPointsContainer(node_i_key), 
                                                                                                                       mAuxPointsContainer(node_j_key));
                        mNegativeInterfaces.push_back(p_intersection_line);

                        // In triangles, a unique face can belong to the interface
                        break; 
                    }
                }
            }
        } else {
            KRATOS_ERROR << "Trying to generate the intersection skin in DivideTriangle2D3::GenerateIntersectionsSkin() for a non-split element.";
        }
    };
        
};
