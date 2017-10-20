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
    DivideTriangle2D3::DivideTriangle2D3(GeometryType& rInputGeometry, Vector& rNodalDistances) :
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

    // TODO: Save this arguments as member variables
    bool DivideTriangle2D3::GenerateDivision(IndexedPointsContainerType& rAuxPoints,
                                             std::vector < IndexedPointGeometryPointerType >& rPositiveSubdivisions,
                                             std::vector < IndexedPointGeometryPointerType >& rNegativeSubdivisions) {
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();

        // Fill the auxiliar points vector set
        if (DivideGeometry::IsSplit()) {

            // Set the triangle geometry parameters
            const unsigned int n_nodes = 3;
            const unsigned int n_edges = 3;

            // Clear the auxiliar vector points set
            rAuxPoints.clear();

            // Add the original geometry points
            for (unsigned int i = 0; i < n_nodes; ++i) {
                const array_1d<double, 3> aux_point_coords = geometry[i].Coordinates();
                IndexedPointPointerType paux_point(new IndexedPoint(aux_point_coords, i));
                rAuxPoints.push_back(paux_point);
            }

            // Decide the splitting pattern
            unsigned int aux_node_id = n_nodes;

            for(unsigned int idedge = 0; idedge < n_edges; ++idedge) {
                const unsigned int edge_node_i = (this->mEdgeNodeI)[idedge];
                const unsigned int edge_node_j = (this->mEdgeNodeJ)[idedge];

                // Check if the edge is split
                if(nodal_distances(edge_node_i) * nodal_distances(edge_node_j) < 0) {
                    // Set the new node id. in the split edge array corresponding slot
                    (this->mSplitEdges)[idedge + n_nodes] = aux_node_id;

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
                    IndexedPointPointerType paux_point(new IndexedPoint(aux_point_coords, aux_node_id));
                    rAuxPoints.push_back(paux_point);
                }

                aux_node_id++;
            }

            // Call the splitting mode computation function
            int edge_ids[3];
            TriangleSplitMode(this->mSplitEdges, edge_ids);

            // Call the splitting function
            int t[12];          // Ids of the generated subdivisions
            int n_int = 0;      // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            Split_Triangle(edge_ids, t, &(this->mDivisionsNumber), &(this->mSplitEdgesNumber), &n_int);

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < this->mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2;
                TriangleGetNewConnectivityGID(idivision, t, this->mSplitEdges, &i0, &i1, &i2);

                // Generate a pointer to an auxiliar triangular geometry made with the subdivision points
                IndexedPointGeometryPointerType p_aux_partition = boost::make_shared<IndexedPointTriangleType>(rAuxPoints(i0), rAuxPoints(i1), rAuxPoints(i2));

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
                    rPositiveSubdivisions.push_back(p_aux_partition);
                } else {
                    rNegativeSubdivisions.push_back(p_aux_partition);
                }
            }

            return true;

        } else {
            (this->mDivisionsNumber) = 1;
            (this->mSplitEdgesNumber) = 0;

            return false;
        }
    };

void DivideTriangle2D3::GenerateIntersectionsSkin(std::vector < IndexedPointGeometryPointerType >& rInterfacesVector,
                                                  IndexedPointsContainerType& rAuxPoints, // TODO: WHY THIS ARGUMENT CANNOT BE const ?¿?¿?¿?¿?¿?¿?¿?¿?¿?¿?¿¿?¿
                                                  const std::vector < IndexedPointGeometryPointerType >& rSubdivisionsVector) {
        // Set some geometry constant parameters
        const int n_nodes = 3;
        const unsigned int n_faces = 3;
        const unsigned int n_subdivision = rSubdivisionsVector.size();
        
        for (unsigned int i_subdivision = 0; i_subdivision < n_subdivision; ++i_subdivision) {
            // Get the subdivision geometry
            const IndexedPointGeometryType& r_subdivision_geom = *rSubdivisionsVector[i_subdivision];

            // Faces iteration
            for (unsigned int i_face = 0; i_face < n_faces; ++i_face) {
                // Get the subdivision face nodal keys
                int node_i_key = r_subdivision_geom[mEdgeNodeI[i_face]].Id();
                int node_j_key = r_subdivision_geom[mEdgeNodeJ[i_face]].Id();
                
                // Check the nodal keys to state which nodes belong to the interface
                // If the indexed keys is larger or equal to the number of nodes means that they are the auxiliar interface points
                if ((node_i_key >= n_nodes) && (node_j_key >= n_nodes)) {
                    // Generate an indexed point line geometry pointer with the two interface nodes
                    IndexedPointGeometryPointerType p_intersection_line = boost::make_shared<IndexedPointLineType>(rAuxPoints(node_i_key), rAuxPoints(node_j_key));
                    rInterfacesVector.push_back(p_intersection_line);

                    // In triangles, a unique face can belong to the interface
                    break; 
                }
            }
        }
    };
        
};
