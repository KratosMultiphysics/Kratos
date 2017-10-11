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
#include "utilities/geometry_splitting_utils.h"

namespace Kratos
{
    /// IndexedPoint class implementation
    /// Constructors
    IndexedPoint::IndexedPoint()
        : Point<3>() , IndexedObject(0) {};

    IndexedPoint::IndexedPoint(const array_1d<double,3>& rCoords, const unsigned int Id)
        : Point<3>(rCoords) , IndexedObject(Id) {};

    /// Destructor
    IndexedPoint::~IndexedPoint() {};

    /// Turn back information as a string.
    std::string IndexedPoint::Info() const {
        std::stringstream info_string;
        info_string << "Indexed point class created as a combination of the Point<3> and IndexedObject classes.\n";
        info_string << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
        return info_string.str();
    };

    /// Print information about this object.
    void IndexedPoint::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Indexed point class created as a combination of the Point<3> and IndexedObject classes.\n";
        rOStream << "This indexed point class is intended to be used to define both the real and auxiliar points of an intersected element.";
    };

    /// Print object's data.
    void IndexedPoint::PrintData(std::ostream& rOStream) const {
        rOStream << "Indexed point object constructed with:\n";
        rOStream << "   Index value: " << this->Id() << "\n";
        const array_1d<double, 3> point_coords = this->Coordinates();
        std::stringstream coordinates_buffer;
        for (unsigned int i = 0; i < 3; ++i) {
            coordinates_buffer << std::to_string(point_coords(i)) << " ";
        }
        rOStream << "   Coordinates: " << coordinates_buffer.str();
    };

    /// GeometrySplittingUtils implementation
    /// Default constructor
    GeometrySplittingUtils::GeometrySplittingUtils(GeometryType& rInputGeometry, Vector& rNodalDistances) :
        mrInputGeometry(rInputGeometry), mrNodalDistances(rNodalDistances) {};

    /// Destructor
    GeometrySplittingUtils::~GeometrySplittingUtils() {};

    /// Turn back information as a string.
    std::string GeometrySplittingUtils::Info() const {
        return "Base class for geometries splitting operations.";
    };

    /// Print information about this object.
    void GeometrySplittingUtils::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations.";
    };

    /// Print object's data.
    void GeometrySplittingUtils::PrintData(std::ostream& rOStream) const {
        rOStream << "Base class for geometries splitting operations constructed with:\n";
        rOStream << "   Geometry type: " << mrInputGeometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            distances_buffer << std::to_string(mrNodalDistances(i)) << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    GeometrySplittingUtils::GeometryType GeometrySplittingUtils::GetInputGeometry() const {
        return mrInputGeometry;
    };

    Vector GeometrySplittingUtils::GetNodalDistances() const {
        return mrNodalDistances;
    };

    bool GeometrySplittingUtils::DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                                std::vector < PointGeometryType >& rPositiveSubdivisions,
                                                std::vector < PointGeometryType >& rNegativeSubdivisions) {
        KRATOS_ERROR << "Calling the base class geometry splitting DivideGeometry method. Call the specific geometry one.";
    };

    bool GeometrySplittingUtils::IsSplit() {
        unsigned int n_pos = 0 , n_neg = 0;
        
        for (unsigned int i = 0; i < mrNodalDistances.size(); ++i) {
            if (mrNodalDistances(i) < 0.0) {
                n_neg++;
            } else {
                n_pos++;
            }
        }

        if ((n_pos > 0) && (n_neg > 0)) {
            return true;
        } else {
            return false;
        }
    };

    void GeometrySplittingUtils::SetIntersectionPointsCondensationMatrix(Matrix& rIntPointCondMatrix,
                                                                         const int rEdgeNodeI[], 
                                                                         const int rEdgeNodeJ[], 
                                                                         const int rSplitEdges[], 
                                                                         const unsigned int splitEdgesNumber) {
        KRATOS_WATCH("Inside SetIntersectionPointsCondensationMatrix")
        KRATOS_WATCH(mrInputGeometry)
        // Initialize intersection points condensation matrix
        // const unsigned int nnodes = mrInputGeometry.PointsNumber();
        // const unsigned int nnodes = mrInputGeometry.size();
        const unsigned int nnodes = 3;
        KRATOS_WATCH(nnodes)
        // const unsigned int nedges = mrInputGeometry.EdgesNumber();
        const unsigned int nedges = 3;
        KRATOS_WATCH(nedges)
        const unsigned int total_nodes = nnodes + splitEdgesNumber;
        KRATOS_WATCH(total_nodes)


        rIntPointCondMatrix = ZeroMatrix(total_nodes, nnodes);

        // Fill the original geometry points main diagonal
        for (unsigned int i = 0; i < nnodes; ++i) {
            rIntPointCondMatrix(i,i) = 1.0;
        }

        // Compute the intersection points contributions
        unsigned int row_count = nnodes;
        for (unsigned int edge_id = 0; edge_id < nedges; ++edge_id) {
            const unsigned int aux_count = edge_id + nnodes;
            KRATOS_WATCH("")
            KRATOS_WATCH(aux_count)
            // Check if the edge has an intersection point
            if (rSplitEdges[aux_count] != -1) {
                // Get the nodes that compose the edge
                const unsigned int edge_node_i = rEdgeNodeI[edge_id];
                const unsigned int edge_node_j = rEdgeNodeJ[edge_id];

                // Compute the relative coordinate of the intersection point over the edge
                const double aux_node_rel_location = std::abs (mrNodalDistances(edge_node_i)/(mrNodalDistances(edge_node_j)-mrNodalDistances(edge_node_i)));

                KRATOS_WATCH(aux_count)
                KRATOS_WATCH(edge_node_i)
                KRATOS_WATCH(edge_node_j)

                // Store the relative coordinate values as the original geometry nodes sh. function value in the intersections
                rIntPointCondMatrix(row_count, edge_node_i) = aux_node_rel_location;
                rIntPointCondMatrix(row_count, edge_node_j) = 1.0 - aux_node_rel_location;

                row_count++;
            }
        }
    }

    /// TriangleSplittingUtils implementation
    /// Default constructor
    TriangleSplittingUtils::TriangleSplittingUtils(GeometryType& rInputGeometry, Vector& rNodalDistances) :
        GeometrySplittingUtils(rInputGeometry, rNodalDistances) {};

    /// Destructor
    TriangleSplittingUtils::~TriangleSplittingUtils() {};

    /// Turn back information as a string.
    std::string TriangleSplittingUtils::Info() const {
        return "Triangle splitting operations utility.";
    };

    /// Print information about this object.
    void TriangleSplittingUtils::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Triangle splitting operations utility.";
    };

    /// Print object's data.
    void TriangleSplittingUtils::PrintData(std::ostream& rOStream) const {
        GeometryType r_geometry = this->GetInputGeometry();
        Vector r_nodal_distances = this->GetNodalDistances();
        rOStream << "Triangle splitting operations utility constructed with:\n";
        rOStream << "   Geometry type: " << r_geometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < r_nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(r_nodal_distances(i)) << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    bool TriangleSplittingUtils::DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                                std::vector < PointGeometryType >& rPositiveSubdivisions,
                                                std::vector < PointGeometryType >& rNegativeSubdivisions) {

        Vector r_nodal_distances = this->GetNodalDistances();
        GeometryType r_geometry = this->GetInputGeometry();
        // Fill the auxiliar points vector set
        if (GeometrySplittingUtils::IsSplit()) {

            // Set the triangle geometry parameters
            const unsigned int n_nodes = 3;
            const unsigned int n_edges = 3;

            // Clear the auxiliar vector points set
            rAuxPoints.clear();

            // Add the original geometry points
            for (unsigned int i = 0; i < n_nodes; ++i) {
                const array_1d<double, 3> aux_point_coords = r_geometry[i].Coordinates();
                IndexedPointPointerType paux_point(new IndexedPoint(aux_point_coords, i));
                rAuxPoints.push_back(paux_point);
            }

            // Decide the splitting pattern
            unsigned int aux_node_id = n_nodes;

            for(unsigned int idedge = 0; idedge < n_edges; ++idedge) {
                const unsigned int edge_node_i = mEdgeNodeI[idedge];
                const unsigned int edge_node_j = mEdgeNodeJ[idedge];

                // Check if the edge is split
                if(r_nodal_distances(edge_node_i) * r_nodal_distances(edge_node_j) < 0) {
                    // Set the new node id. in the split edge array corresponding slot
                    mSplitEdges[idedge + n_nodes] = aux_node_id;

                    // Edge nodes coordinates
                    const array_1d<double, 3> i_node_coords = r_geometry[edge_node_i].Coordinates();
                    const array_1d<double, 3> j_node_coords = r_geometry[edge_node_j].Coordinates();

                    // Compute the coordinates of the point on the edge
                    const double aux_node_rel_location = std::abs (r_nodal_distances(edge_node_i)/(r_nodal_distances(edge_node_j)-r_nodal_distances(edge_node_i)));
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
            TriangleSplitMode(mSplitEdges, edge_ids);

            // Call the splitting function
            int t[12];          // Ids of the generated subdivisions
            int n_int = 0;      // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            Split_Triangle(edge_ids, t, &(this->mDivisionsNumber), &(this->mSplitEdgesNumber), &n_int);

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < this->mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2;
                TriangleGetNewConnectivityGID(idivision, t, mSplitEdges, &i0, &i1, &i2);

                // Get a pointer to the point objects corresponding to the indices above
                boost::shared_ptr< Point<3> > p_point_0(new Point<3>((rAuxPoints.find(i0))->Coordinates()));
                boost::shared_ptr< Point<3> > p_point_1(new Point<3>((rAuxPoints.find(i1))->Coordinates()));
                boost::shared_ptr< Point<3> > p_point_2(new Point<3>((rAuxPoints.find(i2))->Coordinates()));

                // Generate a triangle (with the same type as the input one but generated with points) with the subdivision pts.
	            Geometry< Point<3> >::PointsArrayType PointsArray;
	            PointsArray.reserve(3);
                PointsArray.push_back(p_point_0);
                PointsArray.push_back(p_point_1);
                PointsArray.push_back(p_point_2);
                PointGeometryType aux_partition(PointsArray);

                // Determine if the subdivision is wether in the negative or the positive side
                bool is_positive;
                if ((i0 == 0) || (i0 == 1) || (i0 == 2)) {
                    is_positive = r_nodal_distances(i0) < 0.0 ? false : true;
                } else if ((i1 == 0) || (i1 == 1) || (i1 == 2)) {
                    is_positive = r_nodal_distances(i1) < 0.0 ? false : true;
                } else {
                    is_positive = r_nodal_distances(i2) < 0.0 ? false : true;
                }

                // Add the generated triangle to its corresponding partition arrays
                if (is_positive) {
                    rPositiveSubdivisions.push_back(aux_partition);
                } else {
                    rNegativeSubdivisions.push_back(aux_partition);
                }
            }

            // TODO: TEMPORARY HERE FOR DEBUGGING. MOVE TO THE SHAPE FUNCTIONS COMPUTATION SECTION ASAP
            Matrix p_matrix;
            SetIntersectionPointsCondensationMatrix(p_matrix, mEdgeNodeI, mEdgeNodeJ, mSplitEdges, this->mSplitEdgesNumber);
            KRATOS_WATCH(p_matrix)

            return true;

        } else {

            return false;
        }
    };

};
