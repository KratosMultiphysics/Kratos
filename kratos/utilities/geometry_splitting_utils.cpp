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
    GeometrySplittingUtils::GeometrySplittingUtils(const GeometryType& rInputGeometry, const Vector& rNodalDistances) :
        mrInputGeometry(rInputGeometry),
        mrNodalDistances(rNodalDistances) {
    };

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
        KRATOS_WATCH("Inside GetInputGeometry")
        KRATOS_WATCH(mrInputGeometry)
        return mrInputGeometry;
    };

    Vector GeometrySplittingUtils::GetNodalDistances() const {
        KRATOS_WATCH("Inside GetNodalDistances")
        return mrNodalDistances;
    };

    bool GeometrySplittingUtils::DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                                std::vector < PointGeometryType >& rPositiveSubdivisions,
                                                std::vector < PointGeometryType >& rNegativeSubdivisions) {
        KRATOS_ERROR << "Calling the base class geometry splitting DivideGeometry method. Call the specific geometry one.";
    };

    bool GeometrySplittingUtils::IsSplit() {
        KRATOS_WATCH("IsSplit")
        unsigned int n_pos = 0 , n_neg = 0;
        KRATOS_WATCH(mrNodalDistances)
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
                                                                         const unsigned int splitEdgesNumber) {
        // Initialize intersection points condensation matrix
        const unsigned int nnodes_original = mrInputGeometry.size();
        const unsigned int total_nodes = nnodes_original + splitEdgesNumber;

        rIntPointCondMatrix = ZeroMatrix(total_nodes, nnodes_original);

        // Fill the original geometry points main diagonal
        for (unsigned int i = 0; i < nnodes_original; ++i) {
            rIntPointCondMatrix(i,i) = 1.0;
        }

        KRATOS_WATCH(rIntPointCondMatrix)

        // Compute the intersection points contributions
        for (unsigned int intersection_id = 0; intersection_id < splitEdgesNumber; ++intersection_id) {
            const unsigned int aux_count = intersection_id + nnodes_original;
            KRATOS_WATCH(aux_count)
            KRATOS_WATCH(mSplitEdges[0])
            KRATOS_WATCH(mSplitEdges[1])
            KRATOS_WATCH(mSplitEdges[2])
            KRATOS_WATCH(mSplitEdges[3])
            KRATOS_WATCH(mSplitEdges[4])
            KRATOS_WATCH(mSplitEdges[5])
            // Check if the edge has an intersection point
            if (mSplitEdges[aux_count] != -1) {
                // Get the nodes that compose the edge
                const unsigned int edge_node_i = mEdgeNodeI[intersection_id];
                const unsigned int edge_node_j = mEdgeNodeJ[intersection_id];

                // Compute the relative coordinate of the intersection point over the edge
                const double aux_node_rel_location = std::abs (mrNodalDistances(edge_node_i)/(mrNodalDistances(edge_node_j)-mrNodalDistances(edge_node_i)));

                KRATOS_WATCH(aux_count)
                KRATOS_WATCH(edge_node_i)
                KRATOS_WATCH(edge_node_j)

                // Store the relative coordinate values as the original geometry nodes sh. function value in the intersections
                rIntPointCondMatrix(aux_count, edge_node_i) = aux_node_rel_location;
                rIntPointCondMatrix(aux_count, edge_node_j) = 1.0 - aux_node_rel_location;
            }
        }
    }

    /// TriangleSplittingUtils implementation
    /// Default constructor
    TriangleSplittingUtils::TriangleSplittingUtils(const GeometryType& rInputGeometry, const array_1d< double, 3 >& rNodalDistances) :
        GeometrySplittingUtils(rInputGeometry, rNodalDistances) {

        // Initialize the base class member variables
        int aux;
        for (unsigned int i = 0; i < rInputGeometry.size(); ++i) {
            KRATOS_WATCH(i)
            this->mEdgeNodeI[i] = i;

            if (i != 2) {
                aux = i + 1;
                KRATOS_WATCH("i!=2")
                KRATOS_WATCH(aux)
            } else {
                aux = 0;
                KRATOS_WATCH("i==2")
                KRATOS_WATCH(aux)
            }
            this->mEdgeNodeJ[i] = aux;
        }

        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        // TODO: BUG IN HERE!!! THE EDGE VALUES ARE NOT CORRECT!!!!!!!!!!!
        KRATOS_WATCH(this->mEdgeNodeI[0])
        KRATOS_WATCH(this->mEdgeNodeI[1])
        KRATOS_WATCH(this->mEdgeNodeI[2])
        KRATOS_WATCH(this->mEdgeNodeJ[0])
        KRATOS_WATCH(this->mEdgeNodeJ[1])
        KRATOS_WATCH(this->mEdgeNodeJ[2])
        KRATOS_WATCH("Edge nodes set")

        const unsigned int nedge = 6;
        for (unsigned int i = 0; i < nedge; ++i) {
            (this->mSplitEdges)[i] = i > 2 ? -1 : i; // {0, 1, 2, -1, -1, -1} what is to say {node_0, node_1, node_2, node_edge_01, node_edge_12, node_edge_20}
        }

        KRATOS_WATCH((this->mSplitEdges)[0])
        KRATOS_WATCH((this->mSplitEdges)[1])
        KRATOS_WATCH((this->mSplitEdges)[2])
        KRATOS_WATCH((this->mSplitEdges)[3])
        KRATOS_WATCH((this->mSplitEdges)[4])
        KRATOS_WATCH((this->mSplitEdges)[5])
        KRATOS_WATCH("Split edges set")

        KRATOS_WATCH("Exiting derived class constructor")
    };

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
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Triangle splitting operations utility constructed with:\n";
        rOStream << "   Geometry type: " << geometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(nodal_distances(i)) << " ";
        }
        rOStream << "   Distance values: " << distances_buffer.str();
    };

    bool TriangleSplittingUtils::DivideGeometry(IndexedPointsContainerType& rAuxPoints,
                                                std::vector < PointGeometryType >& rPositiveSubdivisions,
                                                std::vector < PointGeometryType >& rNegativeSubdivisions) {

        KRATOS_WATCH("Inside DivideGeometry")
        GeometryType geometry = this->GetInputGeometry();
        KRATOS_WATCH("GetInputGeometry() works!!!")
        const Vector nodal_distances = this->GetNodalDistances();
        KRATOS_WATCH(nodal_distances)
        KRATOS_WATCH("GetNodalDistances() works!!!")
        // Fill the auxiliar points vector set
        if (GeometrySplittingUtils::IsSplit()) {

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

            KRATOS_WATCH((this->mSplitEdges)[0])
            KRATOS_WATCH((this->mSplitEdges)[1])
            KRATOS_WATCH((this->mSplitEdges)[2])
            KRATOS_WATCH((this->mSplitEdges)[3])
            KRATOS_WATCH((this->mSplitEdges)[4])
            KRATOS_WATCH((this->mSplitEdges)[5])
            KRATOS_WATCH(edge_ids[0])
            KRATOS_WATCH(edge_ids[1])
            KRATOS_WATCH(edge_ids[2])

            // Call the splitting function
            int t[12];          // Ids of the generated subdivisions
            int n_int = 0;      // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            Split_Triangle(edge_ids, t, &(this->mDivisionsNumber), &(this->mSplitEdgesNumber), &n_int);

            for (unsigned int i=0; i<12; ++i) {KRATOS_WATCH(t[i]);}
            KRATOS_WATCH(this->mDivisionsNumber)
            KRATOS_WATCH(this->mSplitEdgesNumber)
            KRATOS_WATCH(n_int)

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < this->mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2;
                TriangleGetNewConnectivityGID(idivision, t, this->mSplitEdges, &i0, &i1, &i2);

                std::cout << "Subdivision " << idivision << " i0: " << i0 << " i1: " << i1 << " i2: " << i2 << std::endl;

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
                    is_positive = nodal_distances(i0) < 0.0 ? false : true;
                } else if ((i1 == 0) || (i1 == 1) || (i1 == 2)) {
                    is_positive = nodal_distances(i1) < 0.0 ? false : true;
                } else {
                    is_positive = nodal_distances(i2) < 0.0 ? false : true;
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
            SetIntersectionPointsCondensationMatrix(p_matrix, this->mSplitEdgesNumber);
            KRATOS_WATCH(p_matrix)

            return true;

        } else {

            return false;
        }
    };

};
