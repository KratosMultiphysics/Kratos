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
#include "geometries/triangle_2d_3.h"
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

    void GeometrySplittingUtils::GetSubdivisionShapeFunctionValues(Matrix& rShapeFunctionValues,
                                                                   const PointGeometryType& rSubdivisionGeom,
                                                                   IntegrationMethodType IntegrationMethod) {
        KRATOS_ERROR << "Calling the base class geometry splitting GetSubdivisionShapeFunctionValues method. Call the specific geometry one.";
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
                                                                         const int mEdgeNodeI[],
                                                                         const int mEdgeNodeJ[],
                                                                         const int mSplitEdges[],
                                                                         const unsigned int splitEdgesNumber) {
        // Initialize intersection points condensation matrix
        const unsigned int nedges = mrInputGeometry.EdgesNumber();
        const unsigned int nnodes = mrInputGeometry.PointsNumber();
        const unsigned int nnodes_total = nnodes + splitEdgesNumber;

        rIntPointCondMatrix = ZeroMatrix(nnodes_total, nnodes);

        // Fill the original geometry points main diagonal
        for (unsigned int i = 0; i < nnodes; ++i) {
            rIntPointCondMatrix(i,i) = 1.0;
        }

        // Compute the intersection points contributions
        unsigned int row = nnodes;
        for (unsigned int idedge = 0; idedge < nedges; ++idedge) {
            // Check if the edge has an intersection point
            if (mSplitEdges[nnodes+idedge] != -1) {
                // Get the nodes that compose the edge
                const unsigned int edge_node_i = mEdgeNodeI[idedge];
                const unsigned int edge_node_j = mEdgeNodeJ[idedge];

                // Compute the relative coordinate of the intersection point over the edge
                const double aux_node_rel_location = std::abs (mrNodalDistances(edge_node_i)/(mrNodalDistances(edge_node_j)-mrNodalDistances(edge_node_i)));

                // Store the relative coordinate values as the original geometry nodes sh. function value in the intersections
                rIntPointCondMatrix(row, edge_node_i) = aux_node_rel_location;
                rIntPointCondMatrix(row, edge_node_j) = 1.0 - aux_node_rel_location;
                row++;
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
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();

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

            // Call the splitting function
            int t[12];          // Ids of the generated subdivisions
            int n_int = 0;      // Number of internal nodes (set to 0 since it is not needed for triangle splitting)
            Split_Triangle(edge_ids, t, &(this->mDivisionsNumber), &(this->mSplitEdgesNumber), &n_int);

            // Fill the subdivisions arrays
            for (int idivision = 0; idivision < this->mDivisionsNumber; ++idivision) {
                // Get the subdivision indices
                int i0, i1, i2;
                TriangleGetNewConnectivityGID(idivision, t, this->mSplitEdges, &i0, &i1, &i2);

                // Get a pointer to the point objects corresponding to the indices above
                boost::shared_ptr< Point<3> > p_point_0(new Point<3>((rAuxPoints.find(i0))->Coordinates()));
                boost::shared_ptr< Point<3> > p_point_1(new Point<3>((rAuxPoints.find(i1))->Coordinates()));
                boost::shared_ptr< Point<3> > p_point_2(new Point<3>((rAuxPoints.find(i2))->Coordinates()));

                // Generate a triangle (with the same type as the input one but generated with points) with the subdivision pts.
	            PointGeometryType::PointsArrayType points_array;
	            points_array.reserve(3);
                points_array.push_back(p_point_0);
                points_array.push_back(p_point_1);
                points_array.push_back(p_point_2);
                PointGeometryType aux_partition(points_array);

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
            SetIntersectionPointsCondensationMatrix(p_matrix, mEdgeNodeI, mEdgeNodeJ, mSplitEdges, this->mSplitEdgesNumber);

            return true;

        } else {
            (this->mDivisionsNumber) = 1;
            (this->mSplitEdgesNumber) = 0;

            return false;
        }
    };

    void TriangleSplittingUtils::GetSubdivisionShapeFunctionValues(Matrix& rShapeFunctionValues,
                                                                   const PointGeometryType& rSubdivisionGeom,
                                                                   IntegrationMethodType IntegrationMethod) {
        // Get the integration points
        IntegrationPointsContainerType all_integration_points = Triangle2D3<PointType>::AllIntegrationPoints(); // TODO: DECIDE WHAT TO DO WITH THE TRIANGLE STATIC METHODS, THEY ARE PRIVATE
        IntegrationPointsArrayType integration_points = all_integration_points[IntegrationMethod];

        // Resize the shape function values matrix
        const unsigned int n_nodes = rSubdivisionGeom.PointsNumber(); // TODO: THIS SHOULD HAVE THE MAX NUMBER OF NODES SIZE
        const unsigned int n_int_pts = integration_points.IntegrationPointsNumber();
        //
        // if(rShapeFunctionValues.size1() != n_int_pts) {
        //     rShapeFunctionValues.resize(n_int_pts, n_nodes, false);
        // } else if(rShapeFunctionValues.size2() != n_nodes) {
        //     rShapeFunctionValues.resize(n_int_pts, n_nodes, false);
        // }
        //
        // // Compute the subdivision shape functions values
        // for (unsigned int i_int_pt = 0; i_int_pt < n_int_pts; ++i_int_pt) {
        //     Vector point_sh_func_values;
        //     CoordinatesArrayType int_pt_local_coords;
        //     int_pt_local_coords[0] = integration_points[i_int_pt].X()
        //     int_pt_local_coords[1] = integration_points[i_int_pt].Y()
        //     int_pt_local_coords[2] = integration_points[i_int_pt].Z()
        //     rSubdivisionGeom.ShapeFunctionsValues(point_sh_func_values, int_pt_local_coords);
        // }

    };

};
