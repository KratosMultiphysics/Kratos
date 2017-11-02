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
#include "modified_shape_functions/modified_shape_functions.h"

namespace Kratos
{

    /// ModifiedShapeFunctions implementation
    /// Default constructor
    ModifiedShapeFunctions::ModifiedShapeFunctions(const GeometryPointerType pInputGeometry, const Vector& rNodalDistances) :
        mpInputGeometry(pInputGeometry), mrNodalDistances(rNodalDistances) {
    };

    /// Destructor
    ModifiedShapeFunctions::~ModifiedShapeFunctions() {};

    /// Turn back information as a string.
    std::string ModifiedShapeFunctions::Info() const {
        return "Modified shape functions computation base class.";
    };

    /// Print information about this object.
    void ModifiedShapeFunctions::PrintInfo(std::ostream& rOStream) const {
        rOStream << "Modified shape functions computation base class.";
    };

    /// Print object's data.
    void ModifiedShapeFunctions::PrintData(std::ostream& rOStream) const {
        const GeometryPointerType p_geometry = this->GetInputGeometry();
        const Vector& nodal_distances = this->GetNodalDistances();
        rOStream << "Modified shape functions computation base class:\n";
        rOStream << "\tGeometry type: " << (*p_geometry).Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(nodal_distances(i)) << " ";
        }
        rOStream << "\tDistance values: " << distances_buffer.str();
    };

    // Returns the input original geometry.
    const ModifiedShapeFunctions::GeometryPointerType ModifiedShapeFunctions::GetInputGeometry() const {
        return mpInputGeometry;
    };

    // Returns the nodal distances vector.
    const Vector& ModifiedShapeFunctions::GetNodalDistances() const {
        return mrNodalDistances;
    };

    // Sets the condensation matrix to transform the subdivsion values to entire element ones.
    void ModifiedShapeFunctions::SetIntersectionPointsCondensationMatrix(
        Matrix& rIntPointCondMatrix,
        const std::vector<int>& rEdgeNodeI,
        const std::vector<int>& rEdgeNodeJ,
        const std::vector<int>& rSplitEdges) {

        const unsigned int nedges = mpInputGeometry->EdgesNumber();
        const unsigned int nnodes = mpInputGeometry->PointsNumber();
            
        // Initialize intersection points condensation matrix
        rIntPointCondMatrix = ZeroMatrix(rSplitEdges.size(), nnodes);

        // Fill the original geometry points main diagonal
        for (unsigned int i = 0; i < nnodes; ++i) {
            rIntPointCondMatrix(i,i) = 1.0;
        }

        // Compute the intersection points contributions
        unsigned int row = nnodes;
        for (unsigned int idedge = 0; idedge < nedges; ++idedge) {
            // Check if the edge has an intersection point
            if (rSplitEdges[nnodes+idedge] != -1) {
                // Get the nodes that compose the edge
                const unsigned int edge_node_i = rEdgeNodeI[idedge];
                const unsigned int edge_node_j = rEdgeNodeJ[idedge];

                // Compute the relative coordinate of the intersection point over the edge
                const double aux_node_rel_location = std::abs (mrNodalDistances(edge_node_i)/(mrNodalDistances(edge_node_j)-mrNodalDistances(edge_node_i)));

                // Store the relative coordinate values as the original geometry nodes sh. function value in the intersections
                rIntPointCondMatrix(row, edge_node_i) = aux_node_rel_location;
                rIntPointCondMatrix(row, edge_node_j) = 1.0 - aux_node_rel_location;
            }
            row++;
        }
    }

    // Given the subdivision pattern of either the positive or negative side, computes the shape function values.
    void ModifiedShapeFunctions::ComputeValuesOnOneSide(
        Matrix &rShapeFunctionsValues,
        std::vector<Matrix> &rShapeFunctionsGradientsValues,
        Vector &rWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
        const Matrix &rPmatrix,
        const IntegrationMethodType IntegrationMethod) {

        // Set some auxiliar constants
        GeometryPointerType p_input_geometry = this->GetInputGeometry();
        const unsigned int n_edges_global = p_input_geometry->EdgesNumber();       // Number of edges in original geometry
        const unsigned int n_nodes_global = p_input_geometry->PointsNumber();      // Number of nodes in original geometry
        const unsigned int split_edges_size = n_edges_global + n_nodes_global;     // Split edges vector size

        const unsigned int n_subdivision = rSubdivisionsVector.size();                                         // Number of positive or negative subdivisions
        const unsigned int n_dim = (*rSubdivisionsVector[0]).Dimension();                                      // Number of dimensions
        const unsigned int n_nodes = (*rSubdivisionsVector[0]).PointsNumber();                                 // Number of nodes per subdivision
        const unsigned int n_int_pts = (*rSubdivisionsVector[0]).IntegrationPointsNumber(IntegrationMethod);   // Number of Gauss pts. per subdivision

        // Resize the shape function values matrix
        const unsigned int n_total_int_pts = n_subdivision * n_int_pts;

        if (rShapeFunctionsValues.size1() != n_total_int_pts) {
            rShapeFunctionsValues.resize(n_total_int_pts, n_nodes, false);
        } else if (rShapeFunctionsValues.size2() != n_nodes) {
            rShapeFunctionsValues.resize(n_total_int_pts, n_nodes, false);
        }

        // Resize the weights vector
        if (rWeightsValues.size() != n_total_int_pts) {
            rWeightsValues.resize(n_total_int_pts, false);
        }

        // Clear the gradients vector
        rShapeFunctionsGradientsValues.clear();
        rShapeFunctionsGradientsValues.reserve(n_total_int_pts);

        // Compute each Gauss pt. shape functions values
        for (unsigned int i_subdivision = 0; i_subdivision < n_subdivision; ++i_subdivision) {

            const IndexedPointGeometryType& r_subdivision_geom = *rSubdivisionsVector[i_subdivision];

            // Get the subdivision shape function values
            const Matrix subdivision_sh_func_values = r_subdivision_geom.ShapeFunctionsValues(IntegrationMethod);
            ShapeFunctionsGradientsType subdivision_sh_func_gradients_values;
            r_subdivision_geom.ShapeFunctionsIntegrationPointsGradients(subdivision_sh_func_gradients_values, IntegrationMethod);

            // Get the subdivision Jacobian values on all Gauss pts.
            Vector subdivision_jacobians_values;
            r_subdivision_geom.DeterminantOfJacobian(subdivision_jacobians_values, IntegrationMethod);

            // Get the subdivision Gauss pts. (x_coord, y_coord, z_coord, weight)
            const IntegrationPointsArrayType subdivision_gauss_points = r_subdivision_geom.IntegrationPoints(IntegrationMethod);

            // Apply the original nodes condensation
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {

                // Store the Gauss pts. weights values
                rWeightsValues(i_subdivision*n_int_pts + i_gauss) = subdivision_jacobians_values(i_gauss) * subdivision_gauss_points[i_gauss].Weight();

                // Condense the shape function local values to obtain the original nodes ones
                Vector sh_func_vect = ZeroVector(split_edges_size);
                for (unsigned int i = 0; i < n_nodes; ++i) {
                    sh_func_vect(r_subdivision_geom[i].Id()) = subdivision_sh_func_values(i_gauss, i);
                }

                const Vector condensed_sh_func_values = prod(sh_func_vect, rPmatrix);

                for (unsigned int i = 0; i < n_nodes; ++i) {
                    rShapeFunctionsValues(i_subdivision*n_int_pts + i_gauss, i) = condensed_sh_func_values(i);
                }

                // Condense the shape function gradients local values to obtain the original nodes ones
                Matrix aux_gauss_gradients = ZeroMatrix(n_nodes, n_dim);

                Matrix sh_func_gradients_mat = ZeroMatrix(n_dim, split_edges_size);
                for (unsigned int dim = 0; dim < n_dim; ++dim ) {
                    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                        sh_func_gradients_mat(dim, r_subdivision_geom[i_node].Id()) = subdivision_sh_func_gradients_values[i_gauss](i_node, dim);
                    }
                }

                rShapeFunctionsGradientsValues.push_back(trans(prod(sh_func_gradients_mat, rPmatrix)));
            }
        }
    };

    // Given the interfaces pattern of either the positive or negative interface side, computes the shape function values.
    void ModifiedShapeFunctions::ComputeInterfaceValuesOnOneSide(
        Matrix &rInterfaceShapeFunctionsValues,
        std::vector<Matrix> &rInterfaceShapeFunctionsGradientsValues,
        Vector &rInterfaceWeightsValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const IntegrationMethodType IntegrationMethod) {

        // Set some auxiliar variables
        const unsigned int n_interfaces = rInterfacesVector.size();                                          // Number of positive or negative subdivisions
        const unsigned int n_nodes = mpInputGeometry->PointsNumber();                                        // Split geometry number of nodes
        const unsigned int n_int_pts = (*rInterfacesVector[0]).IntegrationPointsNumber(IntegrationMethod);   // Number of Gauss pts. per interface
        const unsigned int n_total_int_pts = n_interfaces * n_int_pts;                                       // Total Gauss pts.
        GeometryPointerType p_input_geometry = this->GetInputGeometry();                                     // Pointer to the input geometry

        // Resize the shape function values matrix
        if (rInterfaceShapeFunctionsValues.size1() != n_total_int_pts) {
            rInterfaceShapeFunctionsValues.resize(n_total_int_pts, n_nodes, false);
        } else if (rInterfaceShapeFunctionsValues.size2() != n_nodes) {
            rInterfaceShapeFunctionsValues.resize(n_total_int_pts, n_nodes, false);
        }

        // Resize the weights vector
        if (rInterfaceWeightsValues.size() != n_total_int_pts) {
            rInterfaceWeightsValues.resize(n_total_int_pts, false);
        }

        // Clear the gradients vector
        rInterfaceShapeFunctionsGradientsValues.clear();
        rInterfaceShapeFunctionsGradientsValues.reserve(n_total_int_pts);

        // Compute each Gauss pt. shape functions values
        for (unsigned int i_interface = 0; i_interface < n_interfaces; ++i_interface) {

            const IndexedPointGeometryType& r_interface_geom = *rInterfacesVector[i_interface];

            std::vector < CoordinatesArrayType > interface_gauss_pts_gl_coords, interface_gauss_pts_loc_coords;
            interface_gauss_pts_gl_coords.clear();
            interface_gauss_pts_loc_coords.clear();
            interface_gauss_pts_gl_coords.reserve(n_int_pts);
            interface_gauss_pts_loc_coords.reserve(n_int_pts);

            // Get the intersection Gauss points coordinates
            IntegrationPointsArrayType interface_gauss_pts;
            interface_gauss_pts = r_interface_geom.IntegrationPoints(IntegrationMethod);

            // Get the intersection Jacobians values
            Vector intersection_jacobians;
            r_interface_geom.DeterminantOfJacobian(intersection_jacobians, IntegrationMethod);

            // Store the Gauss points weights
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
                rInterfaceWeightsValues(i_gauss) = intersection_jacobians(i_gauss) * interface_gauss_pts[i_gauss].Weight();
            }

            // Compute the global coordinates of the intersection Gauss pts.
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
                CoordinatesArrayType global_coordinates = ZeroVector(3);
                r_interface_geom.GlobalCoordinates(global_coordinates, interface_gauss_pts[i_gauss].Coordinates());
                interface_gauss_pts_gl_coords.push_back(global_coordinates);
            }

            // Compute the local coordinates of the intersection Gauss pts in the complete geometry.
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
                CoordinatesArrayType local_coordinates = ZeroVector(3);
                local_coordinates = p_input_geometry->PointLocalCoordinates(local_coordinates, interface_gauss_pts_gl_coords[i_gauss]);
                interface_gauss_pts_loc_coords.push_back(local_coordinates);
            }

            // Get the original geometry shape function and gradients values over the intersection
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
                Vector aux_sh_func;
                Matrix aux_grad_sh_func;
                aux_sh_func = p_input_geometry->ShapeFunctionsValues(aux_sh_func, interface_gauss_pts_loc_coords[i_gauss]);
                aux_grad_sh_func = p_input_geometry->ShapeFunctionsLocalGradients(aux_grad_sh_func, interface_gauss_pts_loc_coords[i_gauss]);

                // Store the computed shape function values
                for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
                    rInterfaceShapeFunctionsValues(i_gauss, i_node) = aux_sh_func(i_node);
                }

                // Store the computed gradient values
                rInterfaceShapeFunctionsGradientsValues.push_back(aux_grad_sh_func);
            }
        }
    };

    // Given the interfaces pattern of either the positive or negative interface side, computes the outwards unit normal vector
    void ModifiedShapeFunctions::ComputeInterfaceNormalOnOneSide(
        std::vector<Vector> &rInterfaceUnitNormalValues,
        const std::vector<IndexedPointGeometryPointerType> &rInterfacesVector,
        const IntegrationMethodType IntegrationMethod) {

        // Set some auxiliar variables
        const unsigned int n_interfaces = rInterfacesVector.size();                                          // Number of interfaces
        const unsigned int n_int_pts = (*rInterfacesVector[0]).IntegrationPointsNumber(IntegrationMethod);   // Number of Gauss pts. per interface
        const unsigned int n_total_int_pts = n_interfaces * n_int_pts;                                       // Total Gauss pts.

        rInterfaceUnitNormalValues.clear();
        rInterfaceUnitNormalValues.reserve(n_total_int_pts);

        // Compute each Gauss pt. shape functions values
        for (unsigned int i_interface = 0; i_interface < n_interfaces; ++i_interface) {

            IndexedPointGeometryType& r_interface_geom = *rInterfacesVector[i_interface];
            const unsigned int n_int_pts = r_interface_geom.IntegrationPointsNumber(IntegrationMethod);

            // Get the intersection Gauss points coordinates
            IntegrationPointsArrayType interface_gauss_pts;
            interface_gauss_pts = r_interface_geom.IntegrationPoints(IntegrationMethod);

            // Compute the ouwards unit normal vector values
            for (unsigned int i_gauss = 0; i_gauss < n_int_pts; ++i_gauss) {
                array_1d<double,3> aux_unit_normal = r_interface_geom.UnitNormal(interface_gauss_pts[i_gauss].Coordinates());
                rInterfaceUnitNormalValues.push_back(aux_unit_normal);
            }
        }
    };


}; // namespace Kratos
