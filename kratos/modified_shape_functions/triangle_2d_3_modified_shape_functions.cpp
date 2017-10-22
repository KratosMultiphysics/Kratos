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
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{
    
    /// Triangle2D3ModifiedShapeFunctions implementation
    /// Default constructor
    Triangle2D3ModifiedShapeFunctions::Triangle2D3ModifiedShapeFunctions(GeometryType& rInputGeometry, Vector& rNodalDistances) :
        ModifiedShapeFunctions(rInputGeometry, rNodalDistances) {};

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
        const GeometryType geometry = this->GetInputGeometry();
        const Vector nodal_distances = this->GetNodalDistances();
        rOStream << "Triangle2D3N modified shape functions computation class:\n";
        rOStream << "\tGeometry type: " << geometry.Info() << "\n";
        std::stringstream distances_buffer;
        for (unsigned int i = 0; i < nodal_distances.size(); ++i) {
            distances_buffer << std::to_string(nodal_distances(i)) << " ";
        }
        rOStream << "\tDistance values: " << distances_buffer.str();
    };

    // Internally computes the splitting pattern and returns all the shape function values for both sides.
    void Triangle2D3ModifiedShapeFunctions::GetShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                                                Matrix &rNegativeSideShapeFunctionsValues,
                                                                                std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                                                std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                                                Vector &rPositiveSideWeightsValues,
                                                                                Vector &rNegativeSideWeightsValues,
                                                                                const IntegrationMethodType IntegrationMethod){
        // Build the triangle splitting utility
        Vector nodal_distances = this->GetNodalDistances();
        GeometryType input_geometry = this->GetInputGeometry();
        DivideTriangle2D3 triangle_splitter(input_geometry, nodal_distances);

        // Call the divide geometry method
        DivideTriangle2D3::IndexedPointsContainerType aux_points_set;
        std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > positive_subdivisions, negative_subdivisions;
        bool is_divided = triangle_splitter.GenerateDivision(aux_points_set, positive_subdivisions, negative_subdivisions);

        if (is_divided) {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            SetIntersectionPointsCondensationMatrix(p_matrix, 
                                                    triangle_splitter.mEdgeNodeI, 
                                                    triangle_splitter.mEdgeNodeJ, 
                                                    triangle_splitter.mSplitEdges, 
                                                    triangle_splitter.mSplitEdgesSize);

            // Compute the positive side values
            this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues, 
                                        rPositiveSideShapeFunctionsGradientsValues, 
                                        rPositiveSideWeightsValues,
                                        positive_subdivisions,
                                        p_matrix,
                                        IntegrationMethod);

            // Compute the negative side values
            this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues, 
                                        rNegativeSideShapeFunctionsGradientsValues, 
                                        rNegativeSideWeightsValues,
                                        negative_subdivisions,
                                        p_matrix,
                                        IntegrationMethod);
        } else {
            KRATOS_ERROR << "Using the GetShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the positive side.
    void Triangle2D3ModifiedShapeFunctions::GetPositiveSideShapeFunctionsAndGradientsValues(Matrix &rPositiveSideShapeFunctionsValues,
                                                                                            std::vector<Matrix> &rPositiveSideShapeFunctionsGradientsValues,
                                                                                            Vector &rPositiveSideWeightsValues,
                                                                                            const IntegrationMethodType IntegrationMethod) {
        // Build the triangle splitting utility
        Vector nodal_distances = this->GetNodalDistances();
        GeometryType input_geometry = this->GetInputGeometry();
        DivideTriangle2D3 triangle_splitter(input_geometry, nodal_distances);

        // Call the divide geometry method
        DivideTriangle2D3::IndexedPointsContainerType aux_points_set;
        std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > positive_subdivisions, negative_subdivisions;
        bool is_divided = triangle_splitter.GenerateDivision(aux_points_set, positive_subdivisions, negative_subdivisions);

        if (is_divided) {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            SetIntersectionPointsCondensationMatrix(p_matrix, 
                                                    triangle_splitter.mEdgeNodeI, 
                                                    triangle_splitter.mEdgeNodeJ, 
                                                    triangle_splitter.mSplitEdges, 
                                                    triangle_splitter.mSplitEdgesSize);

            // Compute the positive side values
            this->ComputeValuesOnOneSide(rPositiveSideShapeFunctionsValues, 
                                        rPositiveSideShapeFunctionsGradientsValues, 
                                        rPositiveSideWeightsValues,
                                        positive_subdivisions,
                                        p_matrix,
                                        IntegrationMethod);
        } else {
            KRATOS_ERROR << "Using the GetPositiveSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Internally computes the splitting pattern and returns all the shape function values for the negative side.
    void Triangle2D3ModifiedShapeFunctions::GetNegativeSideShapeFunctionsAndGradientsValues(Matrix &rNegativeSideShapeFunctionsValues,
                                                                                            std::vector<Matrix> &rNegativeSideShapeFunctionsGradientsValues,
                                                                                            Vector &rNegativeSideWeightsValues,
                                                                                            const IntegrationMethodType IntegrationMethod) {
        // Build the triangle splitting utility
        Vector nodal_distances = this->GetNodalDistances();
        GeometryType input_geometry = this->GetInputGeometry();
        DivideTriangle2D3 triangle_splitter(input_geometry, nodal_distances);

        // Call the divide geometry method
        DivideTriangle2D3::IndexedPointsContainerType aux_points_set;
        std::vector < DivideTriangle2D3::IndexedPointGeometryPointerType > positive_subdivisions, negative_subdivisions;
        bool is_divided = triangle_splitter.GenerateDivision(aux_points_set, positive_subdivisions, negative_subdivisions);

        if (is_divided) {
            // Get the intersection points condensation matrix
            Matrix p_matrix;
            SetIntersectionPointsCondensationMatrix(p_matrix, 
                                                    triangle_splitter.mEdgeNodeI, 
                                                    triangle_splitter.mEdgeNodeJ, 
                                                    triangle_splitter.mSplitEdges, 
                                                    triangle_splitter.mSplitEdgesSize);

            // Compute the negative side values
            this->ComputeValuesOnOneSide(rNegativeSideShapeFunctionsValues, 
                                        rNegativeSideShapeFunctionsGradientsValues, 
                                        rNegativeSideWeightsValues,
                                        negative_subdivisions,
                                        p_matrix,
                                        IntegrationMethod);
        } else {
            KRATOS_ERROR << "Using the GetNegativeSideShapeFunctionsAndGradientsValues method for a non divided geometry.";
        }
    };

    // Given the subdivision pattern of either the positive or negative side, computes the shape function values.
    void Triangle2D3ModifiedShapeFunctions::ComputeValuesOnOneSide(Matrix &rShapeFunctionsValues,
                                                                   std::vector<Matrix> &rShapeFunctionsGradientsValues,
                                                                   Vector &rWeightsValues,
                                                                   const std::vector<IndexedPointGeometryPointerType> &rSubdivisionsVector,
                                                                   const Matrix &p_matrix,
                                                                   const IntegrationMethodType IntegrationMethod)
    {
        // Compute some auxiliar constants
        GeometryType input_geometry = this->GetInputGeometry();
        const unsigned int n_edges_global = input_geometry.EdgesNumber();      // Number of edges in original geometry
        const unsigned int n_nodes_global = input_geometry.PointsNumber();    // Number of nodes in original geometry
        const unsigned int split_edges_size = n_nodes_global + n_edges_global; // Split edges vector size 

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
                
                const Vector condensed_sh_func_values = prod(sh_func_vect, p_matrix);
                
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

                rShapeFunctionsGradientsValues.push_back(trans(prod(sh_func_gradients_mat, p_matrix)));
            }
        }
    };
};
