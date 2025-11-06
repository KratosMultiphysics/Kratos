//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//

// System includes

// External includes

// Project includes
#include "radial_basis_function_mapper.h"
#include "mapping_application_variables.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapper_utilities.h"

// External includes
#include <unordered_set>

namespace Kratos {

/* Naming convention for the different matrices involved:
    origin_interpolation_matrix = [ 0  Ps^T ; Ps  Ks ] ---> [beta gamma]^T = inverse_origin_interpolation_matrix * [0 us]^T
    origin_interpolation_vector = [ us ; 0 ]
    destination_interpolation_matrix = [ Pf  Kf ]
*/

template<class TSparseSpace, class TDenseSpace>
RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::RadialBasisFunctionMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters)
    : BaseType(rModelPartOrigin, rModelPartDestination, JsonParameters),
    mLocalMapperSettings(JsonParameters) 
{
    mLocalMapperSettings.ValidateAndAssignDefaults(this->GetMapperDefaultSettings());

    //const bool is_destination_slave = mLocalMapperSettings["is_destination_slave"].GetBool();

    // Define the master and slave interfaces
    mpCouplingInterfaceMaster = &rModelPartOrigin;
    mpCouplingInterfaceSlave = &rModelPartDestination;

    MappingMatrixUtilitiesType::InitializeSystemVector(
        this->pGetInterfaceVectorContainerOrigin()->pGetVector(),
        mpCouplingInterfaceMaster->NumberOfNodes());

    MappingMatrixUtilitiesType::InitializeSystemVector(
        this->pGetInterfaceVectorContainerDestination()->pGetVector(),
        mpCouplingInterfaceSlave->NumberOfNodes());
    
    // This function creates the mapping matrix required for the linear mapping
    this->InitializeInterface();
}


template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // Get a vector with the integration points in the origin and destination coordinates
    std::vector<Condition::Pointer> origin_integration_points;
    std::vector<Condition::Pointer> destination_integration_points;

     // Determine whether we project or not the origin nodes coordinates to the plane of the panel mesh 
    const bool project_origin_nodes_to_destination_domain = mLocalMapperSettings["destination_solver_settings"]["project_origin_nodes_to_destination_domain"].GetBool();
    // Determine whether we map structural displacements to panels angle of attack or not
    const bool map_displacements_to_rotations = mLocalMapperSettings["destination_solver_settings"]["map_displacements_to_rotations"].GetBool();

    // Obtain the polynomial degree and calculate number of polynomial terms
    const unsigned int poly_degree = mLocalMapperSettings["additional_polynomial_degree"].GetInt();
    IndexType n_polynomial = CalculateNumberOfPolynomialTermsFromDegree(poly_degree, project_origin_nodes_to_destination_domain);

    // Determine whether the origin and destination domain are IGA discretizations or not
    const bool is_origin_iga = mLocalMapperSettings["is_origin_iga"].GetBool();
    const bool is_destination_iga = mLocalMapperSettings["is_destination_iga"].GetBool();

    // Throw an error iF the destination domain is discretized with IBRA
    KRATOS_ERROR_IF(is_destination_iga)<< "This mapper is not available when the destination domain is discretized with IBRA" << std::endl;

    // Remember that in the IgaApplication, conditions are represented geometrically by integrtion points
    auto collect_conditions = [](const ModelPart& model_part, auto& integration_points_vector) {
        for(auto condition_it = model_part.ConditionsBegin(); condition_it != model_part.ConditionsEnd(); ++condition_it){
            integration_points_vector.push_back(*condition_it.base());
        }
    };

    if (is_origin_iga) {
        collect_conditions(this->GetOriginModelPart(), origin_integration_points);
    }

    if (is_destination_iga) {
        collect_conditions(this->GetDestinationModelPart(), destination_integration_points);
    }
    
    // Get the sizes of the origin and interpolation interfaces for the mapping
    IndexType n_origin = is_origin_iga
        ? origin_integration_points.size()
        : this->GetOriginModelPart().NumberOfNodes();
    IndexType n_destination = is_destination_iga
        ? destination_integration_points.size()
        : this->GetDestinationModelPart().NumberOfNodes();


    // Create and fill the matrix with the coordinates of the origin interpolation points
    DenseMatrixType origin_coords(n_origin, 3);
    FillCoordinatesMatrix(this->GetOriginModelPart(), origin_integration_points, origin_coords, is_origin_iga);

    // Create and fill the matrix with the coordinates of the destination interpolation points
    DenseMatrixType destination_coords(n_destination, 3);
    FillCoordinatesMatrix(this->GetDestinationModelPart(), destination_integration_points, destination_coords, is_destination_iga);

    // Allocate mapping matrix
    this->mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(n_destination, n_origin);

    // Read string from settings
    std::string rbf_string = mLocalMapperSettings["radial_basis_function_type"].GetString();
    std::transform(rbf_string.begin(), rbf_string.end(), rbf_string.begin(), ::tolower);

    // Convert to enum using a local static map
    static const std::unordered_map<std::string, RBFShapeFunctionsUtility::RBFType> rbf_type_map = {
        {"inverse_multiquadric", RBFShapeFunctionsUtility::RBFType::InverseMultiquadric},
        {"gaussian",             RBFShapeFunctionsUtility::RBFType::Gaussian},
        {"thin_plate_spline",    RBFShapeFunctionsUtility::RBFType::ThinPlateSpline},
        {"wendland_c2",          RBFShapeFunctionsUtility::RBFType::WendlandC2}
    };

    auto it = rbf_type_map.find(rbf_string);
    if (it == rbf_type_map.end()) {
        KRATOS_ERROR << "Unrecognized RBF type: " << rbf_string << std::endl;
    }

    // Now rbf_type is an enum  
    RBFShapeFunctionsUtility::RBFType rbf_type = it->second;

    // Compute shape parameter
    double rbf_shape_parameter = 0.1; // default shape parameter
    if (rbf_type == RBFShapeFunctionsUtility::RBFType::InverseMultiquadric || rbf_type == RBFShapeFunctionsUtility::RBFType::Gaussian) 
    {
        rbf_shape_parameter = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(origin_coords);
    }
    else if (rbf_type == RBFShapeFunctionsUtility::RBFType::WendlandC2) 
    {
        rbf_shape_parameter = RBFShapeFunctionsUtility::CalculateWendlandC2SupportRadius(origin_coords, 2.5);
    }

    double rbf_scale_factor = CalculateScaleFactor(origin_coords);

    // Build origin interpolation matrix and invert it
    const IndexType origin_interpolation_matrix_size = n_origin + n_polynomial;
    DenseMatrixType inverse_origin_interpolation_matrix(origin_interpolation_matrix_size, origin_interpolation_matrix_size);
    CreateAndInvertOriginRBFMatrix(inverse_origin_interpolation_matrix, origin_coords, project_origin_nodes_to_destination_domain,
        poly_degree, rbf_type, rbf_scale_factor, rbf_shape_parameter);

    // Build destination interpolation matrix (considering distances between origin and destination interpolation points)
    DenseMatrixType destination_interpolation_matrix(n_destination, n_origin + n_polynomial);
    CreateDestinationRBFMatrix(destination_interpolation_matrix, origin_coords, destination_coords, project_origin_nodes_to_destination_domain, poly_degree,
        rbf_type, map_displacements_to_rotations, rbf_shape_parameter);
    
    // Compute the mapping matrix
    DenseMatrixType dense_mapping = prod(destination_interpolation_matrix, inverse_origin_interpolation_matrix);

    // Adjusting the size of mpMappingMatrix
    const IndexType n_rows = dense_mapping.size1();
    const IndexType n_cols = dense_mapping.size2() - n_polynomial; // eliminamos columnas polinomiales
    this->mpMappingMatrix->resize(n_rows, n_cols, false);

    // we copy only the required columns for the mapping
    for (IndexType i = 0; i < n_rows; ++i) {
        for (IndexType j = 0; j < n_cols; ++j) {
            if (std::abs(dense_mapping(i, j + n_polynomial)) > 1e-12) { // shift por columnas polinomiales
                (*this->mpMappingMatrix)(i, j) = dense_mapping(i, j + n_polynomial);
            }
        }
    }

    if (is_origin_iga) {
        this->mpMappingMatrix = ComputeMappingMatrixIga(*this->mpMappingMatrix, origin_integration_points, this->GetOriginModelPart());
    }
}

template<class TSparseSpace, class TDenseSpace>
void Kratos::RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::FillCoordinatesMatrix(
    const ModelPart& ModelPart, 
    const std::vector<Condition::Pointer>& Integration_pointsPointerVector, 
    DenseMatrixType& rCoordinatesMatrix, 
    bool IsDomainIGA)
{
    if(IsDomainIGA){
        // If the origin domain is discretized with a NURBS mesh, we consider the integration points (conditions) as the interpolation points
         for (IndexType i = 0; i < Integration_pointsPointerVector.size(); ++i) {
            const auto& p_geometry = Integration_pointsPointerVector[i]->pGetGeometry();
            const auto& integration_point_coords = p_geometry->Center(); 

            rCoordinatesMatrix(i, 0) = integration_point_coords[0];
            rCoordinatesMatrix(i, 1) = integration_point_coords[1];
            rCoordinatesMatrix(i, 2) = integration_point_coords[2];
        }
    } else{
        // If the origin domain is discretized with a FEM mesh, we consider the nodes as the interpolation points
        IndexType i = 0;
        for (const auto& node : ModelPart.Nodes()) {
            const auto& node_coords = node.Coordinates();

            rCoordinatesMatrix(i, 0) = node_coords[0];
            rCoordinatesMatrix(i, 1) = node_coords[1];
            rCoordinatesMatrix(i, 2) = node_coords[2];

            i++;
        }
    } 
}
template<class TSparseSpace, class TDenseSpace>
IndexType RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CalculateNumberOfPolynomialTermsFromDegree(
    IndexType PolyDegree, 
    bool ProjectToAerodynamicPanels)
{
    if (!ProjectToAerodynamicPanels){
        return (PolyDegree + 3) * (PolyDegree + 2) * (PolyDegree + 1) / 6;
    } else{
        return (PolyDegree + 2) * (PolyDegree + 1) / 2;
    }
}

template<class TSparseSpace, class TDenseSpace>
std::vector<double> RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::EvaluatePolynomialBasis(
    const array_1d<double, 3>& rCoords, 
    unsigned int degree, 
    bool ProjectToPanelsPlane) const
{
    std::vector<double> values;
    const double x = rCoords[0];
    const double y = rCoords[1];
    const double z = rCoords[2];

    if (ProjectToPanelsPlane) {
        // 2D polynomial basis in standard order: 1, x, y, x^2, xy, y^2, ...
        for (unsigned int j = 0; j <= degree; ++j) {           // power of y
            for (unsigned int i = 0; i <= degree - j; ++i) {   // power of x
                values.push_back(std::pow(x, i) * std::pow(y, j));
            }
        }
    } else {
        // 3D polynomial basis
        for (unsigned int i = 0; i <= degree; ++i) {
            for (unsigned int j = 0; j <= degree - i; ++j) {
                for (unsigned int k = 0; k <= degree - i - j; ++k) {
                    values.push_back(std::pow(x, i) * std::pow(y, j) * std::pow(z, k));
                }
            }
        }
    }

    return values;
}

template<class TSparseSpace, class TDenseSpace>
double RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CalculateScaleFactor(
    DenseMatrixType& rOriginCoords)
{
    // rOriginCoords is (n_points x 3), columns: 0=x, 1=y, 2=z (ignored)
    const IndexType n_origin = rOriginCoords.size1(); // number of points
    if (n_origin == 0) {
        // handle error or return 0
        return 0.0;
    }

    // Initialize min/max with the first point
    double minx = rOriginCoords(0, 0);
    double maxx = minx;
    double miny = rOriginCoords(0, 1);
    double maxy = miny;

    for (IndexType i = 0; i < n_origin; ++i) {
        double x = rOriginCoords(i, 0);
        double y = rOriginCoords(i, 1);

        if (x > maxx) maxx = x;
        if (x < minx) minx = x;
        if (y > maxy) maxy = y;
        if (y < miny) miny = y;
    }

    // Compute ranges
    double rangeX = maxx - minx;
    double rangeY = maxy - miny;

    // Return the larger of the two ranges
    return (rangeX >= rangeY) ? rangeX : rangeY;
}

template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CreateAndInvertOriginRBFMatrix(
    DenseMatrixType& rInvCMatrix, // OUT: inverted RBF + poly matrix
    const DenseMatrixType& rOriginCoords,
    bool ProjectToAerodynamicPanels,
    IndexType Poly_Degree,
    RBFShapeFunctionsUtility::RBFType RBF_Type,
    double Factor,
    double eps)
{
    // Determine number of points in the origin domain
    const IndexType n_points = rOriginCoords.size1(); // n_points x 3

    IndexType n_polynomial = CalculateNumberOfPolynomialTermsFromDegree(Poly_Degree, ProjectToAerodynamicPanels);
    IndexType n_total = n_points + n_polynomial;

    // Initialize the matrix
    noalias(rInvCMatrix) = ZeroMatrix(n_total, n_total);

    for (IndexType i = 0; i < n_points; ++i) {
        array_1d<double,3> point;
        point[0] = rOriginCoords(i,0);
        point[1] = rOriginCoords(i,1);
        point[2] = rOriginCoords(i,2);

        const auto poly_vals = EvaluatePolynomialBasis(point, Poly_Degree, ProjectToAerodynamicPanels);

        for (IndexType j = 0; j < n_polynomial; ++j) {
            // P block: rows n_poly..n_poly+n_points-1, columns 0..n_poly-1
            rInvCMatrix(i + n_polynomial, j) = poly_vals[j];

            // P^T block: top-left n_poly x n_points
            rInvCMatrix(j, i + n_polynomial) = poly_vals[j]; // transpose
        }
    }
    for (IndexType i = 0; i < n_points; ++i) {
        for (IndexType j = i; j < n_points; ++j) {
            double dx = rOriginCoords(i,0) - rOriginCoords(j,0);
            double dy = rOriginCoords(i,1) - rOriginCoords(j,1);
            double dz = ProjectToAerodynamicPanels ? 0.0 : rOriginCoords(i,2) - rOriginCoords(j,2);
            double r = std::sqrt(dx*dx + dy*dy + dz*dz);

            double k = RBFShapeFunctionsUtility::EvaluateRBF(r, eps, RBF_Type); 

            rInvCMatrix(i + n_polynomial, j + n_polynomial) = k;
            rInvCMatrix(j + n_polynomial, i + n_polynomial) = k; // symmetry
        }
    }

    // Scale the coefficient matrix to improve its conditioning using the semispan or chord as a scaling factor
   double scale_factor = 1.0/Factor;
   double scale_factor_2 = std::pow(scale_factor, 2);

    // Scale the polynomial block (first n_polynomial columns)
    for (IndexType i = 0; i < n_polynomial; ++i) {
        for (IndexType j = 0; j < n_points; ++j) {
            // Rows corresponding to points (offset by n_polynomial), columns = polynomial
            rInvCMatrix(j + n_polynomial, i) *= scale_factor;

            // Symmetric block: top-left
            rInvCMatrix(i, j + n_polynomial) = rInvCMatrix(j + n_polynomial, i);
        }
    }

    // Scale the RBF submatrix (distance part)
    for (IndexType i = 0; i < n_points; ++i) {
        for (IndexType j = i; j < n_points; ++j) {
            rInvCMatrix(i + n_polynomial, j + n_polynomial) *= scale_factor_2;
            rInvCMatrix(j + n_polynomial, i + n_polynomial) = rInvCMatrix(i + n_polynomial, j + n_polynomial); // symmetry
        }
    }

    DenseMatrixType C = rInvCMatrix; // your coefficient matrix, size (n_points+n_polynomial) x (n_points+n_polynomial)
    double det_c_inv;

    // Assuming C is your filled coefficient matrix
    MathUtils<double>::InvertMatrix(C, rInvCMatrix, det_c_inv);

    // First n_polynomial columns/rows (polynomial part)
    for (IndexType i = 0; i < n_polynomial; ++i) {
        for (IndexType j = n_polynomial; j < n_total; ++j) {
            rInvCMatrix(j, i) *= scale_factor;
            rInvCMatrix(i, j) = rInvCMatrix(j, i); // symmetric
        }
    }

    // RBF submatrix block
    for (IndexType i = n_polynomial; i < n_total; ++i) {
        for (IndexType j = i; j < n_total; ++j) {
            rInvCMatrix(j, i) *= scale_factor_2;
            rInvCMatrix(i, j) = rInvCMatrix(j, i); // symmetric
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CreateDestinationRBFMatrix(
    DenseMatrixType& rAMatrix,
    const DenseMatrixType& rOriginCoords, 
    const DenseMatrixType& rDestinationCoords,
    bool ProjectToAerodynamicPanels, 
    IndexType Poly_Degree, 
    RBFShapeFunctionsUtility::RBFType RBF_Type, 
    bool map_displacements_to_rotations,
    double rbf_shape_parameter)
{
    const IndexType n_origin = rOriginCoords.size1();
    const IndexType n_destination = rDestinationCoords.size1();

    IndexType n_polynomial = CalculateNumberOfPolynomialTermsFromDegree(Poly_Degree, ProjectToAerodynamicPanels);

    // === Main assembly ===
    for (IndexType i_dest = 0; i_dest < n_destination; ++i_dest) {
        array_1d<double, 3> dest_point;
        dest_point[0] = rDestinationCoords(i_dest, 0);
        dest_point[1] = rDestinationCoords(i_dest, 1);
        dest_point[2] = rDestinationCoords(i_dest, 2);

        // === Polynomial block ===
        std::vector<double> poly_vals(n_polynomial, 0.0);
        if (!map_displacements_to_rotations) {
            poly_vals = EvaluatePolynomialBasis(dest_point, Poly_Degree, ProjectToAerodynamicPanels);
        } else {
            // Slope/AOA mode: constant derivative in x as per Fortran [0, -1, 0]
            // (extendable if higher polynomial degrees are used)
            if (n_polynomial >= 3) {
                poly_vals[0] = 0.0;
                poly_vals[1] = -1.0;
                poly_vals[2] = 0.0;
            } else {
                poly_vals[0] = 0.0;
            }
        }

        for (IndexType j = 0; j < n_polynomial; ++j) {
            rAMatrix(i_dest, j) = poly_vals[j];
        }

        // === RBF block ===
        for (IndexType j_origin = 0; j_origin < n_origin; ++j_origin) {
            array_1d<double, 3> origin_point;
            origin_point[0] = rOriginCoords(j_origin, 0);
            origin_point[1] = rOriginCoords(j_origin, 1);
            origin_point[2] = rOriginCoords(j_origin, 2);

            double dx = dest_point[0] - origin_point[0];
            double dy = dest_point[1] - origin_point[1];
            double dz = ProjectToAerodynamicPanels ? 0.0 : (dest_point[2] - origin_point[2]);

            const double distance_2 = dx * dx + dy * dy + dz * dz;

            double value = 0.0;
            if (!map_displacements_to_rotations) {
                const double r = std::sqrt(distance_2);
                value = RBFShapeFunctionsUtility::EvaluateRBF(r, rbf_shape_parameter, RBF_Type);
            } else {
                // slope / angle of attack version: -dK/dx = -2*dx*(1 + log(r^2))
                const double dkdx = 2.0 * dx * (1.0 + std::log(distance_2));
                value = -dkdx;
            }

            rAMatrix(i_dest, n_polynomial + j_origin) = value;
        }
    }
}


template<class TSparseSpace, class TDenseSpace>
std::unique_ptr<typename RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::TMappingMatrixType>
RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::ComputeMappingMatrixIga(
    const typename RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::TMappingMatrixType& rMappingMatrixGP,
    const std::vector<Condition::Pointer>& rOriginIntegrationPoints,
    const ModelPart& rOriginModelPart) const
{
    // Determine sizes
    const IndexType n_gp = rOriginIntegrationPoints.size();
    const IndexType n_nodes = rOriginModelPart.GetRootModelPart().NumberOfNodes();

    // Build full N matrix (GP x nodes)
    SparseMatrixType N = ZeroMatrix(n_gp, n_nodes);
    for (IndexType i_gp = 0; i_gp < n_gp; ++i_gp) {
        auto p_geometry = rOriginIntegrationPoints[i_gp]->pGetGeometry();
        const auto& shape_functions_values = p_geometry->ShapeFunctionsValues();

        for (IndexType j = 0; j < p_geometry->size(); ++j) {
            IndexType node_id = (*p_geometry)[j].Id() - 1; // Convert 1-based ID to 0-based
            N(i_gp, node_id) = shape_functions_values(0, j);
        }
    }

    // Determine valid columns corresponding to nodes in the origin model part
    std::unordered_set<IndexType> existing_node_ids;
    for (const auto& node : rOriginModelPart.Nodes()) {
        existing_node_ids.insert(node.Id());
    }

    std::vector<IndexType> valid_columns;
    valid_columns.reserve(n_nodes);
    for (IndexType j = 0; j < n_nodes; ++j) {
        if (existing_node_ids.count(j + 1)) {
            valid_columns.push_back(j);
        }
    }

    // Build reduced N matrix with just the control points belonging to the coupling interface 
    SparseMatrixType N_reduced(n_gp, valid_columns.size());
    for (IndexType i = 0; i < n_gp; ++i) {
        for (IndexType j = 0; j < valid_columns.size(); ++j) {
            N_reduced(i, j) = N(i, valid_columns[j]);
        }
    }


    // Compute final mapping: H = H_tilde @ N_reduced
    auto pMappingMatrix = Kratos::make_unique<TMappingMatrixType>(rMappingMatrixGP.size1(), N_reduced.size2());
    noalias(*pMappingMatrix) = prod(rMappingMatrixGP, N_reduced);

    return pMappingMatrix;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class RadialBasisFunctionMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
