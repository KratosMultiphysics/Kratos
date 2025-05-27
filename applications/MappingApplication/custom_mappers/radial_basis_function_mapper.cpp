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
#include "utilities/variable_utils.h"
#include "factories/linear_solver_factory.h"

// External includes
#include <unordered_set>

namespace Kratos {

/* Naming convention for the different matrices involved:
    OriginInterpolationMatrix = [ Mss  Ps ; Ps^T  0 ] ---> MappingOriginMatrix [gamma beta]^T = [us 0]^T
    OriginInterpolationVector = [ us ; 0 ]
    DestinationEvaluationMatrix = [ Mfs  Pf ]
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

    const bool destination_is_slave = mLocalMapperSettings["destination_is_slave"].GetBool();

    mpCouplingInterfaceMaster = destination_is_slave
        ? &rModelPartOrigin
        : &rModelPartDestination;

    mpCouplingInterfaceSlave = destination_is_slave
        ? &rModelPartDestination
        : &rModelPartOrigin;

    MappingMatrixUtilitiesType::InitializeSystemVector(
        this->pGetInterfaceVectorContainerOrigin()->pGetVector(),
        mpCouplingInterfaceMaster->NumberOfNodes());

    MappingMatrixUtilitiesType::InitializeSystemVector(
        this->pGetInterfaceVectorContainerDestination()->pGetVector(),
        mpCouplingInterfaceSlave->NumberOfNodes());

    this->CreateLinearSolver();
    this->InitializeInterface();
}


template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // Get number of nodes
    std::size_t n_origin = this->GetOriginModelPart().NumberOfNodes();
    std::size_t n_destination = this->GetDestinationModelPart().NumberOfNodes();

    // Get a vector with the integration points in the origin and destination coordinates
    std::vector<Condition::Pointer> origin_integration_points;
    std::vector<Condition::Pointer> destination_integration_points;

    // Calculate number of polynomial terms
    const unsigned int poly_degree = mLocalMapperSettings["additional_polynomial_degree"].GetInt();
    const std::size_t n_polynomial = (poly_degree + 3) * (poly_degree + 2) * (poly_degree + 1) / 6;

    // Determine whether the origin and destination domain are IGA discretizations or not
    const bool is_origin_iga = mLocalMapperSettings["is_origin_iga"].GetBool();
    const bool is_destination_iga = mLocalMapperSettings["is_destination_iga"].GetBool();

    // If the origin is IGA, determine the number of integration points over the coupling boundary
    if (is_origin_iga == 1){
        for(auto condition_it = this->GetOriginModelPart().ConditionsBegin(); condition_it != this->GetOriginModelPart().ConditionsEnd(); ++condition_it){
            std::stringstream condition_info;
            condition_it->PrintInfo(condition_info);
            std::string condition_name = condition_info.str();

             // Skip if it's a Point load Condition
            if (condition_name.find("Point load Condition") != std::string::npos) {
                continue;  // Skip this condition
            }

            origin_integration_points.push_back(*condition_it.base());
        }
    }

    // If the destination is IGA, determine the number of integration points over the coupling boundary
    if (is_destination_iga == 1){
        for(auto condition_it = this->GetDestinationModelPart().ConditionsBegin(); condition_it != this->GetDestinationModelPart().ConditionsEnd(); ++condition_it){
            std::stringstream condition_info;
            condition_it->PrintInfo(condition_info);
            std::string condition_name = condition_info.str();

             // Skip if it's a Point load Condition
            if (condition_name.find("Point load Condition") != std::string::npos) {
                continue;  // Skip this condition
            }

            destination_integration_points.push_back(*condition_it.base());
        }
    }
    

    // Create the coordinate matrix for the origin and destination nodes
    if (is_origin_iga == 1) n_origin = origin_integration_points.size();
    const std::size_t origin_interpolation_matrix_size = n_origin + n_polynomial;
    DenseMatrixType origin_coords(n_origin, 3);

    if(is_origin_iga == 1){
        // If the origin domain is discretized with a NURBS mesh, we consider the integration points (conditions) as the interpolation points
         for (std::size_t i = 0; i < origin_integration_points.size(); ++i) {
            auto p_geometry = origin_integration_points[i]->pGetGeometry();
            auto center = p_geometry->Center(); 

            for (unsigned int dim = 0; dim < 3; ++dim) {
                origin_coords(i, dim) = center[dim];
            }
        }
    } else if(is_origin_iga == 0){
        // If the origin domain is discretized with a FEM mesh, we consider the nodes as the interpolation points
        std::size_t i = 0;
        for (const auto& node : this->GetOriginModelPart().Nodes()) {
            for (unsigned int dim = 0; dim < 3; ++dim) {
                origin_coords(i, dim) = node.Coordinates()[dim];
            }
            ++i;
        }
    } 

    if (is_destination_iga == 1) n_destination = destination_integration_points.size();
    DenseMatrixType destination_coords(n_destination, 3);
    if(is_destination_iga == 1){
        // If the destination domain is discretized with a NURBS mesh, we consider the integration points (conditions) as the interpolation points
         for (std::size_t i = 0; i < destination_integration_points.size(); ++i) {
            auto p_geometry = destination_integration_points[i]->pGetGeometry();
            auto center = p_geometry->Center(); 

            for (unsigned int dim = 0; dim < 3; ++dim) {
                destination_coords(i, dim) = center[dim];
            }
        }
    } else if(is_destination_iga == 0){
        // If the destination domain is discretized with a FEM mesh, we consider the nodes as the interpolation points
        std::size_t i = 0;
        for (const auto& node : this->GetDestinationModelPart().Nodes()) {
            for (unsigned int dim = 0; dim < 3; ++dim) {
                destination_coords(i, dim) = node.Coordinates()[dim];
            }
            ++i;
        }
    } 

    // Allocate mapping matrix
    this->mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(n_destination, n_origin);

    // Compute shape parameter
    const std::string rbf_type = mLocalMapperSettings["radial_basis_function_type"].GetString();
    double h = 0.1;
    if (rbf_type == "inverse_multiquadric" || rbf_type == "gaussian") {
        h = RBFShapeFunctionsUtility::CalculateInverseMultiquadricShapeParameter(origin_coords);
    } else if (rbf_type == "wendland_c2") {
        h = RBFShapeFunctionsUtility::CalculateWendlandC2SupportRadius(origin_coords, 2.5);
    }

    // Build interpolation matrix
    SparseMatrixType origin_interpolation_matrix(origin_interpolation_matrix_size, origin_interpolation_matrix_size);
    for (IndexType r = 0; r < n_origin; ++r) {
        for (IndexType c = 0; c < n_origin; ++c) {
            const double r_ij = norm_2(row(origin_coords, r) - row(origin_coords, c));
            origin_interpolation_matrix(r, c) = RBFShapeFunctionsUtility::EvaluateRBF(r_ij, h, rbf_type);
        }

        const auto poly_vals = EvaluatePolynomialBasis(row(origin_coords, r), poly_degree);
        for (IndexType j = 0; j < n_polynomial; ++j) {
            origin_interpolation_matrix(r, n_origin + j) = poly_vals[j];
            origin_interpolation_matrix(n_origin + j, r) = poly_vals[j];
        }
    }

    // Solve for each destination node
    Vector b(origin_interpolation_matrix_size, 0.0);
    Vector solution(origin_interpolation_matrix_size, 0.0);

    IndexType dest_idx = 0;
    for (IndexType dest_idx = 0; dest_idx < n_destination; ++dest_idx) {
        // Load destination point from matrix
        array_1d<double, 3> destination_point;
        for (unsigned int dim = 0; dim < 3; ++dim) {
            destination_point[dim] = destination_coords(dest_idx, dim);
        }

        // Fill RHS vector b
        for (IndexType j = 0; j < n_origin; ++j) {
            const double r_j = norm_2(destination_point - row(origin_coords, j));
            b[j] = RBFShapeFunctionsUtility::EvaluateRBF(r_j, h, rbf_type);
        }

        const auto destination_poly_vals = EvaluatePolynomialBasis(destination_point, poly_degree);
        for (IndexType j = 0; j < n_polynomial; ++j) {
            b[n_origin + j] = destination_poly_vals[j];
        }

        // Solve linear system
        mpLinearSolver->Solve(origin_interpolation_matrix, solution, b);

        // Fill mapping matrix
        for (IndexType j = 0; j < n_origin; ++j) {
            (*this->mpMappingMatrix)(dest_idx, j) = solution[j];
        }
    }

    // If the origin is iga, we need to transform the gauss point values to control point values ---> H = H_tilde(from gp to fem) @ N
    // N is the matrix of shape functions    
    if(is_origin_iga == 1){
        DenseMatrixType N;
        N = ZeroMatrix(origin_integration_points.size(), this->GetOriginModelPart().GetRootModelPart().NumberOfNodes());
        for (std::size_t i = 0; i < origin_integration_points.size(); ++i) {
            auto p_geometry = origin_integration_points[i]->pGetGeometry();
            const auto& shape_functions = p_geometry->ShapeFunctionsValues();

            for (std::size_t j = 0; j < p_geometry->size(); ++j) {
                std::size_t node_id = (*p_geometry)[j].Id() - 1;  // IDs inician en 1
                N(i, node_id) = shape_functions(0, j);
            }
        }
        // Now we filter the matrix N keeping just those rows and columns which belong to the origin model part 
        std::unordered_set<std::size_t> existing_node_ids;
        for (const auto& node : this->GetOriginModelPart().Nodes()) {
            existing_node_ids.insert(node.Id());
        }

        std::vector<std::size_t> valid_columns;

        for (std::size_t j = 0; j < N.size2(); ++j) {
            std::size_t node_id = j + 1;  
            if (existing_node_ids.count(node_id)) {
                valid_columns.push_back(j);
            }
        }

        DenseMatrixType N_reduced(N.size1(), valid_columns.size());

        for (std::size_t i = 0; i < N.size1(); ++i) {
            for (std::size_t j = 0; j < valid_columns.size(); ++j) {
                N_reduced(i, j) = N(i, valid_columns[j]);
            }
        }

        // Compute the mapping matrix tilde as H_tilde = H @ N
        DenseMatrixType mapping_matrix_tilde;
        mapping_matrix_tilde = ZeroMatrix(this->GetDestinationModelPart().NumberOfNodes(), this->GetOriginModelPart().NumberOfNodes());
        noalias(mapping_matrix_tilde) = prod(*this->mpMappingMatrix, N_reduced);
        
        // Resize the mpMappingMatrix and fill it with the values of mapping_matrix_tilde
        this->mpMappingMatrix->resize(this->GetDestinationModelPart().NumberOfNodes(), this->GetOriginModelPart().NumberOfNodes(), false); 
        noalias(*this->mpMappingMatrix) = mapping_matrix_tilde;
    }
}

template<class TSparseSpace, class TDenseSpace>
std::vector<double> RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::EvaluatePolynomialBasis(
    const array_1d<double, 3>& coords, unsigned int degree) const
{
    std::vector<double> values;
    const double x = coords[0];
    const double y = coords[1];
    const double z = coords[2];

    for (unsigned int i = 0; i <= degree; ++i) {
        for (unsigned int j = 0; j <= degree - i; ++j) {
            for (unsigned int k = 0; k <= degree - i - j; ++k) {
                values.push_back(std::pow(x, i) * std::pow(y, j) * std::pow(z, k));
            }
        }
    }
    return values;
}

template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    if (mLocalMapperSettings["linear_solver_settings"].Has("solver_type")) {
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mLocalMapperSettings["linear_solver_settings"]);
    }
    else {
        Parameters default_settings = this->GetMapperDefaultSettings();
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(default_settings["linear_solver_settings"]);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class RadialBasisFunctionMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
