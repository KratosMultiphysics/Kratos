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
    const std::size_t n_origin = this->GetOriginModelPart().NumberOfNodes();
    const std::size_t n_destination = this->GetDestinationModelPart().NumberOfNodes();

    // Calculate number of polynomial terms
    const unsigned int poly_degree = mLocalMapperSettings["additional_polynomial_degree"].GetInt();
    const std::size_t n_polynomial = (poly_degree + 3) * (poly_degree + 2) * (poly_degree + 1) / 6;
    const std::size_t origin_interpolation_matrix_size = n_origin + n_polynomial;

    // Allocate mapping matrix
    this->mpMappingMatrix = Kratos::make_unique<TMappingMatrixType>(n_destination, n_origin);

    // Create the coordinate matrix for the origin and destination nodes
    DenseMatrixType origin_coords(n_origin, 3);
    std::size_t i = 0;
    for (const auto& node : this->GetOriginModelPart().Nodes()) {
        for (unsigned int dim = 0; dim < 3; ++dim) {
            origin_coords(i, dim) = node.Coordinates()[dim];
        }
        ++i;
    }

    DenseMatrixType destination_coords(n_destination, 3);
    std::size_t j = 0;
    for (const auto& node : this->GetDestinationModelPart().Nodes()) {
        for (unsigned int dim = 0; dim < 3; ++dim) {
            destination_coords(j, dim) = node.Coordinates()[dim];
        }
        ++j;
    }

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
