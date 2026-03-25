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

#include <chrono>
#include <iostream>
#include <iomanip>

// External includes

namespace Kratos {

void RadialBasisFunctionsMapperInterfaceInfoIGA::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    KRATOS_TRY

    const double distance = MapperUtilities::ComputeDistance(this->Coordinates(), rInterfaceObject.Coordinates());

    const auto p_geom = rInterfaceObject.pGetBaseGeometry();

    // Add a geometry candidate to the RBF support accumulator
    mRBFSupportAccumulator.AddGeometryCandidate(p_geom->GetValue(INTERFACE_EQUATION_ID), distance, p_geom->Center());

    // Reorder the support points according to the distance (starting from the closest)
    mRBFSupportAccumulator.ReorderSupportPoints();

    // If the accumulator has enough support points, stop searching 
    if (mRBFSupportAccumulator.HasEnough()) {
        SetLocalSearchWasSuccessful();
    } else if (mRBFSupportAccumulator.HasSome()) {
        SetIsApproximation();
    }

    KRATOS_CATCH("")
}

void RadialBasisFunctionsMapperInterfaceInfoFEM::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    KRATOS_TRY

    const double distance = MapperUtilities::ComputeDistance(this->Coordinates(), rInterfaceObject.Coordinates());

    const Node& r_node = *rInterfaceObject.pGetBaseNode();

    // Add a node candidate to the RBF support accumulator
    mRBFSupportAccumulator.AddNodeCandidate(r_node.GetValue(INTERFACE_EQUATION_ID), distance, r_node.Coordinates());

    // Reorder the support points according to the distance (starting from the closest)
    mRBFSupportAccumulator.ReorderSupportPoints();

    // If the accumulator has enough support points, stop searching 
    if (mRBFSupportAccumulator.HasEnough()) {
        SetLocalSearchWasSuccessful();
    } else if (mRBFSupportAccumulator.HasSome()) {
        SetIsApproximation();
    }

    KRATOS_CATCH("")

}

void RadialBasisFunctionMapperLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    KRATOS_TRY

    // If we don't have any interface info => nothing to do
    if (mInterfaceInfos.empty()) {
        ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
        return;
    }   

    const auto& rbf_support_accumulator = (mBuildOriginInterpolationMatrix && mOriginIsIga)
            ? static_cast<const RadialBasisFunctionsMapperInterfaceInfoIGA&>(*mInterfaceInfos[0]).GetRBFSupportAccumulator()
            : static_cast<const RadialBasisFunctionsMapperInterfaceInfoFEM&>(*mInterfaceInfos[0]).GetRBFSupportAccumulator();
    
    //  Get the support accumulator for this local system. We only have one InterfaceInfo per LocalSystem
    const auto& rbf_support_points = rbf_support_accumulator.Candidates();
    const IndexType rbf_support_points_number = rbf_support_accumulator.NumNeighbors();
    const Matrix rbf_support_points_coordinates = rbf_support_accumulator.GetCandidateCoordinatesMatrix();

    // Initialize the RBF Kernel for this mapping local system
    this->InitializeRBFKernel(rbf_support_points_coordinates);

    if (mBuildOriginInterpolationMatrix){
        // Initialize the local mapping matrix and ids vectors
        const IndexType num_rows = 1 + mNumberOfPolynomialTerms;
        const IndexType num_cols = rbf_support_points_number + mNumberOfPolynomialTerms;

        if (rLocalMappingMatrix.size1() != num_rows || rLocalMappingMatrix.size2() != num_cols) {
            rLocalMappingMatrix.resize(num_rows, num_cols, false);
        }
        // Resize IDs consistent with BuildMatrix
        rDestinationIds.resize(num_rows);
        rOriginIds.resize(num_cols);

        noalias(rLocalMappingMatrix) = ZeroMatrix(num_rows, num_cols);

        // Get the origin ids, first with the rbf part and secondly with the polynomial part
        for (IndexType j = 0; j < rbf_support_points_number; ++j) {
            rOriginIds[j] = rbf_support_points[j].mInterfaceEquationId;
        }
        for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
            const IndexType col = rbf_support_points_number + k;
            const IndexType poly_eq_id =(*mpPolynomialEquationIds)[k];;
            rOriginIds[col] = poly_eq_id;
        }
        
        // Get the destination ids
        if (mpNode)
            rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
        if (mpGeometry)
            rDestinationIds[0] = mpGeometry->GetValue(INTERFACE_EQUATION_ID);
        for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
            rDestinationIds[1 + k] = (*mpPolynomialEquationIds)[k];
        }
        
        for (IndexType j = 0; j < rbf_support_points_number; ++j) {
            const double distance = rbf_support_points[j].mDistance;
            const double phi = std::visit(
                [distance](const auto& kernel) {
                    return kernel(distance); 
                },
                mRBFKernel
            );

            rLocalMappingMatrix(0, j) = phi;
        }   

        const std::vector<double> polynomial = EvaluatePolynomialBasis(this->Coordinates(), mPolynomialDegree, mDimension);

        for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
            rLocalMappingMatrix(0, rbf_support_points_number + k) = polynomial[k];
        }
        
        // Build the polynomial constraint
        for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
            const IndexType row = 1 + k;
    
            // Put P_s[k] in the first polynomial column
            rLocalMappingMatrix(row, 0) = polynomial[k];
        }
        // for (IndexType j = 0; j < rbf_support_points_number; ++j) {
        //     const auto poly_j =
        //         EvaluatePolynomialBasis(rbf_support_points[j].mCoordinates, mPolynomialDegree);

        //     for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
        //         rLocalMappingMatrix(1 + k, j) = poly_j[k];
        //     }
        // }
    } else{
        // Initialize the local mapping matrix and ids vectors
        const IndexType num_rows = 1 ;
        const IndexType num_cols = rbf_support_points_number + mNumberOfPolynomialTerms;

        if (rLocalMappingMatrix.size1() != num_rows || rLocalMappingMatrix.size2() != num_cols) {
            rLocalMappingMatrix.resize(num_rows, num_cols, false);
        }
        // Resize IDs consistent with BuildMatrix
        rDestinationIds.resize(num_rows);
        rOriginIds.resize(num_cols);

        noalias(rLocalMappingMatrix) = ZeroMatrix(num_rows, num_cols);

        // Get the origin ids, first with the rbf part and secondly with the polynomial part
        for (IndexType j = 0; j < rbf_support_points_number; ++j) {
            rOriginIds[j] = rbf_support_points[j].mInterfaceEquationId;
        }
        for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
            const IndexType col = rbf_support_points_number + k;
            const IndexType poly_eq_id =(*mpPolynomialEquationIds)[k];;
            rOriginIds[col] = poly_eq_id;
        }
        
        // Get the destination ids
        if (mpNode)
            rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
        if (mpGeometry)
            rDestinationIds[0] = mpGeometry->GetValue(INTERFACE_EQUATION_ID);
        
        for (IndexType j = 0; j < rbf_support_points_number; ++j) {
            const double distance = rbf_support_points[j].mDistance;
            const double phi = std::visit(
                [distance](const auto& kernel) {
                    return kernel(distance); 
                },
                mRBFKernel
            );

            rLocalMappingMatrix(0, j) = phi;
        }

        const std::vector<double> polynomial = EvaluatePolynomialBasis(this->Coordinates(), mPolynomialDegree, mDimension);

        for (IndexType k = 0; k < mNumberOfPolynomialTerms; ++k) {
            rLocalMappingMatrix(0, rbf_support_points_number + k) = polynomial[k];
        }
    }
    KRATOS_CATCH("")
}

// std::vector<double> RadialBasisFunctionMapperLocalSystem::EvaluatePolynomialBasis(
//     const CoordinatesArrayType& rCoords, 
//     unsigned int degree) const
// {
//     std::vector<double> polynomial_values;
//     polynomial_values.reserve((degree + 3) * (degree + 2) * (degree + 1) / 6); // number of monomials in 3D

//     const double x = rCoords[0];
//     const double y = rCoords[1];
//     const double z = rCoords[2];

//     for (unsigned int n = 0; n <= degree; ++n) {
//         // i = exponente de x
//         for (int i = static_cast<int>(n); i >= 0; --i) {
//             // j = exponente de y
//             for (int j = static_cast<int>(n - i); j >= 0; --j) {
//                 unsigned int k = n - i - j; // exponente de z
//                 polynomial_values.push_back(std::pow(x, i) * std::pow(y, j) * std::pow(z, k));
//             }
//         }
//     }

//     return polynomial_values;
// }
std::vector<double> RadialBasisFunctionMapperLocalSystem::EvaluatePolynomialBasis(
    const CoordinatesArrayType& rCoords,
    unsigned int degree,
    const IndexType dimension) const
{
    std::vector<double> polynomial_values;

    const double x = rCoords[0];
    const double y = rCoords[1];
    const double z = rCoords[2];

    if (dimension == 2) {
        polynomial_values.reserve((degree + 2) * (degree + 1) / 2);

        for (unsigned int n = 0; n <= degree; ++n) {
            for (int i = static_cast<int>(n); i >= 0; --i) {
                const unsigned int j = n - i;
                polynomial_values.push_back(std::pow(x, i) * std::pow(y, j));
            }
        }
    } else if (dimension == 3) {
        polynomial_values.reserve((degree + 3) * (degree + 2) * (degree + 1) / 6);

        for (unsigned int n = 0; n <= degree; ++n) {
            for (int i = static_cast<int>(n); i >= 0; --i) {
                for (int j = static_cast<int>(n - i); j >= 0; --j) {
                    const unsigned int k = n - i - j;
                    polynomial_values.push_back(std::pow(x, i) * std::pow(y, j) * std::pow(z, k));
                }
            }
        }
    } else {
        KRATOS_ERROR << "Unsupported polynomial dimension: " << dimension << std::endl;
    }

    return polynomial_values;
}

template<class TSparseSpace, class TDenseSpace>
RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::RadialBasisFunctionMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters)
    :  mrModelPartOrigin(rModelPartOrigin),
       mrModelPartDestination(rModelPartDestination),
       mMapperSettings(JsonParameters)
{
    mMapperSettings.ValidateAndAssignDefaults(this->GetMapperDefaultSettings());

    mOriginIsIga = mMapperSettings["origin_is_iga"].GetBool();

    // Define the origin and destination interfaces
    mpCouplingInterfaceOrigin = &rModelPartOrigin;
    mpCouplingInterfaceDestination = &rModelPartDestination;
    
    // Construct the InterfaceVectorContainers for the origin and destination domain
    mpInterfaceVectorContainerOrigin = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceOrigin);
    mpInterfaceVectorContainerDestination = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceDestination);
    
    this->ValidateInput();
    this->CreateLinearSolver();
    this->InitializeInterface();
}


template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    const IndexType num_nodes_interface_destination = mpCouplingInterfaceDestination->NumberOfNodes();
    const IndexType num_nodes_interface_origin = mpCouplingInterfaceOrigin->NumberOfNodes();
    const IndexType num_elements_interface_origin = mpCouplingInterfaceOrigin->NumberOfElements();
    mpMappingMatrix = Kratos::make_unique<MappingMatrixType>(num_nodes_interface_destination, num_nodes_interface_origin);

    const bool precompute_mapping_matrix = mMapperSettings["precompute_mapping_matrix"].GetBool();

    // Read the RBF type from the settings
    std::string rbf_type = mMapperSettings["radial_basis_function_type"].GetString();
    std::transform(rbf_type.begin(), rbf_type.end(), rbf_type.begin(), ::tolower);

    static const std::unordered_map<std::string, RadialBasisFunctionsUtilities::RBFType> rbf_type_map = {
        {"inverse_multiquadric", RadialBasisFunctionsUtilities::RBFType::InverseMultiquadric},
        {"multiquadric",         RadialBasisFunctionsUtilities::RBFType::Multiquadric},
        {"gaussian",             RadialBasisFunctionsUtilities::RBFType::Gaussian},
        {"thin_plate_spline",    RadialBasisFunctionsUtilities::RBFType::ThinPlateSpline},
        {"wendland_c2",          RadialBasisFunctionsUtilities::RBFType::WendlandC2}
    };

    const auto it = rbf_type_map.find(rbf_type);
    if (it == rbf_type_map.end()) {
        KRATOS_ERROR << "Unrecognized RBF type: " << rbf_type << std::endl;
    }

    mRBFTypeEnum = it->second;

    // Obtain the polynomial degree and calculate number of polynomial terms
    mPolynomialDegree = mMapperSettings["additional_polynomial_degree"].GetInt();
    mDimension = mMapperSettings["global_space_dimension"].GetInt();
    KRATOS_ERROR_IF(mDimension > 3) << "global_space_dimension cannot be greater than 3" << std::endl;

    // Calculate the number of polynomial terms depending on the dimension
     if (mDimension == 2){
        mNumberOfPolynomialTerms = (mPolynomialDegree + 2) * (mPolynomialDegree + 1) / 2;
    } else {
        mNumberOfPolynomialTerms = (mPolynomialDegree + 3) * (mPolynomialDegree + 2) * (mPolynomialDegree + 1) / 6;
    }

    // Determine whether the destination domain is IGA or not
    const bool destination_is_iga = mMapperSettings["destination_is_iga"].GetBool();

    // Throw an error if the destination domain is discretized with IBRA
    KRATOS_ERROR_IF(destination_is_iga)<< "This mapper is not yet available when the destination domain is discretized with IGA" << std::endl;

    bool use_all_support_points = false;

    if (mMapperSettings.Has("use_all_rbf_support_points")) {
        use_all_support_points =
            mMapperSettings["use_all_rbf_support_points"].GetBool();
    }

    IndexType max_support_points = 0;
    if (!use_all_support_points) {
        KRATOS_ERROR_IF_NOT(mMapperSettings.Has("max_support_points"))
            << "Missing 'max_support_points' in search_settings when 'use_all_rbf_support_points' is false." << std::endl;

        max_support_points = mMapperSettings["max_support_points"].GetInt();
        mRequiredRBFSupportPoints = max_support_points;
    } else if (use_all_support_points && mOriginIsIga) {
        mRequiredRBFSupportPoints = num_elements_interface_origin;
    } else {
        mRequiredRBFSupportPoints = num_nodes_interface_origin;
    }

    if(mOriginIsIga){
        KRATOS_ERROR_IF(mRequiredRBFSupportPoints > num_elements_interface_origin) << 
            "The required support points for the RBF Mapping cannot excede the number of integrations points in the IGA origin domain" << std::endl;
        // KRATOS_ERROR_IF(!precompute_mapping_matrix) << 
        //     "The RBF Mapper for IGA only works by precomputing the mapping matrix" << std::endl;
    } else{
        KRATOS_ERROR_IF(mRequiredRBFSupportPoints > num_nodes_interface_origin) << 
            "The required support points for the RBF Mapping cannot excede the number of nodes in the FEM origin domain" << std::endl;
    }

    // Assign the INTERFACE_EQUATION_IDs required to build the interpolation matrices
    if (!mOriginIsIga){
        AssignInterfaceEquationIds();
    } else{
        AssignInterfaceEquationIdsIga();
    }

    // Initialize the polynomial equation ids
    if (!mOriginIsIga){
        mPolynomialEquationIdsOrigin = InitializePolynomialEquationIds(mpCouplingInterfaceOrigin);
    } else {
        mPolynomialEquationIdsOrigin = InitializePolynomialEquationIdsIga(mpCouplingInterfaceOrigin);
    }

    if (mOriginIsIga) {
        mpNReduced = ComputeNReduced(this->GetOriginModelPart());
    }

    if (!use_all_support_points){
        // Create dummy local systems for the origin and destination domains
        RadialBasisFunctionMapperLocalSystem origin_local_system{nullptr, nullptr, true, mOriginIsIga, mRBFTypeEnum, mPolynomialDegree, mDimension, &mPolynomialEquationIdsOrigin};
        RadialBasisFunctionMapperLocalSystem destination_local_system{nullptr, nullptr, false, mOriginIsIga, mRBFTypeEnum, mPolynomialDegree, mDimension, &mPolynomialEquationIdsOrigin};

        // Create the local systems for the origin and destination domains
        if (mOriginIsIga) {
            MapperUtilities::CreateMapperLocalSystemsFromGeometries(
                origin_local_system,
                mrModelPartOrigin.GetCommunicator(),
                mMapperLocalSystemsOrigin);
        } else {
            MapperUtilities::CreateMapperLocalSystemsFromNodes(
                origin_local_system,
                mrModelPartOrigin.GetCommunicator(),
                mMapperLocalSystemsOrigin);
        }

        // The destination is always FEM
        MapperUtilities::CreateMapperLocalSystemsFromNodes(
            destination_local_system,
            mrModelPartDestination.GetCommunicator(),
            mMapperLocalSystemsDestination);

        auto p_interface_comm_origin_with_origin = Kratos::make_unique<InterfaceCommunicator>(
                mrModelPartOrigin,
                mMapperLocalSystemsOrigin,
                mMapperSettings["search_settings"]);
        
        auto p_interface_comm_origin_with_destination = Kratos::make_unique<InterfaceCommunicator>(
            mrModelPartOrigin,
            mMapperLocalSystemsDestination,
            mMapperSettings["search_settings"]);
        
        MapperInterfaceInfoUniquePointerType p_ref_interface_info_origin_with_origin;
        MapperInterfaceInfoUniquePointerType p_ref_interface_info_origin_with_destination;
        

        p_ref_interface_info_origin_with_origin = GetMapperInterfaceInfo(mRequiredRBFSupportPoints);
        p_ref_interface_info_origin_with_destination = GetMapperInterfaceInfo(mRequiredRBFSupportPoints);
        
        // Perform the local search to find the N nearest support points for each local system (The FEM local systems are always centered on the nodes while the IGA local systems are always
        // centered on the integration points (geometries))
        p_interface_comm_origin_with_origin->ExchangeInterfaceData(mrModelPartOrigin.GetCommunicator(),
                                                        p_ref_interface_info_origin_with_origin);
        p_interface_comm_origin_with_destination->ExchangeInterfaceData(mrModelPartDestination.GetCommunicator(),
                                                        p_ref_interface_info_origin_with_destination);
        
        // Build the origin interpolation matrix
        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrixRBFMapper(
            mpOriginInterpolationMatrix,
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mMapperLocalSystemsOrigin,
            mNumberOfPolynomialTerms,
            true,
            mOriginIsIga,
            0);
        KRATOS_WATCH("I AM DONE WITH THE ORIGIN")
        
        // Build the destination evaluation matrix
        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrixRBFMapper(
            mpDestinationEvaluationMatrix,
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerDestination->pGetVector(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mpInterfaceVectorContainerDestination->GetModelPart(),
            mMapperLocalSystemsDestination,
            mNumberOfPolynomialTerms,
            false,
            mOriginIsIga,
            0);
        KRATOS_WATCH("I AM DONE WITH THE DESTINATION")
    } else {
        using CandidateType = typename RBFSupportAccumulator::Candidate;

        // Collect all origin support points
        std::vector<CandidateType> origin_support_points;
        origin_support_points.reserve(mRequiredRBFSupportPoints);

        if (!mOriginIsIga) {
            for (auto& r_node : mpCouplingInterfaceOrigin->Nodes()) {
                CandidateType candidate;
                candidate.mDistance = 0.0;
                candidate.mInterfaceEquationId = r_node.GetValue(INTERFACE_EQUATION_ID);
                candidate.mCoordinates = r_node.Coordinates();
                origin_support_points.push_back(candidate);
            }
        } else {
            for (auto& r_elem : mpCouplingInterfaceOrigin->Elements()) {
                auto& r_geom = r_elem.GetGeometry();

                CandidateType candidate;
                candidate.mDistance = 0.0;
                candidate.mInterfaceEquationId = r_elem.GetValue(INTERFACE_EQUATION_ID);
                candidate.mCoordinates = r_geom.Center();
                origin_support_points.push_back(candidate);
            }
        }

        KRATOS_ERROR_IF(origin_support_points.size() != mRequiredRBFSupportPoints)
            << "Unexpected number of origin support points. Expected "
            << mRequiredRBFSupportPoints << " but got "
            << origin_support_points.size() << std::endl;

        // Collect all destination points (always FEM nodes)
        std::vector<CandidateType> destination_points;
        destination_points.reserve(mpCouplingInterfaceDestination->NumberOfNodes());

        for (auto& r_node : mpCouplingInterfaceDestination->Nodes()) {
            CandidateType candidate;
            candidate.mDistance = 0.0;
            candidate.mInterfaceEquationId = r_node.GetValue(INTERFACE_EQUATION_ID);
            candidate.mCoordinates = r_node.Coordinates();
            destination_points.push_back(candidate);
        }

        // Create local systems exactly as in the search branch
        RadialBasisFunctionMapperLocalSystem origin_local_system{
            nullptr, nullptr, true, mOriginIsIga, mRBFTypeEnum,
            mPolynomialDegree, mDimension, &mPolynomialEquationIdsOrigin
        };

        RadialBasisFunctionMapperLocalSystem destination_local_system{
            nullptr, nullptr, false, mOriginIsIga, mRBFTypeEnum,
            mPolynomialDegree, mDimension, &mPolynomialEquationIdsOrigin
        };

        if (mOriginIsIga) {
            MapperUtilities::CreateMapperLocalSystemsFromGeometries(
                origin_local_system,
                mrModelPartOrigin.GetCommunicator(),
                mMapperLocalSystemsOrigin);
        } else {
            MapperUtilities::CreateMapperLocalSystemsFromNodes(
                origin_local_system,
                mrModelPartOrigin.GetCommunicator(),
                mMapperLocalSystemsOrigin);
        }

        MapperUtilities::CreateMapperLocalSystemsFromNodes(
            destination_local_system,
            mrModelPartDestination.GetCommunicator(),
            mMapperLocalSystemsDestination);

        // Manually attach all-support interface info to ORIGIN local systems
        for (auto& rp_local_system : mMapperLocalSystemsOrigin) {
            auto p_info = GetMapperInterfaceInfo(mRequiredRBFSupportPoints);

                if (mOriginIsIga) {
                    auto& r_info = static_cast<RadialBasisFunctionsMapperInterfaceInfoIGA&>(*p_info);
                    auto& r_acc = r_info.GetRBFSupportAccumulator();

                    const auto local_coords = rp_local_system->Coordinates();

                    for (const auto& r_candidate : origin_support_points) {
                        const double distance = MapperUtilities::ComputeDistance(
                            local_coords,
                            r_candidate.mCoordinates);

                        r_acc.AddGeometryCandidate(
                            r_candidate.mInterfaceEquationId,
                            distance,
                            r_candidate.mCoordinates);
                    }

                    r_acc.ReorderSupportPoints();
                } else {
                auto& r_info = static_cast<RadialBasisFunctionsMapperInterfaceInfoFEM&>(*p_info);
                auto& r_acc = r_info.GetRBFSupportAccumulator();

                const auto local_coords = rp_local_system->Coordinates();

                for (const auto& r_candidate : origin_support_points) {
                    const double distance = MapperUtilities::ComputeDistance(
                        local_coords,
                        r_candidate.mCoordinates);

                    r_acc.AddNodeCandidate(
                        r_candidate.mInterfaceEquationId,
                        distance,
                        r_candidate.mCoordinates);
                }
                r_acc.ReorderSupportPoints();
            }

            rp_local_system->AddInterfaceInfo(std::move(p_info));
        }

        // Manually attach all-support interface info to DESTINATION local systems
        for (auto& rp_local_system : mMapperLocalSystemsDestination) {
            auto p_info = GetMapperInterfaceInfo(mRequiredRBFSupportPoints);
            const auto local_coords = rp_local_system->Coordinates();

            if (mOriginIsIga) {
                auto& r_info = static_cast<RadialBasisFunctionsMapperInterfaceInfoIGA&>(*p_info);
                auto& r_acc = r_info.GetRBFSupportAccumulator();

                for (const auto& r_candidate : origin_support_points) {
                    const double distance = MapperUtilities::ComputeDistance(
                        local_coords,
                        r_candidate.mCoordinates);

                    r_acc.AddGeometryCandidate(
                        r_candidate.mInterfaceEquationId,
                        distance,
                        r_candidate.mCoordinates);
                }
                r_acc.ReorderSupportPoints();
            } else {
                auto& r_info = static_cast<RadialBasisFunctionsMapperInterfaceInfoFEM&>(*p_info);
                auto& r_acc = r_info.GetRBFSupportAccumulator();

                for (const auto& r_candidate : origin_support_points) {
                    const double distance = MapperUtilities::ComputeDistance(
                        local_coords,
                        r_candidate.mCoordinates);

                    r_acc.AddNodeCandidate(
                        r_candidate.mInterfaceEquationId,
                        distance,
                        r_candidate.mCoordinates);
                }
                r_acc.ReorderSupportPoints();
            }

            rp_local_system->AddInterfaceInfo(std::move(p_info));
        }

        // Reuse your existing matrix builders unchanged
        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrixRBFMapper(
            mpOriginInterpolationMatrix,
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mMapperLocalSystemsOrigin,
            mNumberOfPolynomialTerms,
            true,
            mOriginIsIga,
            0);

        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrixRBFMapper(
            mpDestinationEvaluationMatrix,
            mpInterfaceVectorContainerOrigin->pGetVector(),
            mpInterfaceVectorContainerDestination->pGetVector(),
            mpInterfaceVectorContainerOrigin->GetModelPart(),
            mpInterfaceVectorContainerDestination->GetModelPart(),
            mMapperLocalSystemsDestination,
            mNumberOfPolynomialTerms,
            false,
            mOriginIsIga,
            0);
        

        KRATOS_INFO("RBFMapper") << "All-support path: matrices built without local search." << std::endl;
    }
    
    if (precompute_mapping_matrix){
        CalculateMappingMatrix();
    }
}

template<class TSparseSpace, class TDenseSpace>
std::vector<IndexType> RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::InitializePolynomialEquationIds(const ModelPart* pModelPart)
{
    // Use the local mesh nodes you already grabbed
    const auto& r_nodes = pModelPart->GetCommunicator().LocalMesh().Nodes();

    IndexType max_eq_id = block_for_each<MaxReduction<IndexType>>(r_nodes, 
        [](const Node& rNode){                          
            if (rNode.Has(INTERFACE_EQUATION_ID)) {
                return rNode.GetValue(INTERFACE_EQUATION_ID);
            }
            return 0;
        });

    // Fill polynomial_equation_ids as a simple consecutive range
    std::vector<IndexType> polynomial_equation_ids(mNumberOfPolynomialTerms);
    std::iota(polynomial_equation_ids.begin(),
              polynomial_equation_ids.end(),
              max_eq_id + 1);

    return polynomial_equation_ids;
}

template<class TSparseSpace, class TDenseSpace>
std::vector<IndexType> RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::InitializePolynomialEquationIdsIga(const ModelPart* pModelPart)
{
    // Use the local mesh nodes you already grabbed
    const auto& r_elements = pModelPart->GetCommunicator().LocalMesh().Elements();

    IndexType max_eq_id = block_for_each<MaxReduction<IndexType>>(r_elements, 
        [](const Element& rCondition){                          
            if (rCondition.Has(INTERFACE_EQUATION_ID)) {
                return rCondition.GetValue(INTERFACE_EQUATION_ID);
            }
            return 0;
        });

    // Fill polynomial_equation_ids as a simple consecutive range
    std::vector<IndexType> polynomial_equation_ids(mNumberOfPolynomialTerms);
    std::iota(polynomial_equation_ids.begin(),
              polynomial_equation_ids.end(),
              max_eq_id + 1);

    return polynomial_equation_ids;
}

// template<class TSparseSpace, class TDenseSpace>
// void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MapInternal(
//     const Variable<double>& rOriginVariable,
//     const Variable<double>& rDestinationVariable,
//     Flags MappingOptions)
// {
//     const bool precompute_mapping_matrix =
//         mMapperSettings["precompute_mapping_matrix"].GetBool();

//     // 1. Read origin vector (nodes for FEM, control points for IGA)
//     mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

//     auto& r_origin_vec      = mpInterfaceVectorContainerOrigin->GetVector();
//     auto& r_destination_vec = mpInterfaceVectorContainerDestination->GetVector();

//     if (precompute_mapping_matrix) {
//         KRATOS_ERROR_IF(!mpMappingMatrix) << "MapInternal: mapping matrix is not computed but 'precompute_mapping_matrix' is true.\n";

//         r_destination_vec.resize(mpMappingMatrix->size1(), false);
        
//         // u_dest = M * u_orig
//         TSparseSpace::Mult(
//             *mpMappingMatrix, 
//             r_origin_vec, 
//             r_destination_vec);
//     } else {
//        // solve RBF system on the fly 
//         const IndexType system_size  = mpOriginInterpolationMatrix->size1();
//         const IndexType n_origin_dof = r_origin_vec.size();
//         const IndexType n_poly       = system_size - n_origin_dof;

//         Vector rhs(system_size);
//         for (IndexType i = 0; i < n_origin_dof; ++i)
//             rhs[i] = r_origin_vec[i];
//         for (IndexType i = 0; i < n_poly; ++i)
//             rhs[n_origin_dof + i] = 0.0;

//         Vector solution(system_size);
//         mpLinearSolver->Solve(*mpOriginInterpolationMatrix, solution, rhs);

//         TSparseSpace::Mult(
//             *mpDestinationEvaluationMatrix, 
//             solution, 
//             r_destination_vec);
//     }

//     mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
// }
template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Flags MappingOptions)
{
    const bool precompute_mapping_matrix =
        mMapperSettings["precompute_mapping_matrix"].GetBool();

    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    auto& r_origin_vec      = mpInterfaceVectorContainerOrigin->GetVector();
    auto& r_destination_vec = mpInterfaceVectorContainerDestination->GetVector();

    if (precompute_mapping_matrix) {
        KRATOS_ERROR_IF(!mpMappingMatrix) << "MapInternal: mapping matrix is not computed but 'precompute_mapping_matrix' is true.\n";

        r_destination_vec.resize(mpMappingMatrix->size1(), false);
        TSparseSpace::SetToZero(r_destination_vec);
        TSparseSpace::Mult(*mpMappingMatrix, r_origin_vec, r_destination_vec);
    } else {
        KRATOS_ERROR_IF_NOT(mpOriginInterpolationMatrix) << "mpOriginInterpolationMatrix is null" << std::endl;
        KRATOS_ERROR_IF_NOT(mpDestinationEvaluationMatrix) << "mpDestinationEvaluationMatrix is null" << std::endl;

        const IndexType system_size = mpOriginInterpolationMatrix->size1();

        Vector rhs(system_size, 0.0);

        if (!mOriginIsIga) {
            const IndexType n_origin_dof = r_origin_vec.size();

            KRATOS_ERROR_IF(n_origin_dof + mNumberOfPolynomialTerms != system_size)
                << "FEM on-the-fly mapping: inconsistent sizes. n_origin_dof = "
                << n_origin_dof << ", n_poly = " << mNumberOfPolynomialTerms
                << ", system_size = " << system_size << std::endl;

            for (IndexType i = 0; i < n_origin_dof; ++i) {
                rhs[i] = r_origin_vec[i];
            }
        } else {
            KRATOS_ERROR_IF_NOT(mpNReduced) << "mpNReduced not initialized for IGA on-the-fly mapping." << std::endl;

            KRATOS_ERROR_IF(mpNReduced->size2() != r_origin_vec.size())
                << "IGA on-the-fly mapping: size mismatch in N_reduced * u_cp. "
                << "N_reduced size2 = " << mpNReduced->size2()
                << ", origin vector size = " << r_origin_vec.size() << std::endl;

            Vector gp_values(mpNReduced->size1(), 0.0);
            TSparseSpace::Mult(*mpNReduced, r_origin_vec, gp_values);

            KRATOS_ERROR_IF(gp_values.size() + mNumberOfPolynomialTerms != system_size)
                << "IGA on-the-fly mapping: inconsistent sizes. n_gp = "
                << gp_values.size() << ", n_poly = " << mNumberOfPolynomialTerms
                << ", system_size = " << system_size << std::endl;

            for (IndexType i = 0; i < gp_values.size(); ++i) {
                rhs[i] = gp_values[i];
            }
        }

        Vector solution(system_size, 0.0);
        KRATOS_WATCH("STARTING THE SOLVE")
        mpLinearSolver->Solve(*mpOriginInterpolationMatrix, solution, rhs);
        KRATOS_WATCH("I AM DONE WITH SOLVE")

        r_destination_vec.resize(mpDestinationEvaluationMatrix->size1(), false);
        TSparseSpace::SetToZero(r_destination_vec);
        TSparseSpace::Mult(*mpDestinationEvaluationMatrix, solution, r_destination_vec);
    }

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Flags MappingOptions)
{
    const bool precompute_mapping_matrix =
        mMapperSettings["precompute_mapping_matrix"].GetBool();

    // Read destination variable from the FEM nodes
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    auto& r_dest_vec   = mpInterfaceVectorContainerDestination->GetVector();
    auto& r_origin_vec = mpInterfaceVectorContainerOrigin->GetVector();

    if (precompute_mapping_matrix) {
        r_origin_vec.resize(mpMappingMatrix->size2(), false);

        // f_orig = M^T * f_dest
        TSparseSpace::TransposeMult(
            *mpMappingMatrix, 
            r_dest_vec, 
            r_origin_vec);
    } else {
        // solve RBF system on the fly 
        const IndexType system_size = mpOriginInterpolationMatrix->size1();

        Vector temp_vector_1(system_size);
        TSparseSpace::TransposeMult(
            *mpDestinationEvaluationMatrix,
            r_dest_vec,
            temp_vector_1
        );

        Vector temp_vector_2(system_size);
        mpLinearSolver->Solve(*mpOriginInterpolationMatrix, temp_vector_2, temp_vector_1);

        const IndexType n_origin_dof = r_origin_vec.size();
        r_origin_vec.resize(n_origin_dof, false);
        for (IndexType i = 0; i < n_origin_dof; ++i)
            r_origin_vec[i] = temp_vector_2[i];
    }

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
}


template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    Flags MappingOptions)
{
    for (const auto var_ext : {"_X", "_Y", "_Z"}) {
        const auto& var_origin = KratosComponents<Variable<double>>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<Variable<double>>::Get(rDestinationVariable.Name() + var_ext);

        MapInternal(var_origin, var_destination, MappingOptions);
    }
}

template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<array_1d<double,3>>& rOriginVariable,
    const Variable<array_1d<double,3>>& rDestinationVariable,
    Flags MappingOptions)
{
    for (const auto var_ext : {"_X", "_Y", "_Z"}) {
        const auto& var_origin = KratosComponents<Variable<double>>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<Variable<double>>::Get(rDestinationVariable.Name() + var_ext);

        MapInternalTranspose(var_origin, var_destination, MappingOptions);
    }
}

#include <chrono>
#include <iomanip>

// ...

template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CalculateMappingMatrix()
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpOriginInterpolationMatrix)
        << "CalculateMappingMatrix: mpOriginInterpolationMatrix is nullptr." << std::endl;
    KRATOS_DEBUG_ERROR_IF_NOT(mpDestinationEvaluationMatrix)
        << "CalculateMappingMatrix: mpDestinationEvaluationMatrix is nullptr." << std::endl;

    using Clock = std::chrono::steady_clock;
    using Duration = std::chrono::duration<double>;

    // Dense matrix type for multi-RHS solve
    //using DenseMatrixType = boost::numeric::ublas::matrix<double>;
    using DenseMatrixType = typename TDenseSpace::MatrixType;

    const auto t_total_begin = Clock::now();

    const IndexType system_size = mpOriginInterpolationMatrix->size1(); // n_support + n_poly
    const IndexType n_dest      = mpDestinationEvaluationMatrix->size1();

    IndexType n_support = 0;
    if (!mOriginIsIga) {
        n_support = mpInterfaceVectorContainerOrigin->GetVector().size();
    } else {
        n_support = mpCouplingInterfaceOrigin->NumberOfElements();
    }

    KRATOS_INFO("RBFMapper") << "system_size = " << system_size
                             << ", n_support = " << n_support
                             << ", n_dest = " << n_dest << std::endl;

    double alloc_time        = 0.0;
    double rhs_build_time    = 0.0;
    double solve_time        = 0.0;
    double mult_time         = 0.0;
    double convert_time      = 0.0;
    double postprocess_time  = 0.0;

    // -------------------------------------------------
    // 1) Allocate dense RHS and solution matrices
    // -------------------------------------------------
    auto t0 = Clock::now();

    DenseMatrixType rhs_mat(system_size, n_support, 0.0);
    DenseMatrixType sol_mat(system_size, n_support, 0.0);

    auto t1 = Clock::now();
    alloc_time = Duration(t1 - t0).count();

    // -------------------------------------------------
    // 2) Build RHS = [I; 0]
    // -------------------------------------------------
    auto t2 = Clock::now();

    for (IndexType j = 0; j < n_support; ++j) {
        rhs_mat(j, j) = 1.0;
    }

    auto t3 = Clock::now();
    rhs_build_time = Duration(t3 - t2).count();

    // -------------------------------------------------
    // 3) Solve A * X = RHS
    // -------------------------------------------------
    auto t4 = Clock::now();

    // IMPORTANT:
    // This requires that the concrete linear solver supports
    // Solve(SparseMatrixType, DenseMatrixType, DenseMatrixType)
    mpLinearSolver->Solve(*mpOriginInterpolationMatrix, sol_mat, rhs_mat);
    KRATOS_WATCH("SOLVE DONE")

    // // -----------------------------------------------------------------
    // using SparseMatrixType = boost::numeric::ublas::compressed_matrix<double>;
    // using DenseMatrixType  = boost::numeric::ublas::matrix<double>;

    // const SparseMatrixType& A_sparse = *mpOriginInterpolationMatrix;

    // KRATOS_ERROR_IF(A_sparse.size1() != A_sparse.size2())
    //     << "Interpolation matrix is not square! "
    //     << A_sparse.size1() << " x " << A_sparse.size2() << std::endl;

    // DenseMatrixType A(A_sparse.size1(), A_sparse.size2(), 0.0);

    // for (auto it1 = A_sparse.begin1(); it1 != A_sparse.end1(); ++it1) {
    //     for (auto it2 = it1.begin(); it2 != it1.end(); ++it2) {
    //         A(it2.index1(), it2.index2()) = *it2;
    //     }
    // }

    // DenseMatrixType A_inv(A.size1(), A.size2());
    // double det = 0.0;

    // MathUtils<double>::InvertMatrix(A, A_inv, det);

    // const double cond_est = norm_frobenius(A) * norm_frobenius(A_inv);

    // KRATOS_INFO("RBF Mapper") << "det(A)    = " << det << std::endl;
    // KRATOS_INFO("RBF Mapper") << "cond_F(A) = " << cond_est << std::endl;
    //---------------------------------------------------------------

    auto t5 = Clock::now();
    solve_time = Duration(t5 - t4).count();

    // -------------------------------------------------
    // 4) Compute H_tilde = B * X
    // -------------------------------------------------
    auto t6 = Clock::now();

    auto pMappingMatrixGP = Kratos::make_unique<MappingMatrixType>(n_dest, n_support);

    // Since MappingMatrixType is compressed_matrix in your setup,
    // compute first into a dense matrix and then convert.
    // DenseMatrixType mapping_dense(n_dest, n_support, 0.0);

    // noalias(mapping_dense) = prod(*mpDestinationEvaluationMatrix, sol_mat);
    DenseMatrixType mapping_dense(n_dest, n_support, 0.0);

    typename TSparseSpace::VectorType sol_col;
    typename TSparseSpace::VectorType dest_col;

    TSparseSpace::Resize(sol_col, system_size);
    TSparseSpace::Resize(dest_col, n_dest);

    KRATOS_WATCH("STARTING MULTIPLICATION")

    for (IndexType j = 0; j < n_support; ++j) {
        for (IndexType i = 0; i < system_size; ++i) {
            sol_col[i] = sol_mat(i, j);
        }

        TSparseSpace::SetToZero(dest_col);
        TSparseSpace::Mult(*mpDestinationEvaluationMatrix, sol_col, dest_col);

        for (IndexType i = 0; i < n_dest; ++i) {
            mapping_dense(i, j) = dest_col[i];
        }
    }
    KRATOS_WATCH("FINISHING MULTIPLICATION")

    auto t7 = Clock::now();
    mult_time = Duration(t7 - t6).count();

    // -------------------------------------------------
    // 5) Convert dense mapping to MappingMatrixType
    // -------------------------------------------------
    auto t8 = Clock::now();

    for (IndexType i = 0; i < n_dest; ++i) {
        for (IndexType j = 0; j < n_support; ++j) {
            const double value = mapping_dense(i, j);
            if (std::abs(value) > 1e-16) {
                (*pMappingMatrixGP)(i, j) = value;
            }
        }
    }

    auto t9 = Clock::now();
    convert_time = Duration(t9 - t8).count();

    KRATOS_INFO("RBFMapper") << std::fixed << std::setprecision(6)
        << "\n=== Timing summary for CalculateMappingMatrix (matrix RHS) ===\n"
        << "Allocation time   : " << alloc_time       << " s\n"
        << "RHS build time    : " << rhs_build_time   << " s\n"
        << "Solve time        : " << solve_time       << " s\n"
        << "B * X time        : " << mult_time        << " s\n"
        << "Convert time      : " << convert_time     << " s"
        << std::endl;

    // -------------------------------------------------
    // 6) Final postprocessing
    // -------------------------------------------------
    auto t10 = Clock::now();

    if (mOriginIsIga) {
        mpMappingMatrix = ComputeMappingMatrixIgaOnControlPoints(
            *pMappingMatrixGP,
            this->GetOriginModelPart()
        );
    } else {
        mpMappingMatrix = std::move(pMappingMatrixGP);
    }

    auto t11 = Clock::now();
    postprocess_time = Duration(t11 - t10).count();

    const auto t_total_end = Clock::now();

    KRATOS_INFO("RBFMapper") << "Postprocessing time : "
                             << postprocess_time << " s" << std::endl;

    KRATOS_INFO("RBFMapper") << "Total CalculateMappingMatrix time : "
                             << Duration(t_total_end - t_total_begin).count() << " s" << std::endl;
}

template<class TSparseSpace, class TDenseSpace>
std::unique_ptr<typename RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MappingMatrixType>
RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::ComputeMappingMatrixIgaOnControlPoints(
    const MappingMatrixType& rMappingMatrixGP,
    const ModelPart& rOriginModelPart) const
{
    using MappingMatrixType = typename RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MappingMatrixType;

    // Number of "origin interface objects" = number of GPs (here: 1 GP per condition)
    const IndexType n_gp = mpCouplingInterfaceOrigin->NumberOfElements();

    // 1) Build map: NodeId -> local control-point index (column index in N_reduced)
    std::unordered_map<IndexType, IndexType> node_id_to_cp_index;
    node_id_to_cp_index.reserve(mpCouplingInterfaceOrigin->NumberOfNodes());

    IndexType n_cp = 0;
    for (auto& r_node : mpCouplingInterfaceOrigin->Nodes()) {
        const IndexType node_id = r_node.Id();
        node_id_to_cp_index.emplace(node_id, n_cp++);
    }

    // 2) Build N_reduced (n_gp x n_cp) with shape functions of interface GPs
    Matrix N_reduced(n_gp, n_cp);
    noalias(N_reduced) = ZeroMatrix(n_gp, n_cp);

    IndexType gp_idx = 0;
    for (auto& r_cond : mpCouplingInterfaceOrigin->Elements())
    {
        auto& r_geom = r_cond.GetGeometry();

        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues();

        const auto& r_N = row(Ncontainer, 0);
        const IndexType num_points = r_geom.PointsNumber();

        for (IndexType j = 0; j < num_points; ++j) {
            const auto& r_node   = r_geom[j];
            const IndexType node_id = r_node.Id();

            auto it = node_id_to_cp_index.find(node_id);
            const IndexType cp_col = it->second;

            N_reduced(gp_idx, cp_col) = r_N[j];
        }

        ++gp_idx;
    }

    // 3) Compute final mapping: M = M_tilde * N_reduced
    //    rMappingMatrixGP: (n_dest x n_gp)
    //    N_reduced      : (n_gp  x n_cp)
    //    -> H           : (n_dest x n_cp)
    auto pMappingMatrix = Kratos::make_unique<MappingMatrixType>(
        rMappingMatrixGP.size1(),  // n_dest
        n_cp                       // number of control points on origin interface
    );

    noalias(*pMappingMatrix) = prod(rMappingMatrixGP, N_reduced);

    return pMappingMatrix;
}



template<class TSparseSpace, class TDenseSpace>
void RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    // Create the linear solver factory
    LinearSolverFactory<TSparseSpace, TDenseSpace> solver_factory;

    if (mMapperSettings["linear_solver_settings"].Has("solver_type")) {
        const std::string solver_type =
            mMapperSettings["linear_solver_settings"]["solver_type"].GetString();

        KRATOS_INFO("Radial Basis Functions Mapper")
            << "Using specified linear solver: " << solver_type << std::endl;

        mpLinearSolver = solver_factory.Create(mMapperSettings["linear_solver_settings"]);
    }
    else {
        if (solver_factory.Has("pardiso_lu")) {
            KRATOS_INFO("Radial Basis Functions Mapper")
                << "No linear solver specified. Using PardisoLU." << std::endl;

            Parameters default_settings(R"(
            {
                "solver_type": "pardiso_lu"
            })");

            mpLinearSolver = solver_factory.Create(default_settings);
        }
        else {
            KRATOS_WARNING("Radial Basis Functions Mapper")
                << "PardisoLU not available. Falling back to skyline_lu_factorization." << std::endl;

            Parameters default_settings(R"(
            {
                "solver_type": "skyline_lu_factorization"
            })");

            mpLinearSolver = solver_factory.Create(default_settings);
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
std::unique_ptr<typename RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MappingMatrixType>
RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::ComputeNReduced(const ModelPart& rOriginModelPart) const
{
    using MappingMatrixType = typename RadialBasisFunctionMapper<TSparseSpace, TDenseSpace>::MappingMatrixType;

    const IndexType n_gp = mpCouplingInterfaceOrigin->NumberOfElements();

    std::unordered_map<IndexType, IndexType> node_id_to_cp_index;
    node_id_to_cp_index.reserve(mpCouplingInterfaceOrigin->NumberOfNodes());

    IndexType n_cp = 0;
    for (auto& r_node : mpCouplingInterfaceOrigin->Nodes()) {
        node_id_to_cp_index.emplace(r_node.Id(), n_cp++);
    }

    auto pNReduced = Kratos::make_unique<MappingMatrixType>(n_gp, n_cp);

    IndexType gp_idx = 0;
    for (auto& r_elem : mpCouplingInterfaceOrigin->Elements())
    {
        auto& r_geom = r_elem.GetGeometry();
        const Matrix& Ncontainer = r_geom.ShapeFunctionsValues();
        const auto& r_N = row(Ncontainer, 0);
        const IndexType num_points = r_geom.PointsNumber();

        for (IndexType j = 0; j < num_points; ++j) {
            const IndexType node_id = r_geom[j].Id();
            const IndexType cp_col = node_id_to_cp_index.at(node_id);
            (*pNReduced)(gp_idx, cp_col) = r_N[j];
        }

        ++gp_idx;
    }

    return pNReduced;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class RadialBasisFunctionMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

}  // namespace Kratos.