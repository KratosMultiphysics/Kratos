//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Andrea Gorgi
//                   Nicoló Antonelli
//
// System includes
// External includes

// Project includes
#include "snake_gap_sbm_process.h"
#include "integration/integration_point_utilities.h"
#include "iga_application_variables.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{
namespace
{
} // unnamed namespace

//-------------------------------

SnakeGapSbmProcess::SnakeGapSbmProcess(
    Model& rModel, Parameters ThisParameters) : 
    SnakeSbmProcess(
        rModel,
        // Reject user-defined lambdas before the base constructor adds their default values.
        (ThisParameters.Has("lambda_inner") || ThisParameters.Has("lambda_outer"))
            ? KRATOS_ERROR << "[SnakeGapSbmProcess]: \"lambda_inner\" and \"lambda_outer\" "
                           << "must not be provided in ProjectParameters. GapSBM assigns these "
                           << "values internally." << std::endl
            : ThisParameters)
{
    
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("gap_element_name")) << "::[SnakeGapSbmProcess]::" 
                    << "Missing \"gap_element_name\" section." << std::endl;
    KRATOS_ERROR_IF_NOT(ThisParameters.Has("gap_interface_condition_name")) << "::[SnakeGapSbmProcess]::" 
                    << "Missing \"gap_interface_condition_name\" section." << std::endl;

    ThisParameters.AddMissingParameters(this->GetDefaultParameters());
    
    mpGapElementsSubModelPart = &(mpIgaModelPart->CreateSubModelPart("GapElements"));
    mpGapInterfaceSubModelPart = &(mpIgaModelPart->CreateSubModelPart("GapInterfaces"));
    mGapElementName = ThisParameters["gap_element_name"].GetString();
    mGapInterfaceConditionName = ThisParameters["gap_interface_condition_name"].GetString();
    mGapSbmType = ThisParameters["gap_sbm_type"].GetString(); 

    if (mGapSbmType != "default" && mGapSbmType != "interpolation" && mGapSbmType != "sbm") {
        KRATOS_ERROR << "::[SnakeGapSbmProcess]::"
                     << "The gap_sbm_type \"" << mGapSbmType << "\" is not supported. Available options are: "
                     << "default, interpolation, sbm." << std::endl;
    }

    if (ThisParameters.Has("gap_approximation_order"))
        mGapApproximationOrder = ThisParameters["gap_approximation_order"].GetInt();
    if (ThisParameters.Has("gap_relative_tolerance_for_subdivisions"))
        mGapRelativeToleranceForSubdivisions = ThisParameters["gap_relative_tolerance_for_subdivisions"].GetDouble();
    if (ThisParameters.Has("number_of_interpolation_levels"))
        mNumberOfInterpolationLevels = ThisParameters["number_of_interpolation_levels"].GetInt();
    if (ThisParameters.Has("polynomial_order"))
        mGapInterpolationOrder = ThisParameters["polynomial_order"].GetVector()[0];

    mLambdaInner = 0.0;
    mLambdaOuter = 1.0;
    mInternalDivisions = ThisParameters["number_internal_divisions"].GetInt();

}

//-------------------------------

void SnakeGapSbmProcess::CreateSbmExtendedGeometries()
{
    CreateSbmExtendedGeometries2D();
}

//-------------------------------

const Parameters SnakeGapSbmProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 1,
        "model_part_name" : "IgaModelPart",
        "lower_point_xyz": [0.0, 0.0, 0.0],
        "upper_point_xyz": [1.0, 1.0, 0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [10, 10],
        "gap_relative_tolerance_for_subdivisions": 0.1,
        "number_of_interpolation_levels": 3
    })");
}

//-------------------------------

const Parameters SnakeGapSbmProcess::GetValidParameters() const
{
    return Parameters(R"(
    {
        "echo_level": 0,
        "lower_point_xyz": [-0.5, -0.5,0.0],
        "upper_point_xyz": [0.5,0.5,0.0],
        "polynomial_order" : [2, 2],
        "number_of_knot_spans" : [7, 7],
        "number_of_inner_loops": 0,
        "number_initial_points_if_importing_nurbs": 1000,
        "number_internal_divisions": 0,
        "gap_relative_tolerance_for_subdivisions": 0.1,
        "number_of_interpolation_levels": 3,
        "gap_sbm_type": "default",
        "skin_model_part_outer_initial_name": "initial_skin_model_part_out",    
        "skin_model_part_inner_initial_name": "initial_skin_model_part_in",           
        "skin_model_part_name": "skin_model_part",
        "gap_element_name": "CutSbmSolidElement",
        "gap_interface_condition_name": "CutSbmSolidInterfaceCondition"
    })");
}



//-------------------------------

void SnakeGapSbmProcess::ConditionConfigure::CalculateBoundingBox(
            const PointerType& rObject,
            PointType& rLowPoint,
            PointType& rHighPoint)
{
            const auto& r_geometry = rObject->GetGeometry();
            if (r_geometry.size() == 0) {
                rLowPoint = PointType(0.0, 0.0, 0.0);
                rHighPoint = rLowPoint;
                return;
            }
            rLowPoint = r_geometry[0];
            rHighPoint = r_geometry[0];
            for (IndexType i = 1; i < r_geometry.size(); ++i) {
                const auto& r_point = r_geometry[i];
                for (std::size_t d = 0; d < Dimension; ++d) {
                    rLowPoint[d] = std::min(rLowPoint[d], r_point[d]);
                    rHighPoint[d] = std::max(rHighPoint[d], r_point[d]);
                }
            }
        }

//-------------------------------

void SnakeGapSbmProcess::ConditionConfigure::CalculateBoundingBox(
            const PointerType& rObject,
            PointType& rLowPoint,
            PointType& rHighPoint,
            const double Radius)
{
            CalculateBoundingBox(rObject, rLowPoint, rHighPoint);
            for (std::size_t d = 0; d < Dimension; ++d) {
                rLowPoint[d] -= Radius;
                rHighPoint[d] += Radius;
            }
        }

//-------------------------------

void SnakeGapSbmProcess::ConditionConfigure::CalculateCenter(
            const PointerType& rObject,
            PointType& rCentralPoint)
{
            rCentralPoint = rObject->GetGeometry().Center();
        }

//-------------------------------

bool SnakeGapSbmProcess::ConditionConfigure::Intersection(
            const PointerType& rObj_1,
            const PointerType& rObj_2)
{
            PointType low_1, high_1, low_2, high_2;
            CalculateBoundingBox(rObj_1, low_1, high_1);
            CalculateBoundingBox(rObj_2, low_2, high_2);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high_1[d] < low_2[d] - Epsilon || low_1[d] > high_2[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

//-------------------------------

bool SnakeGapSbmProcess::ConditionConfigure::Intersection(
            const PointerType& rObj_1,
            const PointerType& rObj_2,
            const double Radius)
{
            PointType low_1, high_1, low_2, high_2;
            CalculateBoundingBox(rObj_1, low_1, high_1, Radius);
            CalculateBoundingBox(rObj_2, low_2, high_2, Radius);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high_1[d] < low_2[d] - Epsilon || low_1[d] > high_2[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

//-------------------------------

bool SnakeGapSbmProcess::ConditionConfigure::IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint)
{
            PointType low, high;
            CalculateBoundingBox(rObject, low, high);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high[d] < rLowPoint[d] - Epsilon || low[d] > rHighPoint[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

//-------------------------------

bool SnakeGapSbmProcess::ConditionConfigure::IntersectionBox(
            const PointerType& rObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint,
            const double Radius)
{
            PointType low, high;
            CalculateBoundingBox(rObject, low, high, Radius);
            for (std::size_t d = 0; d < Dimension; ++d) {
                if (high[d] < rLowPoint[d] - Epsilon || low[d] > rHighPoint[d] + Epsilon) {
                    return false;
                }
            }
            return true;
        }

//-------------------------------

void SnakeGapSbmProcess::ConditionConfigure::Distance(
            const PointerType& rObj_1,
            const PointerType& rObj_2,
            double& distance)
{
            PointType c1, c2;
            CalculateCenter(rObj_1, c1);
            CalculateCenter(rObj_2, c2);
            double pwdDistance = 0.0;
            for (std::size_t d = 0; d < Dimension; ++d) {
                const double delta = c1[d] - c2[d];
                pwdDistance += delta * delta;
            }
            distance = std::sqrt(pwdDistance);
        }

//-------------------------------

std::size_t SnakeGapSbmProcess::FindNnzIndex(const SparseMatrixType& A, std::size_t i, std::size_t j)
{
        if (i >= A.size1() || j >= A.size2()) return static_cast<std::size_t>(-1);

        const auto& rp = A.index1_data();   // row_ptr, size = size1+1
        const auto& ci = A.index2_data();   // col_ind, size = nnz

        if (rp.size() < A.size1() + 1) return static_cast<std::size_t>(-1);
        std::size_t b = rp[i];
        std::size_t e = rp[i + 1];

        // guard rails against corrupted/unfinalized structure
        const std::size_t nnz = static_cast<std::size_t>(ci.size());
        if (b > e) return static_cast<std::size_t>(-1);
        if (b > nnz) b = nnz;
        if (e > nnz) e = nnz;

        const auto first = ci.begin() + b;
        const auto last  = ci.begin() + e;
        const auto it = std::lower_bound(first, last, j);
        if (it != last && *it == j) return static_cast<std::size_t>(it - ci.begin());
        return static_cast<std::size_t>(-1);
    }

//-------------------------------

std::size_t SnakeGapSbmProcess::FlattenSpanColumnIndex(
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const std::size_t SpanIndexY,
        const std::size_t SpanIndexZ)
{
        return SpanIndexY + rSkinBinsPerSpan.NumberOfSpansY * SpanIndexZ;
    }

//-------------------------------

std::size_t SnakeGapSbmProcess::FindSpanNnzIndex(
        const KnotSpanSkinBinsCSR& rSkinBinsPerSpan,
        const std::size_t SpanIndexX,
        const std::size_t SpanIndexY,
        const std::size_t SpanIndexZ)
{
        return FindNnzIndex(
            rSkinBinsPerSpan.Occupancy,
            SpanIndexX,
            FlattenSpanColumnIndex(rSkinBinsPerSpan, SpanIndexY, SpanIndexZ));
    }

//-------------------------------

void SnakeGapSbmProcess::CommitPayload(std::vector<std::vector<IndexType>>& tmp,
                                    KnotSpanIdsCSR& out)
{
        out.nnz_off.resize(tmp.size());
        out.nnz_len.resize(tmp.size());

        for (std::size_t k = 0; k < tmp.size(); ++k) {
            auto& v = tmp[k];
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());

            out.nnz_off[k] = static_cast<std::size_t>(out.pool.size());
            out.nnz_len[k] = static_cast<std::size_t>(v.size());
            out.pool.insert(out.pool.end(), v.begin(), v.end());
        }
    }

//-------------------------------

SnakeGapSbmProcess::IdsView SnakeGapSbmProcess::CellIdsByK(const KnotSpanIdsCSR& S, const std::size_t k)
{
        if (k >= S.nnz_off.size() || k >= S.nnz_len.size()) {
            KRATOS_WARNING_ONCE("SnakeGapSbmProcess")
                << "CellIdsByK received index " << k << " but nnz_off has size "
                << S.nnz_off.size() << " and nnz_len has size " << S.nnz_len.size() << ". "
                << "Returning empty view to avoid accessing invalid data." << std::endl;
            return {};
        }

        const std::size_t offset = S.nnz_off[k];
        const std::size_t length = S.nnz_len[k];
        const std::size_t pool_size = S.pool.size();

        if (offset >= pool_size) {
            KRATOS_WARNING_ONCE("SnakeGapSbmProcess")
                << "CellIdsByK computed offset " << offset << " beyond pool size "
                << pool_size << ". Returning empty view to avoid invalid access." << std::endl;
            return {};
        }

        const std::size_t max_length = pool_size - offset;
        const std::size_t safe_length = std::min(length, max_length);
        if (safe_length != length) {
            KRATOS_WARNING_ONCE("SnakeGapSbmProcess")
                << "CellIdsByK clipped entry length from " << length << " to "
                << safe_length << " to stay within pool bounds (" << pool_size << ")." << std::endl;
        }

        return IdsView{ S.pool.data() + offset, safe_length };
    }

//-------------------------------

SnakeGapSbmProcess::IdsView SnakeGapSbmProcess::CellIds(const KnotSpanIdsCSR& S, const std::size_t i, const std::size_t j)
{
        const std::size_t k = FindNnzIndex(S.Occupancy, i, j);
        if (k == static_cast<std::size_t>(-1)) return {};
        return CellIdsByK(S, k);
    }

//-------------------------------

SnakeGapSbmProcess::NurbsCurveGeometryType::Pointer SnakeGapSbmProcess::CreateBrepCurve(
        const Node::Pointer pFirstBrepPoint,
        const Node::Pointer pSecondBrepPoint,
        const Vector& rActiveRangeKnotVector)
{
            // Create the data for the trimming curves
            PointerVector<Node> control_points;
            control_points.push_back(pFirstBrepPoint);
            control_points.push_back(pSecondBrepPoint);
            const int polynomial_degree = 1;
            Vector knot_vector = ZeroVector(4) ;
            knot_vector[0] = rActiveRangeKnotVector[0] ;
            knot_vector[1] = rActiveRangeKnotVector[0] ;
            knot_vector[2] = rActiveRangeKnotVector[1] ;
            knot_vector[3] = rActiveRangeKnotVector[1] ;
            // Create the trimming curves
            typename NurbsCurveGeometry<3, PointerVector<Node>>::Pointer p_trimming_curve(
                new NurbsCurveGeometry<3, PointerVector<Node>>(
                    control_points,
                    polynomial_degree,
                    knot_vector));
            return p_trimming_curve;
        }

}  // namespace Kratos.
