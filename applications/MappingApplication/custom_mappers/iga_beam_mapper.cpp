//    Kratos Multiphysics
//
//  License: BSD License
//           Kratos default license: kratos/license.txt

// System includes
#include <cmath>
#include <unordered_map>

// Project includes
#include "iga_beam_mapper.h"
#include "mappers/mapper_define.h"
#include "utilities/math_utils.h"
#include "utilities/nurbs_utilities/projection_nurbs_geometry_utilities.h"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace>
IgaBeamMapper<TSparseSpace, TDenseSpace>::IgaBeamMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination)
    : mrModelPartOrigin(rModelPartOrigin),
      mrModelPartDestination(rModelPartDestination)
{
}

template<class TSparseSpace, class TDenseSpace>
IgaBeamMapper<TSparseSpace, TDenseSpace>::IgaBeamMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters)
    : IgaBeamMapper(rModelPartOrigin, rModelPartDestination)
{
    KRATOS_TRY

    JsonParameters.ValidateAndAssignDefaults(Parameters(R"({
        "echo_level" : 0,
        "projection_tolerance" : 1e-8,
        "projection_max_iterations" : 50
    })"));

    mProjectionTolerance = JsonParameters["projection_tolerance"].GetDouble();
    mProjectionMaxIterations = JsonParameters["projection_max_iterations"].GetInt();
    KRATOS_ERROR_IF(!std::isfinite(mProjectionTolerance) || mProjectionTolerance <= 0.0)
        << "projection_tolerance must be finite and positive." << std::endl;
    KRATOS_ERROR_IF(mProjectionMaxIterations <= 0)
        << "projection_max_iterations must be positive." << std::endl;

    InitializeInput();
    InitializeReferenceAttachments();

    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::InitializeReferenceAttachments()
{
    KRATOS_TRY
    // Never move the structural nodes to evaluate the reference curve.
    typename NurbsCurveType::PointsArrayType reference_points;
    for (const auto& r_node : *mpReferenceCurve) {
        reference_points.push_back(Kratos::make_intrusive<Node>(
            r_node.Id(), r_node.X0(), r_node.Y0(), r_node.Z0()));
    }
    NurbsCurveType reference_curve = mpReferenceCurve->IsRational()
        ? NurbsCurveType(reference_points, mpReferenceCurve->PolynomialDegree(0),
            mpReferenceCurve->Knots(), mpReferenceCurve->Weights())
        : NurbsCurveType(reference_points, mpReferenceCurve->PolynomialDegree(0),
            mpReferenceCurve->Knots());
    KRATOS_ERROR_IF(reference_curve.PolynomialDegree(0) < 2)
        << "IgaBeamMapper requires a beam curve of degree at least two." << std::endl;

    std::vector<double> spans;
    const auto& r_knots = reference_curve.Knots();
    const std::size_t degree = reference_curve.PolynomialDegree(0);
    for (std::size_t i = degree - 1; i <= r_knots.size() - degree; ++i) {
        if (spans.empty() || r_knots[i] > spans.back()) spans.push_back(r_knots[i]);
    }
    KRATOS_ERROR_IF(spans.size() < 2)
        << "IgaBeamMapper: minimum-distance projection failed; empty curve parameter domain." << std::endl;
    const double start = spans.front();
    const double end = spans.back();
    const std::size_t samples_per_span = 2 * reference_curve.PolynomialDegree(0) + 1;
    std::vector<ReferenceAttachment> attachments;
    attachments.reserve(mrModelPartDestination.NumberOfNodes());

    for (auto& r_node : mrModelPartDestination.Nodes()) {
        ReferenceAttachment attachment;
        attachment.pSurfaceNode = mrModelPartDestination.pGetNode(r_node.Id());
        attachment.ReferencePosition = r_node.GetInitialPosition();
        double best_distance = std::numeric_limits<double>::max();
        double sampled_distance = std::numeric_limits<double>::max();

        // Check stationarity (or the one-sided endpoint condition), independently
        // of the iterative solver's convergence flag. Reject stationary maxima.
        const auto consider_candidate = [&](const double Parameter) {
            if (!std::isfinite(Parameter) || Parameter < start || Parameter > end) return;
            array_1d<double, 3> local = ZeroVector(3);
            local[0] = Parameter;
            std::vector<array_1d<double, 3>> derivatives;
            reference_curve.GlobalSpaceDerivatives(derivatives, local, 2);
            const array_1d<double, 3> offset = derivatives[0] - attachment.ReferencePosition;
            const double tangent_length = norm_2(derivatives[1]);
            if (!std::isfinite(tangent_length) || tangent_length <= std::numeric_limits<double>::epsilon()) return;
            const double gradient = inner_prod(offset, derivatives[1]);
            const bool endpoint_minimum = (Parameter == start && gradient >= 0.0) ||
                (Parameter == end && gradient <= 0.0);
            if (!endpoint_minimum && (std::abs(gradient) / tangent_length > mProjectionTolerance ||
                inner_prod(derivatives[1], derivatives[1]) + inner_prod(offset, derivatives[2]) < 0.0)) return;
            const double distance = norm_2(offset);
            if (std::isfinite(distance) && distance < best_distance) {
                best_distance = distance;
                attachment.Parameter = Parameter;
                attachment.CenterlinePosition = derivatives[0];
            }
        };

        consider_candidate(start);
        consider_candidate(end);
        for (std::size_t span = 0; span + 1 < spans.size(); ++span) {
            if (spans[span + 1] <= spans[span]) continue;
            for (std::size_t sample = 0; sample <= samples_per_span; ++sample) {
                array_1d<double, 3> local = ZeroVector(3);
                local[0] = spans[span] + (spans[span + 1] - spans[span]) * sample / samples_per_span;
                array_1d<double, 3> projected = ZeroVector(3);
                reference_curve.GlobalCoordinates(projected, local);
                sampled_distance = std::min(sampled_distance, norm_2(projected - attachment.ReferencePosition));
                consider_candidate(local[0]);
                if (ProjectionNurbsGeometryUtilities::NewtonRaphsonCurve(
                    local, attachment.ReferencePosition, projected, reference_curve,
                    mProjectionMaxIterations, mProjectionTolerance)) {
                    consider_candidate(local[0]);
                }
            }
        }
        KRATOS_ERROR_IF(best_distance == std::numeric_limits<double>::max() ||
            best_distance > sampled_distance + mProjectionTolerance)
            << "IgaBeamMapper: minimum-distance projection failed for surface node "
            << r_node.Id() << "." << std::endl;

        array_1d<double, 3> local = ZeroVector(3);
        local[0] = attachment.Parameter;
        std::vector<std::size_t> control_point_ids;
        DenseVector<Matrix> derivatives;
        reference_curve.ShapeFunctionsValuesAndCPIndices(
            local, control_point_ids, attachment.ShapeFunctions, 1, &derivatives);
        attachment.ShapeFunctionDerivatives = derivatives[0];
        for (const auto id : control_point_ids) {
            attachment.ControlPoints.push_back(mrModelPartOrigin.pGetNode(id));
        }

        // Select an element with the same active support, then use its section
        // properties. No IgaApplication dependency is needed.
        Element* p_source_element = nullptr;
        double closest_parameter_distance = std::numeric_limits<double>::max();
        for (auto& r_element : mrModelPartOrigin.Elements()) {
            const auto& r_geometry = r_element.GetGeometry();
            if (r_geometry.size() != control_point_ids.size()) continue;
            bool same_support = true;
            for (std::size_t i = 0; i < control_point_ids.size(); ++i) {
                same_support = same_support && r_geometry[i].Id() == control_point_ids[i];
            }
            if (!same_support) continue;
            const double distance = std::abs(
                r_geometry.IntegrationPoints()[0].X() - attachment.Parameter);
            if (distance < closest_parameter_distance) {
                closest_parameter_distance = distance;
                p_source_element = &r_element;
            }
        }
        KRATOS_ERROR_IF_NOT(p_source_element)
            << "IgaBeamMapper: no beam element covers the projected location of surface node "
            << r_node.Id() << "." << std::endl;
        typename NurbsCurveType::IntegrationPointsArrayType integration_points;
        integration_points.emplace_back(attachment.Parameter, 0.0, 0.0, 1.0);
        typename NurbsCurveType::GeometriesArrayType quadrature_points;
        auto integration_info = reference_curve.GetDefaultIntegrationInfo();
        reference_curve.CreateQuadraturePointGeometries(
            quadrature_points, 3, integration_points, integration_info);
        auto p_frame_element = p_source_element->Create(
            0, quadrature_points(0), p_source_element->pGetProperties());
        std::vector<Matrix> frames;
        p_frame_element->CalculateOnIntegrationPoints(
            LOCAL_AXES_MATRIX, frames, mrModelPartOrigin.GetProcessInfo());
        KRATOS_ERROR_IF(frames.size() != 1 || frames[0].size1() != 3 || frames[0].size2() != 3)
            << "IgaBeamMapper requires reference frame rows T, N, V from LOCAL_AXES_MATRIX." << std::endl;
        attachment.ReferenceFrame = frames[0];
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double dot = inner_prod(row(frames[0], i), row(frames[0], j));
                KRATOS_ERROR_IF(!std::isfinite(dot) || std::abs(dot - (i == j ? 1.0 : 0.0)) > 1e-8)
                    << "IgaBeamMapper: invalid reference frame at surface node " << r_node.Id() << "." << std::endl;
            }
        }
        std::vector<array_1d<double, 3>> curve_derivatives;
        reference_curve.GlobalSpaceDerivatives(curve_derivatives, local, 1);
        const array_1d<double, 3> tangent = curve_derivatives[1] / norm_2(curve_derivatives[1]);
        const array_1d<double, 3> frame_tangent = row(frames[0], 0);
        const array_1d<double, 3> frame_normal = row(frames[0], 1);
        const array_1d<double, 3> frame_binormal = row(frames[0], 2);
        const array_1d<double, 3> cross = MathUtils<double>::CrossProduct(frame_tangent, frame_normal);
        KRATOS_ERROR_IF(norm_2(tangent - frame_tangent) > 1e-8 || norm_2(cross - frame_binormal) > 1e-8)
            << "IgaBeamMapper: reference frame must follow the beam tangent and be right-handed." << std::endl;
        const array_1d<double, 3> offset = attachment.ReferencePosition - attachment.CenterlinePosition;
        attachment.NormalOffset = inner_prod(offset, row(frames[0], 1));
        attachment.BinormalOffset = inner_prod(offset, row(frames[0], 2));
        const array_1d<double, 3> residual = offset - attachment.NormalOffset * row(frames[0], 1)
            - attachment.BinormalOffset * row(frames[0], 2);
        KRATOS_ERROR_IF(norm_2(residual) > mProjectionTolerance)
            << "IgaBeamMapper: reference attachment for surface node " << r_node.Id()
            << " has a tangential offset; the projected cross-section cannot reconstruct the point." << std::endl;
        attachments.push_back(std::move(attachment));
    }
    mReferenceAttachments = std::move(attachments);
    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::InitializeInput()
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mrModelPartOrigin.IsDistributed() || mrModelPartDestination.IsDistributed())
        << "IgaBeamMapper currently supports serial model parts only." << std::endl;
    KRATOS_ERROR_IF(mrModelPartOrigin.NumberOfElements() == 0)
        << "IgaBeamMapper requires origin beam elements with a parent NURBS curve." << std::endl;
    KRATOS_ERROR_IF(mrModelPartDestination.NumberOfNodes() == 0)
        << "IgaBeamMapper requires destination surface nodes." << std::endl;

    for (const auto& r_element : mrModelPartOrigin.Elements()) {
        const auto& r_geometry = r_element.GetGeometry();
        KRATOS_ERROR_IF(r_geometry.GetGeometryFamily() !=
            GeometryData::KratosGeometryFamily::Kratos_Quadrature_Geometry)
            << "IgaBeamMapper expects quadrature-point beam elements; element "
            << r_element.Id() << " has an unsupported geometry." << std::endl;

        const auto* p_curve = dynamic_cast<const NurbsCurveType*>(&r_geometry.GetGeometryParent(0));
        KRATOS_ERROR_IF_NOT(p_curve)
            << "IgaBeamMapper requires a NurbsCurveGeometry3D parent for element "
            << r_element.Id() << "." << std::endl;
        KRATOS_ERROR_IF(mpReferenceCurve != nullptr && mpReferenceCurve != p_curve)
            << "IgaBeamMapper currently supports exactly one parent NURBS curve." << std::endl;
        mpReferenceCurve = p_curve;
    }

    for (const auto& r_node : *mpReferenceCurve) {
        KRATOS_ERROR_IF_NOT(mrModelPartOrigin.HasNode(r_node.Id()))
            << "IgaBeamMapper: control point " << r_node.Id()
            << " is missing from the origin model part." << std::endl;
        KRATOS_ERROR_IF(&mrModelPartOrigin.GetNode(r_node.Id()) != &r_node)
            << "IgaBeamMapper: control point " << r_node.Id()
            << " does not match the origin node." << std::endl;
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(DISPLACEMENT))
            << "IgaBeamMapper: control point " << r_node.Id()
            << " is missing historical DISPLACEMENT." << std::endl;
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(ROTATION_X))
            << "IgaBeamMapper: control point " << r_node.Id()
            << " is missing historical ROTATION_X (scalar twist)." << std::endl;
    }

    for (const auto& r_node : mrModelPartDestination.Nodes()) {
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(DISPLACEMENT))
            << "IgaBeamMapper: surface node " << r_node.Id()
            << " is missing historical DISPLACEMENT." << std::endl;
    }

    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
typename IgaBeamMapper<TSparseSpace, TDenseSpace>::MapperUniquePointerType
IgaBeamMapper<TSparseSpace, TDenseSpace>::Clone(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters) const
{
    return Kratos::make_unique<IgaBeamMapper>(
        rModelPartOrigin, rModelPartDestination, JsonParameters);
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::UpdateInterface(
    Kratos::Flags MappingOptions,
    double SearchRadius)
{
    KRATOS_ERROR << "IgaBeamMapper::UpdateInterface is not implemented yet." << std::endl;
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    KRATOS_ERROR << "IgaBeamMapper does not support scalar field mapping." << std::endl;
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    KRATOS_TRY
    KRATOS_ERROR_IF(rOriginVariable != DISPLACEMENT || rDestinationVariable != DISPLACEMENT)
        << "IgaBeamMapper maps DISPLACEMENT to DISPLACEMENT and reads ROTATION_X as scalar twist."
        << std::endl;
    KRATOS_ERROR_IF(MappingOptions != Flags())
        << "IgaBeamMapper currently supports total historical displacement mapping without flags." << std::endl;

    // Compute all results before writing, so an invalid section leaves the
    // destination unchanged. Positions come from X0 + u, never mesh coordinates.
    std::vector<array_1d<double, 3>> displacements;
    displacements.reserve(mReferenceAttachments.size());
    for (const auto& r_attachment : mReferenceAttachments) {
        displacements.push_back(EvaluateAttachment(r_attachment, nullptr));
    }
    for (std::size_t i = 0; i < mReferenceAttachments.size(); ++i) {
        mReferenceAttachments[i].pSurfaceNode->FastGetSolutionStepValue(DISPLACEMENT) = displacements[i];
    }
    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
array_1d<double, 3> IgaBeamMapper<TSparseSpace, TDenseSpace>::EvaluateAttachment(
    const ReferenceAttachment& r_attachment, Matrix* pTangent) const
{
    array_1d<double, 3> centerline = ZeroVector(3);
    array_1d<double, 3> derivative = ZeroVector(3);
    double twist = 0.0;
    for (std::size_t i = 0; i < r_attachment.ControlPoints.size(); ++i) {
        const auto& r_node = *r_attachment.ControlPoints[i];
        const array_1d<double, 3> position = r_node.GetInitialPosition()
            + r_node.FastGetSolutionStepValue(DISPLACEMENT);
        centerline += r_attachment.ShapeFunctions[i] * position;
        derivative += r_attachment.ShapeFunctionDerivatives(i, 0) * position;
        twist += r_attachment.ShapeFunctions[i] * r_node.FastGetSolutionStepValue(ROTATION_X);
    }
    const double length = norm_2(derivative);
    KRATOS_ERROR_IF(!std::isfinite(length) || length <= std::numeric_limits<double>::epsilon() ||
        !std::isfinite(twist) || !std::isfinite(norm_2(centerline)))
        << "IgaBeamMapper: invalid current beam state at surface node "
        << r_attachment.pSurfaceNode->Id() << "." << std::endl;
    const array_1d<double, 3> tangent = derivative / length;
    const array_1d<double, 3> reference_tangent = row(r_attachment.ReferenceFrame, 0);
    const double cosine = inner_prod(reference_tangent, tangent);
    KRATOS_ERROR_IF(1.0 + cosine <= 1e-12)
        << "IgaBeamMapper: opposite reference and current tangents at surface node "
        << r_attachment.pSurfaceNode->Id() << "; the smallest-rotation frame is singular." << std::endl;
    const array_1d<double, 3> axis = MathUtils<double>::CrossProduct(reference_tangent, tangent);
    const array_1d<double, 3> reference_offset =
        r_attachment.NormalOffset * row(r_attachment.ReferenceFrame, 1)
        + r_attachment.BinormalOffset * row(r_attachment.ReferenceFrame, 2);
    // Same composition as IsogeometricBeamElement: first align T with t,
    // then twist about t. The cached frame already includes reference twist.
    const array_1d<double, 3> transported = cosine * reference_offset
        + MathUtils<double>::CrossProduct(axis, reference_offset)
        + axis * (inner_prod(axis, reference_offset) / (1.0 + cosine));
    const array_1d<double, 3> current_offset = std::cos(twist) * transported
        + std::sin(twist) * MathUtils<double>::CrossProduct(tangent, transported);

    if (pTangent != nullptr) {
        pTangent->resize(3, 4 * r_attachment.ControlPoints.size(), false);
        const double c = std::cos(twist);
        const double s = std::sin(twist);
        const double denominator = 1.0 + cosine;
        const double projection = inner_prod(axis, reference_offset);
        for (std::size_t i = 0; i < r_attachment.ControlPoints.size(); ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                // Differentiate the normalized tangent, minimal rotation and twist.
                array_1d<double, 3> dt = -tangent[j] * tangent;
                dt[j] += 1.0;
                dt *= r_attachment.ShapeFunctionDerivatives(i, 0) / length;
                const double dc = inner_prod(reference_tangent, dt);
                const array_1d<double, 3> da = MathUtils<double>::CrossProduct(reference_tangent, dt);
                const array_1d<double, 3> db = dc * reference_offset
                    + MathUtils<double>::CrossProduct(da, reference_offset)
                    + da * (projection / denominator)
                    + axis * (inner_prod(da, reference_offset) / denominator
                        - projection * dc / (denominator * denominator));
                array_1d<double, 3> dx = c * db + s * (
                    MathUtils<double>::CrossProduct(dt, transported)
                    + MathUtils<double>::CrossProduct(tangent, db));
                dx[j] += r_attachment.ShapeFunctions[i];
                for (std::size_t k = 0; k < 3; ++k) (*pTangent)(k, 4 * i + j) = dx[k];
            }
            const array_1d<double, 3> dx = r_attachment.ShapeFunctions[i] * (
                -s * transported + c * MathUtils<double>::CrossProduct(tangent, transported));
            for (std::size_t k = 0; k < 3; ++k) (*pTangent)(k, 4 * i + 3) = dx[k];
        }
    }
    return centerline + current_offset - r_attachment.ReferencePosition;
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::InverseMap(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    KRATOS_ERROR << "IgaBeamMapper does not support scalar field mapping." << std::endl;
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::InverseMap(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    KRATOS_TRY
    // Resolve StructuralMechanics variables at runtime, without an application dependency.
    const auto& r_load = KratosComponents<Variable<array_1d<double, 3>>>::Get("POINT_LOAD");
    const auto& r_moment = KratosComponents<Variable<array_1d<double, 3>>>::Get("POINT_MOMENT");
    KRATOS_ERROR_IF(rOriginVariable != r_load)
        << "IgaBeamMapper::InverseMap requires POINT_LOAD as origin output variable." << std::endl;
    KRATOS_ERROR_IF(MappingOptions != Flags())
        << "IgaBeamMapper force mapping supports historical nodal forces without flags." << std::endl;
    std::unordered_map<std::size_t, std::size_t> indices;
    std::vector<array_1d<double, 3>> loads(mpReferenceCurve->size(), ZeroVector(3));
    std::vector<double> moments(mpReferenceCurve->size(), 0.0);
    for (std::size_t i = 0; i < mpReferenceCurve->size(); ++i) {
        const auto& r_node = (*mpReferenceCurve)[i];
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(r_load) && r_node.SolutionStepsDataHas(r_moment))
            << "IgaBeamMapper requires historical POINT_LOAD and POINT_MOMENT on control points." << std::endl;
        indices.emplace(r_node.Id(), i);
    }
    Matrix tangent;
    for (const auto& r_attachment : mReferenceAttachments) {
        KRATOS_ERROR_IF_NOT(r_attachment.pSurfaceNode->SolutionStepsDataHas(rDestinationVariable))
            << "IgaBeamMapper: missing historical surface force variable." << std::endl;
        const auto& r_force = r_attachment.pSurfaceNode->FastGetSolutionStepValue(rDestinationVariable);
        KRATOS_ERROR_IF(!std::isfinite(norm_2(r_force)))
            << "IgaBeamMapper: non-finite surface force." << std::endl;
        EvaluateAttachment(r_attachment, &tangent);
        for (std::size_t i = 0; i < r_attachment.ControlPoints.size(); ++i) {
            const auto index = indices.at(r_attachment.ControlPoints[i]->Id());
            for (std::size_t k = 0; k < 3; ++k) {
                for (std::size_t j = 0; j < 3; ++j) loads[index][j] += tangent(k, 4 * i + j) * r_force[k];
                moments[index] += tangent(k, 4 * i + 3) * r_force[k];
            }
        }
    }
    // Commit only after every attachment and input force has been validated.
    for (std::size_t i = 0; i < mpReferenceCurve->size(); ++i) {
        auto& r_node = mrModelPartOrigin.GetNode((*mpReferenceCurve)[i].Id());
        r_node.FastGetSolutionStepValue(r_load) = loads[i];
        auto& r_value = r_node.FastGetSolutionStepValue(r_moment);
        r_value = ZeroVector(3);
        r_value[0] = moments[i];
    }
    KRATOS_CATCH("")
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::Map(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rRotationVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    KRATOS_ERROR_IF(rRotationVariable != ROTATION)
        << "IgaBeamMapper: invalid secondary rotation variable." << std::endl;
    Map(rOriginVariable, rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::InverseMap(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rMomentVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    KRATOS_ERROR_IF((rMomentVariable != KratosComponents<Variable<array_1d<double, 3>>>::Get("POINT_MOMENT")))
        << "IgaBeamMapper: invalid secondary moment variable." << std::endl;
    InverseMap(rOriginVariable, rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
ModelPart& IgaBeamMapper<TSparseSpace, TDenseSpace>::GetInterfaceModelPartOrigin()
{
    return mrModelPartOrigin;
}

template<class TSparseSpace, class TDenseSpace>
ModelPart& IgaBeamMapper<TSparseSpace, TDenseSpace>::GetInterfaceModelPartDestination()
{
    return mrModelPartDestination;
}

template<class TSparseSpace, class TDenseSpace>
std::string IgaBeamMapper<TSparseSpace, TDenseSpace>::Info() const
{
    return "IgaBeamMapper";
}

template<class TSparseSpace, class TDenseSpace>
void IgaBeamMapper<TSparseSpace, TDenseSpace>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

template class IgaBeamMapper<MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType>;

} // namespace Kratos
