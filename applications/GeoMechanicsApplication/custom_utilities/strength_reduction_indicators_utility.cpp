// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "custom_utilities/strength_reduction_indicators_utility.h"

#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <tuple>

namespace
{

using namespace Kratos;

constexpr auto number_of_normal_components = std::size_t{3};

double DisplacementNorm(const Node& rNode)
{
    return norm_2(rNode.FastGetSolutionStepValue(DISPLACEMENT));
}

double IncrementalDisplacementNorm(const Node& rNode)
{
    return norm_2(rNode.FastGetSolutionStepValue(DISPLACEMENT, 0) -
                  rNode.FastGetSolutionStepValue(DISPLACEMENT, 1));
}

template <typename TNodalNormFunction>
double CalculateMaxOfNodalNorms(const ModelPart& rModelPart, TNodalNormFunction NodalNorm)
{
    const auto local_max = block_for_each<MaxReduction<double>>(
        rModelPart.GetCommunicator().LocalMesh().Nodes(), [&NodalNorm](const Node& rNode) {
        return NodalNorm(rNode);
    });
    return rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_max);
}

template <typename TNodalNormFunction>
double CalculateSumOfNodalNorms(const ModelPart& rModelPart, TNodalNormFunction NodalNorm)
{
    const auto local_sum = block_for_each<SumReduction<double>>(
        rModelPart.GetCommunicator().LocalMesh().Nodes(), [&NodalNorm](const Node& rNode) {
        return NodalNorm(rNode);
    });
    return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_sum);
}

bool IsContinuumElement(const Element& rElement, std::size_t ContinuumDimension)
{
    return rElement.IsActive() && rElement.GetGeometry().LocalSpaceDimension() == ContinuumDimension;
}

std::size_t GetContinuumDimension(const ModelPart& rModelPart)
{
    // The continuum elements are the ones with the highest local dimension (this excludes
    // interface, line and beam elements)
    const auto local_dimension = block_for_each<MaxReduction<int>>(
        rModelPart.GetCommunicator().LocalMesh().Elements(), [](const Element& rElement) {
        return static_cast<int>(rElement.GetGeometry().LocalSpaceDimension());
    });
    return static_cast<std::size_t>(rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(local_dimension));
}

Vector CalculateNodalWeights(const Geometry<Node>& rGeometry)
{
    // Integral of each shape function over the element. Where one of them is not positive (e.g.
    // the corner nodes of 6-noded triangles and 10-noded tetrahedra), the element measure is
    // distributed equally among its nodes instead.
    const auto  integration_method   = rGeometry.GetDefaultIntegrationMethod();
    const auto& r_integration_points = rGeometry.IntegrationPoints(integration_method);
    const auto& r_shape_functions    = rGeometry.ShapeFunctionsValues(integration_method);
    Vector      determinants_of_jacobian;
    rGeometry.DeterminantOfJacobian(determinants_of_jacobian, integration_method);

    auto result = Vector{ZeroVector{rGeometry.PointsNumber()}};
    for (std::size_t g = 0; g < r_integration_points.size(); ++g) {
        result += row(r_shape_functions, g) * r_integration_points[g].Weight() * determinants_of_jacobian[g];
    }

    if (std::ranges::any_of(result, [](double Weight) { return Weight <= 0.0; })) {
        const auto element_measure = sum(result);
        std::ranges::fill(result, element_measure / static_cast<double>(result.size()));
    }
    return result;
}

template <typename TNodalNormFunction>
double CalculateAreaWeightedMeanOfNodalNorms(const ModelPart& rModelPart, TNodalNormFunction NodalNorm)
{
    // sum_i A_i ||u_i|| is assembled element by element, so the nodal areas A_i are never stored
    using WeightedSumReduction = CombinedReduction<SumReduction<double>, SumReduction<double>>;

    const auto continuum_dimension = GetContinuumDimension(rModelPart);
    const auto [local_weighted_sum, local_measure] = block_for_each<WeightedSumReduction>(
        rModelPart.GetCommunicator().LocalMesh().Elements(),
        [continuum_dimension, &NodalNorm](const Element& rElement) {
        if (!IsContinuumElement(rElement, continuum_dimension)) return std::make_tuple(0.0, 0.0);

        const auto& r_geometry    = rElement.GetGeometry();
        const auto  nodal_weights = CalculateNodalWeights(r_geometry);
        auto        weighted_sum  = 0.0;
        for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
            weighted_sum += nodal_weights[i] * NodalNorm(r_geometry[i]);
        }
        return std::make_tuple(weighted_sum, sum(nodal_weights));
    });

    const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();
    const auto  measure             = r_data_communicator.SumAll(local_measure);
    KRATOS_ERROR_IF(measure <= 0.0) << "StrengthReductionIndicatorsUtility: the total area/volume "
                                       "of the continuum elements is not positive."
                                    << std::endl;
    return r_data_communicator.SumAll(local_weighted_sum) / measure;
}

bool IsContinuumVoigtSize(std::size_t StressVectorSize)
{
    // Plane strain [xx, yy, zz, xy] and 3D [xx, yy, zz, xy, yz, xz]
    return StressVectorSize == VOIGT_SIZE_2D_PLANE_STRAIN || StressVectorSize == VOIGT_SIZE_3D;
}

Matrix CalculateInverseElasticMatrix(const Element& rElement, std::size_t StressVectorSize)
{
    const auto& r_properties = rElement.GetProperties();
    KRATOS_ERROR_IF_NOT(r_properties.Has(YOUNG_MODULUS) && r_properties.Has(POISSON_RATIO))
        << "StrengthReductionIndicatorsUtility: the properties of element " << rElement.Id()
        << " must contain YOUNG_MODULUS and POISSON_RATIO to compute the elastic energy." << std::endl;

    const auto [youngs_modulus, poissons_ratio] =
        ConstitutiveLawUtilities::GetOrCalculateElasticProperties(r_properties);
    const auto elastic_matrix = ConstitutiveLawUtilities::MakeContinuumElasticConstitutiveTensor(
        youngs_modulus, poissons_ratio, StressVectorSize, number_of_normal_components);

    auto   result      = Matrix{};
    double determinant = 0.0;
    MathUtils<double>::InvertMatrix(elastic_matrix, result, determinant);
    return result;
}

} // namespace

namespace Kratos
{

StrengthReductionIndicatorsUtility::StrengthReductionIndicatorsUtility(ModelPart& rModelPart, Parameters Settings)
    : mrModelPart(rModelPart)
{
    const auto default_settings = Parameters(R"(
    {
        "output_file_name": "srf_indicators.txt"
    })");
    Settings.ValidateAndAssignDefaults(default_settings);
    mOutputFileName = Settings["output_file_name"].GetString();

    KRATOS_ERROR_IF(mrModelPart.GetBufferSize() < 2)
        << "StrengthReductionIndicatorsUtility: a buffer size of at least 2 is required to compute "
           "the incremental displacements."
        << std::endl;

    if (!IsWritingRank()) return;

    mOutputFile.open(mOutputFileName, std::ios::out | std::ios::trunc);
    KRATOS_ERROR_IF_NOT(mOutputFile.is_open())
        << "StrengthReductionIndicatorsUtility: could not open " << mOutputFileName << std::endl;
    mOutputFile << "step\tSRF\tconverged\tNR_iterations\tresidual_norm\tU_max\tsum_U\tU_mean_"
                   "area\tdU_max\tsum_dU\tdU_mean_area\tW_elastic\n";
    mOutputFile.flush();
}

double StrengthReductionIndicatorsUtility::CalculateMaxDisplacement() const
{
    return CalculateMaxOfNodalNorms(mrModelPart, DisplacementNorm);
}

double StrengthReductionIndicatorsUtility::CalculateSumOfDisplacementNorms() const
{
    return CalculateSumOfNodalNorms(mrModelPart, DisplacementNorm);
}

double StrengthReductionIndicatorsUtility::CalculateAreaWeightedMeanDisplacement() const
{
    return CalculateAreaWeightedMeanOfNodalNorms(mrModelPart, DisplacementNorm);
}

double StrengthReductionIndicatorsUtility::CalculateMaxIncrementalDisplacement() const
{
    return CalculateMaxOfNodalNorms(mrModelPart, IncrementalDisplacementNorm);
}

double StrengthReductionIndicatorsUtility::CalculateSumOfIncrementalDisplacementNorms() const
{
    return CalculateSumOfNodalNorms(mrModelPart, IncrementalDisplacementNorm);
}

double StrengthReductionIndicatorsUtility::CalculateAreaWeightedMeanIncrementalDisplacement() const
{
    return CalculateAreaWeightedMeanOfNodalNorms(mrModelPart, IncrementalDisplacementNorm);
}

double StrengthReductionIndicatorsUtility::CalculateTotalElasticEnergy() const
{
    const auto  continuum_dimension = GetContinuumDimension(mrModelPart);
    const auto& r_process_info      = mrModelPart.GetProcessInfo();

    const auto local_energy = block_for_each<SumReduction<double>>(
        mrModelPart.GetCommunicator().LocalMesh().Elements(),
        [continuum_dimension, &r_process_info](Element& rElement) {
        if (!IsContinuumElement(rElement, continuum_dimension)) return 0.0;

        std::vector<Vector> stress_vectors;
        rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress_vectors, r_process_info);
        if (stress_vectors.empty() || !IsContinuumVoigtSize(stress_vectors.front().size())) return 0.0;

        const auto& r_geometry           = rElement.GetGeometry();
        const auto  integration_method   = rElement.GetIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(integration_method);
        KRATOS_ERROR_IF_NOT(r_integration_points.size() == stress_vectors.size())
            << "StrengthReductionIndicatorsUtility: element " << rElement.Id() << " returns "
            << stress_vectors.size() << " stress vectors for " << r_integration_points.size()
            << " integration points." << std::endl;
        Vector determinants_of_jacobian;
        r_geometry.DeterminantOfJacobian(determinants_of_jacobian, integration_method);

        const auto inverse_elastic_matrix =
            CalculateInverseElasticMatrix(rElement, stress_vectors.front().size());

        auto element_energy = 0.0;
        for (std::size_t g = 0; g < r_integration_points.size(); ++g) {
            const Vector elastic_strain = prod(inverse_elastic_matrix, stress_vectors[g]);
            element_energy += 0.5 * inner_prod(stress_vectors[g], elastic_strain) *
                              r_integration_points[g].Weight() * determinants_of_jacobian[g];
        }
        return element_energy;
    });

    return mrModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_energy);
}

void StrengthReductionIndicatorsUtility::RecordIncrement(std::size_t Step, double StrengthReductionFactor, bool IsConverged)
{
    const auto& r_process_info = mrModelPart.GetProcessInfo();
    const auto  nr_iterations  = r_process_info[NL_ITERATION_NUMBER];
    const auto  residual_norm  = r_process_info.Has(RESIDUAL_NORM)
                                     ? r_process_info[RESIDUAL_NORM]
                                     : std::numeric_limits<double>::quiet_NaN();

    // All ranks take part in the reductions; only the writing rank writes
    const auto max_displacement               = CalculateMaxDisplacement();
    const auto sum_of_norms                   = CalculateSumOfDisplacementNorms();
    const auto mean_displacement              = CalculateAreaWeightedMeanDisplacement();
    const auto max_incremental_displacement   = CalculateMaxIncrementalDisplacement();
    const auto sum_of_incremental_norms       = CalculateSumOfIncrementalDisplacementNorms();
    const auto mean_incremental_displacement  = CalculateAreaWeightedMeanIncrementalDisplacement();
    const auto elastic_energy                 = CalculateTotalElasticEnergy();

    if (!IsWritingRank()) return;

    mOutputFile << Step << '\t' << std::scientific << std::setprecision(10)
                << StrengthReductionFactor << '\t' << (IsConverged ? 1 : 0) << '\t' << nr_iterations
                << '\t' << residual_norm << '\t' << max_displacement << '\t' << sum_of_norms << '\t'
                << mean_displacement << '\t' << max_incremental_displacement << '\t'
                << sum_of_incremental_norms << '\t' << mean_incremental_displacement << '\t'
                << elastic_energy << '\n';
    mOutputFile.flush();
}

const std::string& StrengthReductionIndicatorsUtility::GetOutputFileName() const
{
    return mOutputFileName;
}

bool StrengthReductionIndicatorsUtility::IsWritingRank() const
{
    return mrModelPart.GetCommunicator().GetDataCommunicator().Rank() == 0;
}

} // namespace Kratos
