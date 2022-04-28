//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/fluid_adjoint_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"

// Include base h
#include "domain_integrated_3d_vector_magnitude_square_p_mean_response_function.h"

namespace Kratos
{

template<>
template<>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<2>::CalculateDomainSizeDerivative(
    const Element& rElement,
    const IndexType DerivativeNodeIndex,
    const IndexType DerivativeDirectionIndex) const
{
    KRATOS_TRY

    KRATOS_ERROR << "2D element domain size derivative is not implemented.\n";

    KRATOS_CATCH("");
}

template<>
template<>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<3>::CalculateDomainSizeDerivative(
    const Element& rElement,
    const IndexType DerivativeNodeIndex,
    const IndexType DerivativeDirectionIndex) const
{
    KRATOS_TRY

    KRATOS_ERROR << "3D element domain size derivative is not implemented.\n";

    KRATOS_CATCH("");
}

template<>
template<>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<2>::CalculateDomainSizeDerivative(
    const Condition& rCondition,
    const IndexType DerivativeNodeIndex,
    const IndexType DerivativeDirectionIndex) const
{
    KRATOS_TRY

    KRATOS_ERROR << "2D condition domain size derivative is not implemented.\n";

    KRATOS_CATCH("");
}

template<>
template<>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<3>::CalculateDomainSizeDerivative(
    const Condition& rCondition,
    const IndexType DerivativeNodeIndex,
    const IndexType DerivativeDirectionIndex) const
{
    KRATOS_TRY

    const auto& r_geometry = rCondition.GetGeometry();
    const IndexType number_of_nodes = r_geometry.PointsNumber();

    if (number_of_nodes == 3) {
        return FluidAdjointUtilities<3>::CalculateTriangleAreaDerivative(r_geometry, DerivativeNodeIndex, DerivativeDirectionIndex);
    } else {
        KRATOS_ERROR << "3D condition domain size derivative is only defined "
                        "for triangle conditions only. [ requested geometry "
                        "number of nodes = "
                     << number_of_nodes << " ].\n";
        return 0.0;
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim>
template<class TEntityType>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateResponseGradientContribution(
    const TEntityType& rEntity,
    const Matrix& rResidualGradient,
    Vector& rResponseGradient,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1());
    }

    rResponseGradient.clear();

    const double current_time = rProcessInfo[TIME];

    if (mStartTime <= current_time) {
        if (rEntity.Is(*mpDomainFlag)) {
            const int k = OpenMPUtils::ThisThread();
            auto& r_shape_function = mShapeFunctions[k];

            const auto& r_geometry = rEntity.GetGeometry();
            const IndexType number_of_nodes = r_geometry.PointsNumber();
            const IndexType block_size = rResidualGradient.size1() / number_of_nodes;
            const IndexType skip_size = block_size - TDim;

            const double domain_size = r_geometry.DomainSize();

            this->CalculateGeometryData(r_geometry, r_shape_function);
            const Vector& rN = row(r_shape_function, 0);

            array_1d<double, 3> phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, rN, std::tie(phi, *mpVariable));

            const double phi_norm_square_p_1 = std::pow(inner_prod(phi, phi), mPower - 1);

            IndexType local_index = mDofPosition;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                for (IndexType k = 0; k < TDim; ++k) {
                    rResponseGradient[local_index++] =  domain_size * 2.0 * mPower * phi_norm_square_p_1 * rN[c] * phi[k] / mIntegrationDomainSize;
                }
                local_index += skip_size;
            }
        }
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim>
template<class TEntityType>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateResponsePartialSensitivity(
        const TEntityType& rEntity,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
        rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
    }

    rSensitivityGradient.clear();

    const double current_time = rProcessInfo[TIME];

    if (mStartTime <= current_time) {
        if (rEntity.Is(*mpDomainFlag)) {
            const int k = OpenMPUtils::ThisThread();
            auto& r_shape_function = mShapeFunctions[k];

            const auto& r_geometry = rEntity.GetGeometry();
            const IndexType number_of_nodes = r_geometry.PointsNumber();

            this->CalculateGeometryData(r_geometry, r_shape_function);
            const Vector& N = row(r_shape_function, 0);

            array_1d<double, 3> phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(phi, *mpVariable));

            const double inner_phi_power = std::pow(inner_prod(phi, phi), mPower);

            IndexType local_index = 0;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                for (IndexType k = 0; k < TDim; ++k) {
                    rSensitivityGradient[local_index++] =
                        CalculateDomainSizeDerivative(rEntity, c, k) *
                        (inner_phi_power - mDomainIntegratedSquareMean) / mIntegrationDomainSize;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim>
DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction(
    Parameters Settings,
    ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
    KRATOS_TRY;

    Parameters default_settings(R"(
    {
        "model_part_name"          : "",
        "variable_name"            : "PLEASE_SPECIFY_VARIABLE_NAME",
        "dof_position"             : 0,
        "magnitude_square_to_power": 1,
        "flag_to_be_used"          : "STRUCTURE",
        "entities_to_consider"     : ["conditions"],
        "start_time"               : 0.0
    })");

    Settings.ValidateAndAssignDefaults(default_settings);

    mIntegrationDomainModelPartName = Settings["model_part_name"].GetString();
    mpVariable = &KratosComponents<Variable<array_1d<double, 3>>>::Get(Settings["variable_name"].GetString());
    mpDomainFlag = &KratosComponents<Flags>::Get(Settings["flag_to_be_used"].GetString());
    mPower = Settings["magnitude_square_to_power"].GetInt();
    mDofPosition = Settings["dof_position"].GetInt();
    mStartTime = Settings["start_time"].GetDouble();

    const auto& r_entity_types = Settings["entities_to_consider"].GetStringArray();
    for (const auto& r_entity_type : r_entity_types) {
        if (r_entity_type == "elements") {
            mIsElementsConsidered = true;
        } else if (r_entity_type == "conditions") {
            mIsConditionsConsidered = true;
        } else {
            KRATOS_ERROR
                << "Unsupported entity type provided in "
                   "\"entities_to_consider\". Only supports \"elements\" and "
                   "\"conditions\" [ \"entities_to_consider\" = "
                << r_entity_type << " ].\n";
        }
    }

    KRATOS_ERROR_IF(mPower <= 0) << "\"magnitude_square_to_power\" should be greater than zero.";

    KRATOS_CATCH("");
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::Initialize()
{
    KRATOS_TRY;

    const int num_threads = ParallelUtilities::GetNumThreads();
    mShapeFunctions.resize(num_threads);

    VariableUtils().SetFlag(*mpDomainFlag, false, mrModelPart.Elements());
    VariableUtils().SetFlag(*mpDomainFlag, false, mrModelPart.Conditions());

    if (mIntegrationDomainModelPartName != "") {
        auto& r_integration_domain_model_part = mrModelPart.GetSubModelPart(mIntegrationDomainModelPartName);
        if (mIsElementsConsidered) {
            VariableUtils().SetFlag(*mpDomainFlag, true, r_integration_domain_model_part.Elements());
        }

        if (mIsConditionsConsidered) {
            VariableUtils().SetFlag(*mpDomainFlag, true, r_integration_domain_model_part.Conditions());
        }
    }

    mIntegrationDomainSize = CalculateIntegrationDomainSize();

    KRATOS_CATCH("");
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::InitializeSolutionStep()
{
    KRATOS_TRY

    mDomainIntegratedSquareMean = CalculateValue(mrModelPart);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
{
    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
{
    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
{
    CalculateResponseGradientContribution(rAdjointElement, rResidualGradient, rResponseGradient, rProcessInfo);
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
{
    CalculateResponseGradientContribution(rAdjointCondition, rResidualGradient, rResponseGradient, rProcessInfo);
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
{
    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo)
{
    if (rResponseGradient.size() != rResidualGradient.size1()) {
        rResponseGradient.resize(rResidualGradient.size1(), false);
    }

    rResponseGradient.clear();
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == SHAPE_SENSITIVITY) {
        CalculateResponsePartialSensitivity(rAdjointElement, rSensitivityMatrix, rSensitivityGradient, rProcessInfo);
    } else {
        KRATOS_ERROR << "Unsupported partial sensitivity requested. [ rVariable = " << rVariable.Name() << " ].\n";
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY;

    if (rVariable == SHAPE_SENSITIVITY) {
        CalculateResponsePartialSensitivity(rAdjointCondition, rSensitivityMatrix, rSensitivityGradient, rProcessInfo);
    } else {
        KRATOS_ERROR << "Unsupported partial sensitivity requested. [ rVariable = " << rVariable.Name() << " ].\n";
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateValue(ModelPart&)
{
    KRATOS_TRY

    const auto& r_model_part = mrModelPart.GetSubModelPart(mIntegrationDomainModelPartName);
    const auto& r_communicator = r_model_part.GetCommunicator();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();

    const double integration_domain_size = CalculateIntegrationDomainSize();

    double area_weighted_phi = 0.0;
    if (mIsElementsConsidered) {
        area_weighted_phi += IndexPartition<IndexType>(r_model_part.NumberOfElements())
                           .for_each<SumReduction<double>>(Matrix(), [&](const IndexType iElement, Matrix& rTLS) {
                                return CalculateGeometryValueContribution((r_model_part.ElementsBegin() + iElement)->GetGeometry(), rTLS);
                           });
    }

    if (mIsConditionsConsidered) {
        area_weighted_phi += IndexPartition<IndexType>(r_model_part.NumberOfConditions())
                           .for_each<SumReduction<double>>(Matrix(), [&](const IndexType iCondition, Matrix& rTLS) {
                                return CalculateGeometryValueContribution((r_model_part.ConditionsBegin() + iCondition)->GetGeometry(), rTLS);
                           });
    }

    area_weighted_phi = r_data_communicator.SumAll(area_weighted_phi);
    area_weighted_phi /= integration_domain_size;

    return area_weighted_phi;

    KRATOS_CATCH("");
}

template<unsigned int TDim>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateIntegrationDomainSize() const
{
    KRATOS_TRY

    const auto& r_model_part = mrModelPart.GetSubModelPart(mIntegrationDomainModelPartName);
    const auto& r_communicator = r_model_part.GetCommunicator();
    const auto& r_data_communicator = r_communicator.GetDataCommunicator();

    double integration_domain_size = 0.0;
    if (mIsElementsConsidered) {
        integration_domain_size += IndexPartition<IndexType>(r_model_part.NumberOfElements())
                           .for_each<SumReduction<double>>([&r_model_part](const IndexType iElement) {
                               return (r_model_part.ElementsBegin() + iElement)->GetGeometry().DomainSize();
                           });
    }

    if (mIsConditionsConsidered) {
        integration_domain_size += IndexPartition<IndexType>(r_model_part.NumberOfConditions())
                           .for_each<SumReduction<double>>([&r_model_part](const IndexType iCondition) {
                               return (r_model_part.ConditionsBegin() + iCondition)->GetGeometry().DomainSize();
                           });
    }

    return r_data_communicator.SumAll(integration_domain_size);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
double DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateGeometryValueContribution(
    const GeometryType& rGeometry,
    Matrix& rShapeFunctions) const
{
    KRATOS_TRY

    this->CalculateGeometryData(rGeometry, rShapeFunctions);
    const Vector& rN = row(rShapeFunctions, 0);

    array_1d<double, 3> phi;
    FluidCalculationUtilities::EvaluateInPoint(
        rGeometry, rN, std::tie(phi, *mpVariable));

    return rGeometry.DomainSize() * std::pow(inner_prod(phi, phi), mPower);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
void DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<TDim>::CalculateGeometryData(
        const GeometryType& rGeometry,
        Matrix& rNContainer) const
{
    const auto r_integrations_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    const unsigned int number_of_gauss_points =
        rGeometry.IntegrationPointsNumber(r_integrations_method);

    const std::size_t number_of_nodes = rGeometry.PointsNumber();

    if (rNContainer.size1() != number_of_gauss_points ||
        rNContainer.size2() != number_of_nodes) {
        rNContainer.resize(number_of_gauss_points, number_of_nodes, false);
    }
    rNContainer = rGeometry.ShapeFunctionsValues(r_integrations_method);
}

// template instantiations
template class DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<2>;
template class DomainIntegrated3DArrayMagnitudeSquarePMeanResponseFunction<3>;

}