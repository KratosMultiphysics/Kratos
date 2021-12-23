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

#if !defined(KRATOS_VELOCITY_PRESSURE_NORM_SQUARE_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_VELOCITY_PRESSURE_NORM_SQUARE_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "containers/model.h"
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/element_size_calculator.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// A response function for norm square calculation.
/**
 * The response function is defined as:
 *
 * \f[
 * \bar{D} = \Sigma_{n=1}^N F_u\left(A_n|\underline{u}_n|\right)^2 + F_p\left(A_nP\right)^2 \Delta t
 * \f]
 *
 * if "integrate_in_time" is true.
 *
 * F_u is the velocity norm factor
 * F_p is the pressure norm factor
 */
class VelocityPressureNormSquareResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using ElementType = ModelPart::ElementType;

    using GeometryType = typename ElementType::GeometryType;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(VelocityPressureNormSquareResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VelocityPressureNormSquareResponseFunction(
        Parameters Settings,
        Model& rModel)
        : mrModel(rModel)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "main_model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
            "norm_model_part_name": "PLEASE_SPECIFY_MODEL_PART_NAME",
            "velocity_norm_factor": 1.0,
            "pressure_norm_factor": 1.0,
            "entities"            : ["elements"]
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        mMainModelPartName = Settings["main_model_part_name"].GetString();
        mNormModelPartName = Settings["norm_model_part_name"].GetString();
        mVelocityNormFactor = Settings["velocity_norm_factor"].GetDouble();
        mPressureNormFactor = Settings["pressure_norm_factor"].GetDouble();

        mIsElements = false;
        mIsConditions = false;
        const auto& entities = Settings["entities"].GetStringArray();
        for (const auto& entity : entities) {
            if (entity == "elements") {
                mIsElements = true;
            } else if (entity == "conditions") {
                mIsConditions = true;
            } else {
                KRATOS_ERROR << "Unsupported entity type provided under "
                                "\"entities\". Supported entity types are:\n"
                             << "\t elements\n"
                             << "\t conditions\n.";
            }
        }

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~VelocityPressureNormSquareResponseFunction() override
    {
    }

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        auto& r_main_model_part = mrModel.GetModelPart(mMainModelPartName);
        auto& r_norm_model_part = mrModel.GetModelPart(mNormModelPartName);

        if (mIsElements && !mIsConditions) {
            const IndexType number_of_entities =
                r_norm_model_part.GetCommunicator().GlobalNumberOfElements();

            if (number_of_entities == 0) {
                KRATOS_WARNING("VelocityPressureNormSquareResponseFunction")
                    << "Only \"elements\" are chosen as entities "
                       "and no elements are found in "
                    << r_norm_model_part.Name() << " model part. Therefore trying to compute response function on conditions.\n";
                mIsConditions = true;
                mIsElements = false;
            }
        }

        VariableUtils().SetFlag(STRUCTURE, false, r_main_model_part.Elements());
        VariableUtils().SetFlag(STRUCTURE, false, r_main_model_part.Conditions());

        IndexType total_entities = 0;
        if (mIsElements) {
            total_entities += r_norm_model_part.GetCommunicator().GlobalNumberOfElements();
            VariableUtils().SetFlag(STRUCTURE, true, r_norm_model_part.Elements());
        }

        if (mIsConditions) {
            total_entities += r_norm_model_part.GetCommunicator().GlobalNumberOfConditions();
            VariableUtils().SetFlag(STRUCTURE, true, r_norm_model_part.Conditions());
        }

        KRATOS_ERROR_IF(total_entities == 0)
            << "No entities were found in "
            << r_norm_model_part.Name() << " to calculate response function value. Please check model part.\n";

        KRATOS_CATCH("");
    }

    void CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculateGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntityFirstDerivatives<Element>(
            rAdjointElement, rResidualGradient, rResponseGradient, rProcessInfo);
    }

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntityFirstDerivatives<Condition>(
            rAdjointCondition, rResidualGradient, rResponseGradient, rProcessInfo);
    }

    void CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntitySensitivityDerivatives(rAdjointElement, rSensitivityGradient, rProcessInfo);
    }

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntitySensitivityDerivatives(rAdjointCondition, rSensitivityGradient, rProcessInfo);
    }

    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        double norm_square = 0.0;
        if (mIsElements) {
            norm_square += block_for_each<SumReduction<double>>(
                rModelPart.Elements(), Matrix(),
                [&](ModelPart::ElementType& rElement, Matrix& rMatrix) -> double {
                    return CalculateEntityValue<Element>(rElement, rMatrix);
                });
        }

        if (mIsConditions) {
            norm_square += block_for_each<SumReduction<double>>(
                rModelPart.Conditions(), Matrix(),
                [&](ModelPart::ConditionType& rCondition, Matrix& rMatrix) -> double {
                    return CalculateEntityValue<Condition>(rCondition, rMatrix);
                });
        }

        return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(norm_square);

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    Model& mrModel;

    std::string mMainModelPartName;
    std::string mNormModelPartName;

    double mVelocityNormFactor;
    double mPressureNormFactor;
    bool mIsElements;
    bool mIsConditions;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    double CalculateEntityValue(
        const TEntityType& rEntity,
        Matrix& rShapeFunctions) const
    {
        if (rEntity.Is(STRUCTURE)) {
            const auto& r_geometry = rEntity.GetGeometry();

            this->CalculateGeometryData(r_geometry, rShapeFunctions);
            const Vector& N = row(rShapeFunctions, 0);

            const double volume = r_geometry.DomainSize();

            double pressure;
            array_1d<double, 3> velocity;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(velocity, VELOCITY), std::tie(pressure, PRESSURE));

            return mVelocityNormFactor * std::pow(volume * norm_2(velocity), 2) +
                   mPressureNormFactor * std::pow(volume * pressure, 2);
        } else {
            return 0;
        }
    }

    template<class TEntityType>
    void CalculateEntityFirstDerivatives(
        const TEntityType& rEntity,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_TRY

        const auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType domain_size = rProcessInfo[DOMAIN_SIZE];
        const IndexType block_size = rResidualGradient.size2() / number_of_nodes;
        const IndexType skip_size = block_size - domain_size - 1;

        if (rResponseGradient.size() != rResidualGradient.size2()) {
            rResponseGradient.resize(rResidualGradient.size2());
        }

        if (rEntity.Is(STRUCTURE)) {
            Matrix shape_functions;
            this->CalculateGeometryData(r_geometry, shape_functions);
            const Vector& N = row(shape_functions, 0);

            array_1d<double, 3> velocity;
            double pressure;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(velocity, VELOCITY), std::tie(pressure, PRESSURE));

            const double volume_2 = std::pow(r_geometry.DomainSize(), 2);

            const double coeff_1 = volume_2 * 2.0 * mVelocityNormFactor;
            const double coeff_2 = volume_2 * 2.0 * mPressureNormFactor * pressure;

            IndexType local_index = 0;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                // adding velocity derivatives
                for (IndexType k = 0; k < domain_size; ++k) {
                    rResponseGradient[local_index++] = coeff_1 * N[c] * velocity[k];
                }

                // adding pressure derivatives
                rResponseGradient[local_index] = coeff_2 * N[c];

                // skipping rest of the derivatives if they are used
                local_index += skip_size + 1;
            }
        } else {
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntitySensitivityDerivatives(
        const TEntityType& rEntity,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_TRY

        const auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType domain_size = rProcessInfo[DOMAIN_SIZE];
        const IndexType local_size = domain_size * number_of_nodes;

        if (rResponseGradient.size() != local_size) {
            rResponseGradient.resize(local_size);
        }

        if (rEntity.Is(STRUCTURE)) {
            Matrix shape_functions;
            this->CalculateGeometryData(r_geometry, shape_functions);
            const Vector& N = row(shape_functions, 0);

            const double volume = r_geometry.DomainSize();

            double pressure;
            array_1d<double, 3> velocity;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(velocity, VELOCITY), std::tie(pressure, PRESSURE));

            const double lx = r_geometry[0].X() - r_geometry[1].X();
            const double ly = r_geometry[0].Y() - r_geometry[1].Y();

            for (IndexType c = 0; c < number_of_nodes; ++c) {
                for (IndexType k = 0; k < domain_size; ++k) {
                    const double lx_derivative =
                        (c == 0 && k == 0) * (1.0) + (c == 1 && k == 0) * (-1.0);
                    const double ly_derivative =
                        (c == 0 && k == 1) * (1.0) + (c == 1 && k == 1) * (-1.0);
                    const double volume_derivative =
                        (1.0 / volume) * (lx * lx_derivative + ly * ly_derivative);

                    double value = 0.0;
                    value += mVelocityNormFactor * std::pow(norm_2(velocity), 2) *
                             2 * volume * volume_derivative;
                    value += mPressureNormFactor * std::pow(pressure, 2) * 2 *
                             volume * volume_derivative;

                    rResponseGradient[c * domain_size + k] = value;
                }
            }
        } else {
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void CalculateGeometryData(
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

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_VELOCITY_PRESSURE_NORM_SQUARE_RESPONSE_FUNCTION_H_INCLUDED defined */
