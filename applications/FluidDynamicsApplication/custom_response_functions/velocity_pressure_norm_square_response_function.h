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
#include <tuple>

// External includes

// Project includes
#include "geometries/geometry_data.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
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

/// A response function for drag.
/**
 * The response function is defined as:
 *
 * \f[
 * \bar{D} = \Sigma_{n=1}^N D^n \Delta t
 * \f]
 *
 * if "integrate_in_time" is true.
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
        ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "norm_model_part_name": "PLEASE_SPECIFY_STRUCTURE_MODEL_PART",
            "velocity_norm_factor": 1.0,
            "pressure_norm_factor": 1.0
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        mNormModelPartName = Settings["norm_model_part_name"].GetString();
        mVelocityNormFactor = Settings["velocity_norm_factor"].GetDouble();
        mPressureNormFactor = Settings["pressure_norm_factor"].GetDouble();

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

        KRATOS_ERROR_IF(!mrModelPart.HasSubModelPart(mNormModelPartName))
            << mNormModelPartName << " model part not found in "
            << mrModelPart.Name() << ".\n";

        auto& r_norm_model_part = mrModelPart.GetSubModelPart(mNormModelPartName);

        VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Elements());
        VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Conditions());
        VariableUtils().SetFlag(STRUCTURE, true, r_norm_model_part.Elements());
        VariableUtils().SetFlag(STRUCTURE, true, r_norm_model_part.Conditions());

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntityFirstDerivatives<Element>(
            rAdjointElement, rResponseGradient, rProcessInfo);
    }

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateEntityFirstDerivatives<Condition>(
            rAdjointCondition, rResponseGradient, rProcessInfo);
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
        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

        rSensitivityGradient.clear();
    }

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override
    {
        if (rSensitivityGradient.size() != rSensitivityMatrix.size1())
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);

        rSensitivityGradient.clear();
    }

    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        const double element_norm_square = block_for_each<SumReduction<double>>(
            rModelPart.Elements(), Matrix(),
            [&](ModelPart::ElementType& rElement, Matrix& rMatrix) -> double {
                return CalculateEntityValue<Element>(rElement, rMatrix);
            });

        const double condition_norm_square = block_for_each<SumReduction<double>>(
            rModelPart.Conditions(), Matrix(),
            [&](ModelPart::ConditionType& rCondition, Matrix& rMatrix) -> double {
                return CalculateEntityValue<Condition>(rCondition, rMatrix);
            });

        const double norm_square = element_norm_square + condition_norm_square;
        return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(norm_square);

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mNormModelPartName;
    double mVelocityNormFactor;
    double mPressureNormFactor;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    double CalculateEntityValue(
        const TEntityType& rEntity,
        Matrix& rShapeFunctions) const
    {
        if (rEntity.Is(STRUCTURE)) {
            this->CalculateGeometryData(rEntity.GetGeometry(), rShapeFunctions);
            const Vector& N = row(rShapeFunctions, 0);

            double pressure;
            array_1d<double, 3> velocity;
            FluidCalculationUtilities::EvaluateInPoint(
                rEntity.GetGeometry(), N, std::tie(velocity, VELOCITY),
                std::tie(pressure, PRESSURE));

            return mVelocityNormFactor * std::pow(norm_2(velocity), 2) +
                   mPressureNormFactor * std::pow(pressure, 2);
        } else {
            return 0;
        }
    }

    template<class TEntityType>
    void CalculateEntityFirstDerivatives(
        const TEntityType& rEntity,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) const
    {
        KRATOS_TRY

        const auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType domain_size = rProcessInfo[DOMAIN_SIZE];
        const IndexType block_size = domain_size + 1;
        const IndexType local_size = block_size * number_of_nodes;

        if (rResponseGradient.size() != local_size) {
            rResponseGradient.resize(local_size);
        }

        if (rEntity.Is(STRUCTURE)) {
            Matrix shape_functions;
            this->CalculateGeometryData(r_geometry, shape_functions);
            const Vector& N = row(shape_functions, 0);

            array_1d<double, 3> velocity;
            double pressure;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(velocity, VELOCITY), std::tie(pressure, PRESSURE));

            const double coeff_1 = 2.0 * mVelocityNormFactor;
            const double coeff_2 = 2.0 * mPressureNormFactor * pressure;

            IndexType local_index = 0;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                // adding velocity derivatives
                for (IndexType k = 0; k < domain_size; ++k) {
                    rResponseGradient[local_index++] = coeff_1 * N[c] * velocity[k];
                }

                // adding pressure derivatives
                rResponseGradient[local_index++] = coeff_2 * N[c];
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
