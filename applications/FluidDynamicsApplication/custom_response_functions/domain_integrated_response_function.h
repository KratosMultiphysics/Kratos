//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_DOMAIN_INTEGRATED_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_DOMAIN_INTEGRATED_RESPONSE_FUNCTION_H_INCLUDED

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
class DomainIntegratedResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using ElementType = ModelPart::ElementType;

    using GeometryType = typename ElementType::GeometryType;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(DomainIntegratedResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DomainIntegratedResponseFunction(Parameters Settings, ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "model_part_name": "",
            "variable_name"  : "PLEASE_SPECIFY_VARIABLE_NAME"
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        mStructureModelPartName = Settings["model_part_name"].GetString();

        pVariable = &KratosComponents<Variable<double>>::Get(Settings["variable_name"].GetString());

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~DomainIntegratedResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        if (mStructureModelPartName != "") {
            VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Elements());
            VariableUtils().SetFlag(STRUCTURE, true, mrModelPart.GetSubModelPart(mStructureModelPartName).Elements());
        }

        KRATOS_CATCH("");
    }

    void CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        if (rResponseGradient.size() != rResidualGradient.size1()) {
            rResponseGradient.resize(rResidualGradient.size1(), false);
        }

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
        KRATOS_TRY

        const auto& r_geometry = rAdjointElement.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType block_size = rResidualGradient.size2() / number_of_nodes;
        const IndexType skip_size = block_size - 1;

        if (rResponseGradient.size() != rResidualGradient.size2()) {
            rResponseGradient.resize(rResidualGradient.size2());
        }

        rResponseGradient.clear();

        if (rAdjointElement.Is(STRUCTURE)) {
            Matrix shape_functions;
            this->CalculateGeometryData(r_geometry, shape_functions);
            const Vector& N = row(shape_functions, 0);

            double phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(phi, *pVariable));

            const double volume = r_geometry.DomainSize();

            IndexType local_index = 0;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                // adding velocity potential derivative
                rResponseGradient[local_index] = volume * N[c];

                // skipping rest of the derivatives if they are used
                local_index += skip_size + 1;
            }
        } else {
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        if (rResponseGradient.size() != rResidualGradient.size1()) {
            rResponseGradient.resize(rResidualGradient.size1(), false);
        }

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
        KRATOS_TRY;

        const auto& r_geometry = rAdjointElement.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType domain_size = rProcessInfo[DOMAIN_SIZE];
        const IndexType local_size = domain_size * number_of_nodes;

        if (rSensitivityGradient.size() != local_size) {
            rSensitivityGradient.resize(local_size);
        }

        if (rAdjointElement.Is(STRUCTURE)) {
            Matrix shape_functions;
            this->CalculateGeometryData(r_geometry, shape_functions);
            const Vector& N = row(shape_functions, 0);

            double phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(phi, *pVariable));

            const double x10 = r_geometry[1].X() - r_geometry[0].X();
            const double y10 = r_geometry[1].Y() - r_geometry[0].Y();

            const double x20 = r_geometry[2].X() - r_geometry[0].X();
            const double y20 = r_geometry[2].Y() - r_geometry[0].Y();

            for (IndexType c = 0; c < number_of_nodes; ++c) {
                for (IndexType k = 0; k < domain_size; ++k) {
                    const double x10_derivative = ((c == 1) - (c == 0)) * (k == 0);
                    const double y10_derivative = ((c == 1) - (c == 0)) * (k == 1);
                    const double x20_derivative = ((c == 2) - (c == 0)) * (k == 0);
                    const double y20_derivative = ((c == 2) - (c == 0)) * (k == 1);

                    double detJ_derivative = 0.0;
                    detJ_derivative += x10 * y20_derivative;
                    detJ_derivative += x10_derivative * y20;
                    detJ_derivative -= y10 * x20_derivative;
                    detJ_derivative -= y10_derivative * x20;

                    rSensitivityGradient[c * domain_size + k] = 0.5 * detJ_derivative * phi;
                }
            }
        } else {
            rSensitivityGradient.clear();
        }

        KRATOS_CATCH("");
    }

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override
    {
        if (rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        }

        rSensitivityGradient.clear();
    }

    double CalculateValue(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const double local_value = block_for_each<SumReduction<double>>(rModelPart.Elements(), Matrix(), [&](ModelPart::ElementType& rElement, Matrix& rNs) {
            if (rElement.Is(STRUCTURE)) {
                const auto& r_geometry = rElement.GetGeometry();

                this->CalculateGeometryData(r_geometry, rNs);
                const Vector& N = row(rNs, 0);

                const double volume = r_geometry.DomainSize();

                double phi;
                FluidCalculationUtilities::EvaluateInPoint(
                    r_geometry, N, std::tie(phi, *pVariable));

                return volume * phi;
            } else {
                return 0.0;
            }
        });

        return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_value);

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mStructureModelPartName;
    const Variable<double> *pVariable;

    ///@}
    ///@name Protected Operations
    ///@{

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

#endif /* KRATOS_DOMAIN_INTEGRATED_RESPONSE_FUNCTION_H_INCLUDED defined */
