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

#if !defined(KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

class GaussPointKreisselmeierAggregationResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(GaussPointKreisselmeierAggregationResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    GaussPointKreisselmeierAggregationResponseFunction(
        Parameters Settings,
        ModelPart& rModelPart);

    /// Destructor.
    ~GaussPointKreisselmeierAggregationResponseFunction() override = default;
    ///@}
    ///@name Operations
    ///@{

    void Initialize() override;

    double CalculateValue(ModelPart& rModelPart) override;

    void CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculateSecondDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Element& rAdjointElement,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    void CalculatePartialSensitivity(
        Condition& rAdjointCondition,
        const Variable<array_1d<double, 3>>& rVariable,
        const Matrix& rSensitivityMatrix,
        Vector& rSensitivityGradient,
        const ProcessInfo& rProcessInfo) override;

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mCriticalModelPartName;
    std::string mPerturbationVariableName;
    double mGaussPointValueScalingFactor;
    double mStepSize;

    const Variable<double>* mpGaussPointValueScalarVariable;
    std::vector<const Variable<double>*> mPrimalStateScalarVariablePointersList;
    std::vector<const Variable<double>*> mPrimalStateFirstDerivativeScalarVariablePointersList;
    std::vector<const Variable<double>*> mPrimalStateSecondDerivativeScalarVariablePointersList;

    std::map<int,double> mKSPrefactors;
    double mSumKSPrefactors;
    double mRho;

    int mEchoLevel;

    bool mAreKSPrefactorsInitialized;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CalculateFiniteDifferenceStateVariableSensitivities(
        Vector& rOutput,
        ModelPart::ElementType& rElement,
        const std::vector<const Variable<double>*>& rDerivativeVariablePointersList,
        const ProcessInfo& rProcessInfo);

    void CalculateFiniteDifferenceShapeVariableSensitivities(
        Vector& rOutput,
        ModelPart::ElementType& rElement,
        const ProcessInfo& rProcessInfo);

    void CalculateStateDerivative(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const std::vector<const Variable<double>*>& rDerivativeVariablesList,
        const ProcessInfo& rProcessInfo);

    template<class TEntityType, class PerturbationMethod, class ValueMethod>
    void CalculateFiniteDifferenceSensitivity(
        Matrix& rOutput,
        TEntityType& rEntity,
        ValueMethod&& rValueMethod,
        PerturbationMethod&& rPerturbationMethod,
        const std::vector<const Variable<double>*>& rDerivativeVariablesList,
        const std::vector<double>& rPerturbationScalingFactor,
        const double Tolerance,
        const IndexType MaxNumberOfIterations)
    {
        KRATOS_TRY

        using value_type = Vector;

        value_type ref_value;
        rValueMethod(ref_value);

        auto& r_geometry = rEntity.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType number_of_variables = rDerivativeVariablesList.size();
        const IndexType number_of_dofs = number_of_nodes * number_of_variables;
        const IndexType number_of_values = ref_value.size();

        if (rOutput.size1() != number_of_dofs || rOutput.size2() !=  number_of_values) {
            rOutput.resize(number_of_dofs, number_of_values, false);
        }

        #pragma omp critical
        {
            value_type perturbed_value(number_of_values), previous_derivative(number_of_values, 0.0), current_derivative(number_of_values);
            for (IndexType i_node = 0; i_node < number_of_nodes; ++i_node) {
                auto& r_node = r_geometry[i_node];
                for (IndexType i_var = 0; i_var < number_of_variables; ++i_var) {
                    const auto p_variable = rDerivativeVariablesList[i_var];

                    double delta = 1e-12;
                    IndexType number_of_iterations = 0;
                    double dx = 1.0;
                    while (dx > Tolerance * norm_2(current_derivative) && number_of_iterations++ < MaxNumberOfIterations) {
                        rPerturbationMethod(r_node, *p_variable, +delta);
                        rValueMethod(perturbed_value);
                        rPerturbationMethod(r_node, *p_variable, -delta);

                        noalias(current_derivative) = (perturbed_value - ref_value) / delta;
                        noalias(previous_derivative) = previous_derivative - current_derivative;
                        dx = std::pow(inner_prod(previous_derivative, previous_derivative), 0.5);
                        noalias(previous_derivative) = current_derivative;

                        delta *= rPerturbationScalingFactor[i_var];
                        // KRATOS_WATCH(perturbed_value);
                        // KRATOS_WATCH(ref_value);
                        // KRATOS_WATCH(current_derivative);
                        // KRATOS_WATCH(number_of_iterations);
                        // KRATOS_WATCH(delta);
                        // KRATOS_WATCH(dx);

                    }

                    KRATOS_WATCH("Done var");

                    row(rOutput, i_node * number_of_variables + i_var) = current_derivative;
                    current_derivative.clear();
                    // std::exit(-1);
                    KRATOS_WARNING_IF("GaussPointKreisselmeierAggregationResponseFunction", number_of_iterations >= MaxNumberOfIterations) << "Max number of iterations reached";
                }

                KRATOS_WATCH("Done node");
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_GAUSS_POINT_KREISSELMEIER_AGGREGATION_RESPONSE_FUNCTION_H_INCLUDED defined */
