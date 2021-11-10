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

#if !defined(KRATOS_DOMAIN_INTEGRATED_SQUARE_MEAN_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_DOMAIN_INTEGRATED_SQUARE_MEAN_RESPONSE_FUNCTION_H_INCLUDED

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
#include "utilities/element_size_calculator.h"
#include "utilities/openmp_utils.h"

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
class DomainIntegratedSquareMeanResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    using ConditionType = ModelPart::ConditionType;

    using GeometryType = typename ConditionType::GeometryType;

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(DomainIntegratedSquareMeanResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DomainIntegratedSquareMeanResponseFunction(Parameters Settings, ModelPart& rModelPart)
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
        mpVariable = &KratosComponents<Variable<array_1d<double, 3>>>::Get(Settings["variable_name"].GetString());

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~DomainIntegratedSquareMeanResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        const int num_threads = ParallelUtilities::GetNumThreads();
        mShapeFunctions.resize(num_threads);

        if (mStructureModelPartName != "") {
            VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Conditions());
            VariableUtils().SetFlag(STRUCTURE, true, mrModelPart.GetSubModelPart(mStructureModelPartName).Conditions());
        }

        CalculateDomainArea();

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep() override
    {
        KRATOS_TRY

        mSquareMean = CalculateValue(mrModelPart);

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
        if (rResponseGradient.size() != rResidualGradient.size2()) {
            rResponseGradient.resize(rResidualGradient.size2());
        }

        rResponseGradient.clear();
    }

    void CalculateFirstDerivativesGradient(
        const Condition& rAdjointCondition,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY

        const auto& r_geometry = rAdjointCondition.GetGeometry();
        const IndexType number_of_nodes = r_geometry.PointsNumber();
        const IndexType block_size = rResidualGradient.size2() / number_of_nodes;
        const IndexType domain_size = rProcessInfo[DOMAIN_SIZE];
        const IndexType skip_size = block_size - domain_size;

        if (rResponseGradient.size() != rResidualGradient.size2()) {
            rResponseGradient.resize(rResidualGradient.size2());
        }

        if (rAdjointCondition.Is(STRUCTURE)) {
            const int k = OpenMPUtils::ThisThread();
            auto& r_shape_function = mShapeFunctions[k];

            const auto& r_geometry = rAdjointCondition.GetGeometry();

            this->CalculateGeometryData(r_geometry, r_shape_function);
            const Vector& N = row(r_shape_function, 0);

            const double area = std::pow(
                ElementSizeCalculator<2, 3>::AverageElementSize(r_geometry), 2);

            array_1d<double, 3> phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(phi, *mpVariable));


            IndexType local_index = 0;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                for (IndexType k = 0; k < domain_size; ++k) {
                    rResponseGradient[local_index++] = 2.0 * area * phi[k] * N[c] / mDomainArea;
                }
                local_index += skip_size;
            }

        } else {
            rResponseGradient.clear();
        }

        KRATOS_CATCH("");
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

        if (rSensitivityGradient.size() != rSensitivityMatrix.size1()) {
            rSensitivityGradient.resize(rSensitivityMatrix.size1(), false);
        }

        rSensitivityGradient.clear();


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

        if (rAdjointCondition.Is(STRUCTURE)) {
            const int k = OpenMPUtils::ThisThread();
            auto& r_shape_function = mShapeFunctions[k];

            const auto& r_geometry = rAdjointCondition.GetGeometry();
            const IndexType number_of_nodes = r_geometry.PointsNumber();
            const IndexType domain_size = rProcessInfo[DOMAIN_SIZE];

            this->CalculateGeometryData(r_geometry, r_shape_function);
            const Vector& N = row(r_shape_function, 0);

            array_1d<double, 3> phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(phi, *mpVariable));

            const double inner_phi = inner_prod(phi, phi);

            IndexType local_index = 0;
            for (IndexType c = 0; c < number_of_nodes; ++c) {
                for (IndexType k = 0; k < domain_size; ++k) {
                    rSensitivityGradient[local_index++] = CalculateGeometrySpecificDomainAreaDerivative<2, 3>(rAdjointCondition, c, k) * (inner_phi - mSquareMean) / mDomainArea;
                }
            }

        } else {
            rSensitivityGradient.clear();
        }
    }

    double CalculateValue(ModelPart&) override
    {
        KRATOS_TRY

        auto& r_model_part = mrModelPart.GetSubModelPart(mStructureModelPartName);
        auto& r_communicator = r_model_part.GetCommunicator();
        auto& r_data_communicator = r_communicator.GetDataCommunicator();

        KRATOS_ERROR_IF(r_communicator.GlobalNumberOfConditions() == 0)
            << r_model_part.FullName() << " has no conditions.\n";

        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

        // get first condition to decide number of nodes in the condition
        int number_of_nodes = 0;
        if (r_communicator.LocalMesh().Conditions().size() != 0) {
            number_of_nodes = r_communicator.LocalMesh().Conditions().begin()->GetGeometry().PointsNumber();
        }
        number_of_nodes = r_data_communicator.MaxAll(number_of_nodes);

        double area_weighted_phi = 0.0;
        if (domain_size == 3) {
            if (number_of_nodes == 3) {
                area_weighted_phi = CalculateGeometrySpecificValue<2, 3>();
            } else if (number_of_nodes == 4) {
                area_weighted_phi = CalculateGeometrySpecificValue<2, 4>();
            } else {
                KRATOS_ERROR << "Unsupported geometry type found in "
                             << r_model_part.FullName() << " in 3D. Only supports triangle and quads as surface conditions in 3D. [ number of nodes found in geometry = "
                             << number_of_nodes << " ].\n";
            }
            area_weighted_phi = r_data_communicator.SumAll(area_weighted_phi);
        } else {
            KRATOS_ERROR << "Unsupported domain size provided. Only supports 3D. [ DOMAIN_SIZE = " << domain_size << " ].\n";
        }

        area_weighted_phi /= mDomainArea;

        return area_weighted_phi;

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mStructureModelPartName;
    const Variable<array_1d<double, 3>> *mpVariable;

    double mDomainArea;
    double mSquareMean;

    std::vector<Matrix> mShapeFunctions;

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

    void CalculateDomainArea()
    {
        KRATOS_TRY

        auto& r_model_part = mrModelPart.GetSubModelPart(mStructureModelPartName);
        auto& r_communicator = r_model_part.GetCommunicator();
        auto& r_data_communicator = r_communicator.GetDataCommunicator();

        KRATOS_ERROR_IF(r_communicator.GlobalNumberOfConditions() == 0)
            << r_model_part.FullName() << " has no conditions.\n";

        const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

        // get first condition to decide number of nodes in the condition
        int number_of_nodes = 0;
        if (r_communicator.LocalMesh().Conditions().size() != 0) {
            number_of_nodes = r_communicator.LocalMesh().Conditions().begin()->GetGeometry().PointsNumber();
        }
        number_of_nodes = r_data_communicator.MaxAll(number_of_nodes);

        if (domain_size == 3) {
            if (number_of_nodes == 3) {
                mDomainArea = CalculateGeometrySpecificDomainArea<2, 3>();
            } else if (number_of_nodes == 4) {
                mDomainArea = CalculateGeometrySpecificDomainArea<2, 4>();
            } else {
                KRATOS_ERROR << "Unsupported geometry type found in "
                             << r_model_part.FullName() << " in 3D. Only supports triangle and quads as surface conditions in 3D. [ number of nodes found in geometry = "
                             << number_of_nodes << " ].\n";
            }
            mDomainArea = r_data_communicator.SumAll(mDomainArea);
        } else {
            KRATOS_ERROR << "Unsupported domain size provided. Only supports 3D. [ DOMAIN_SIZE = " << domain_size << " ].\n";
        }

        KRATOS_CATCH("");
    }

    template<unsigned int TDim, unsigned int TNumNodes>
    double CalculateGeometrySpecificDomainArea()
    {
        KRATOS_TRY

        auto& r_model_part = mrModelPart.GetSubModelPart(mStructureModelPartName);

        return block_for_each<SumReduction<double>>(r_model_part.Conditions(), [](const ConditionType& rCondition) {
            return std::pow(ElementSizeCalculator<TDim, TNumNodes>::AverageElementSize(rCondition.GetGeometry()), 2);
        });

        KRATOS_CATCH("");
    }

    template<unsigned int TDim, unsigned int TNumNodes>
    double CalculateGeometrySpecificDomainAreaDerivative(
        ConditionType& rCondition,
        const unsigned int DerivativeNodeIndex,
        const unsigned int DerivativeDirectionIndex)
    {
        KRATOS_TRY

        using element_size_calculator = ElementSizeCalculator<TDim, TNumNodes>;

        const auto& r_geometry = rCondition.GetGeometry();

        return 2.0 * element_size_calculator::AverageElementSize(r_geometry) *
               element_size_calculator::AverageElementSizeDerivative(
                   DerivativeNodeIndex, DerivativeDirectionIndex, r_geometry);

        KRATOS_CATCH("");
    }

    template<unsigned int TDim, unsigned int TNumNodes>
    double CalculateGeometrySpecificValue()
    {
        KRATOS_TRY;

        auto& r_model_part = mrModelPart.GetSubModelPart(mStructureModelPartName);

        return block_for_each<SumReduction<double>>(r_model_part.Conditions(), Matrix(), [&](ConditionType& rCondition, Matrix& rNs) {
            const auto& r_geometry = rCondition.GetGeometry();

            this->CalculateGeometryData(r_geometry, rNs);
            const Vector& N = row(rNs, 0);

            const double area = std::pow(ElementSizeCalculator<TDim, TNumNodes>::AverageElementSize(r_geometry), 2);

            array_1d<double, 3> phi;
            FluidCalculationUtilities::EvaluateInPoint(
                r_geometry, N, std::tie(phi, *mpVariable));

            return area * inner_prod(phi, phi);
        });

        KRATOS_CATCH("");
    }

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_DOMAIN_INTEGRATED_SQUARE_MEAN_RESPONSE_FUNCTION_H_INCLUDED defined */
