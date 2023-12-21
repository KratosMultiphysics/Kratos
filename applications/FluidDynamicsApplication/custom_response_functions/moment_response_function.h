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

#if !defined(KRATOS_MOMENT_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_MOMENT_RESPONSE_FUNCTION_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/ublas_interface.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

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
template <unsigned int TDim>
class MomentResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(MomentResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    MomentResponseFunction(Parameters Settings, ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "structure_model_part_name": "PLEASE_SPECIFY_STRUCTURE_MODEL_PART",
            "moment_direction": [1.0, 0.0, 0.0],
            "moment_point"    : [0.0, 0.0, 0.0],
            "start_time"    : 0.0
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        mStructureModelPartName = Settings["structure_model_part_name"].GetString();
        mStartTime = Settings["start_time"].GetDouble();

        if (Settings["moment_direction"].IsArray() == false || Settings["moment_direction"].size() != 3) {
            KRATOS_ERROR << "Invalid \"moment_direction\"." << std::endl;
        }
        mMomentDirection = Settings["moment_direction"].GetVector();

        if (std::abs(norm_2(mMomentDirection) - 1.0) > 1e-3) {
            const double magnitude = norm_2(mMomentDirection);
            KRATOS_ERROR_IF(magnitude == 0.0) << "\"moment_direction\" is zero." << std::endl;
            KRATOS_WARNING("MomentResponseFunction") << "Non unit magnitude in \"moment_direction\"." << std::endl;
            KRATOS_WARNING("MomentResponseFunction") << "Normalizing ..." << std::endl;
            mMomentDirection /= magnitude;
        }

        if (Settings["moment_point"].IsArray() == false || Settings["moment_point"].size() != 3) {
            KRATOS_ERROR << "Invalid \"moment_point\"." << std::endl;
        }
        mMomentPoint = Settings["moment_point"].GetVector();

        KRATOS_CATCH("");
    }

    MomentResponseFunction(
        ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~MomentResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        Check();

        VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Nodes());
        VariableUtils().SetFlag(STRUCTURE, true, mrModelPart.GetSubModelPart(mStructureModelPartName).Nodes());

        KRATOS_CATCH("");
    }

    void CalculateGradient(
        const Element& rAdjointElement,
        const Matrix& rResidualGradient,
        Vector& rResponseGradient,
        const ProcessInfo& rProcessInfo) override
    {
        CalculateMomentContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient, rProcessInfo);
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
        CalculateMomentContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient, rProcessInfo);
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
        CalculateMomentContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient, rProcessInfo);
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

        CalculateMomentShapeContribution(
            rSensitivityMatrix, rAdjointElement.GetGeometry().Points(), rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
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
        KRATOS_TRY;

        const double current_time = rModelPart.GetProcessInfo()[TIME];

        if (mStartTime <= current_time) {
            const double local_moment = block_for_each<SumReduction<double>>(rModelPart.GetCommunicator().LocalMesh().Nodes(), [&](const ModelPart::NodeType &rNode) {
                const double coeff = static_cast<double>(rNode.Is(STRUCTURE));
                const array_1d<double, 3>& reaction = rNode.FastGetSolutionStepValue(REACTION);
                const array_1d<double, 3>& r = rNode.Coordinates() - mMomentPoint;

                array_1d<double, 3> moment;
                moment[0] = r[1]*reaction[2] - r[2]*reaction[1];
                moment[1] = r[2]*reaction[0] - r[0]*reaction[2];
                moment[2] = r[0]*reaction[1] - r[1]*reaction[0];

                double nodal_moment_contribution = 0.0;
                for (unsigned int i = 0; i < 3; ++i) {
                    nodal_moment_contribution += mMomentDirection[i] * moment[i] * coeff;
                }

                return nodal_moment_contribution; });
            return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_moment);
        } else {
            return 0.0;
        }

        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mStructureModelPartName;
    array_1d<double, 3> mMomentDirection;
    array_1d<double, 3> mMomentPoint;
    double mStartTime;

    ///@}
    ///@name Protected Operations
    ///@{

    void Check()
    {
        KRATOS_TRY;

        if (mrModelPart.HasSubModelPart(mStructureModelPartName) == false)
            KRATOS_ERROR << "No sub model part \"" << mStructureModelPartName
                         << "\"" << std::endl;

        KRATOS_CATCH("");
    }

    void CalculateMomentContribution(
        const Matrix& rDerivativesOfResidual,
        const Element::NodesArrayType& rNodes,
        Vector& rDerivativesOfMoment,
        const ProcessInfo& rProcessInfo) const
    {
        constexpr std::size_t max_size = 50;
        if (rDerivativesOfMoment.size() != rDerivativesOfResidual.size1())
            rDerivativesOfMoment.resize(rDerivativesOfResidual.size1(), false);

        const double current_time = rProcessInfo[TIME];

        if (mStartTime <= current_time) {
            const unsigned num_nodes = rNodes.size();
            const unsigned local_size = rDerivativesOfResidual.size2() / num_nodes;
            const unsigned modified_local_size = (TDim == 3) ? local_size : (local_size + 1);

            BoundedMatrix<double, max_size, max_size> moment_matrix(rDerivativesOfResidual.size1(), modified_local_size * num_nodes);
            BoundedVector<double, max_size> moment_flag_vector(modified_local_size * num_nodes);
            moment_matrix.clear();
            moment_flag_vector.clear();

            for (unsigned i_node = 0; i_node < num_nodes; ++i_node) {
                const auto& r_node = rNodes[i_node];
                if (r_node.Is(STRUCTURE)) {
                    const array_1d<double, 3>& r = r_node.Coordinates() - mMomentPoint;
                    for (unsigned c = 0; c < rDerivativesOfResidual.size1(); ++c) {
                        if constexpr(TDim == 3) {
                            moment_matrix(c, i_node * modified_local_size + 0) = rDerivativesOfResidual(c, i_node * local_size + 2) * r[1]  - rDerivativesOfResidual(c, i_node * local_size + 1) * r[2];
                            moment_matrix(c, i_node * modified_local_size + 1) = rDerivativesOfResidual(c, i_node * local_size + 0) * r[2]  - rDerivativesOfResidual(c, i_node * local_size + 2) * r[0];
                        }
                        moment_matrix(c, i_node * modified_local_size + 2) = rDerivativesOfResidual(c, i_node * local_size + 1) * r[0]  - rDerivativesOfResidual(c, i_node * local_size + 0) * r[1];
                    }

                    for (unsigned d = 0; d < 3; ++d) {
                        moment_flag_vector[i_node * modified_local_size + d] = mMomentDirection[d];
                    }
                }
            }

            noalias(rDerivativesOfMoment) = prod(moment_matrix, moment_flag_vector);
        } else {
            rDerivativesOfMoment.clear();
        }
    }

    void CalculateMomentShapeContribution(
        const Matrix& rDerivativesOfResidual,
        const Element::NodesArrayType& rNodes,
        Vector& rDerivativesOfMoment,
        const ProcessInfo& rProcessInfo) const
    {
        constexpr std::size_t max_size = 50;
        if (rDerivativesOfMoment.size() != rDerivativesOfResidual.size1())
            rDerivativesOfMoment.resize(rDerivativesOfResidual.size1(), false);

        const double current_time = rProcessInfo[TIME];

        if (mStartTime <= current_time) {
            const unsigned num_nodes = rNodes.size();
            const unsigned residual_local_size = rDerivativesOfResidual.size2() / num_nodes;
            const unsigned modified_residual_local_size = (TDim == 3) ? residual_local_size : (residual_local_size + 1);
            const unsigned shape_local_size = rDerivativesOfResidual.size1() / num_nodes;

            BoundedMatrix<double, max_size, max_size> moment_matrix(rDerivativesOfResidual.size1(), modified_residual_local_size * num_nodes);
            BoundedVector<double, max_size> moment_flag_vector(modified_residual_local_size * num_nodes);
            moment_matrix.clear();
            moment_flag_vector.clear();

            for (unsigned i_node = 0; i_node < num_nodes; ++i_node) {
                const auto& r_node = rNodes[i_node];
                if (r_node.Is(STRUCTURE)) {
                    const array_1d<double, 3>& r = r_node.Coordinates() - mMomentPoint;
                    const array_1d<double, 3>& reaction = r_node.FastGetSolutionStepValue(REACTION);

                    for (unsigned c = 0; c < rDerivativesOfResidual.size1(); ++c) {
                        if constexpr(TDim == 3) {
                            moment_matrix(c, i_node * modified_residual_local_size + 0) = rDerivativesOfResidual(c, i_node * residual_local_size + 2) * r[1]  - rDerivativesOfResidual(c, i_node * residual_local_size + 1) * r[2];
                            moment_matrix(c, i_node * modified_residual_local_size + 1) = rDerivativesOfResidual(c, i_node * residual_local_size + 0) * r[2]  - rDerivativesOfResidual(c, i_node * residual_local_size + 2) * r[0];
                        }
                        moment_matrix(c, i_node * modified_residual_local_size + 2) = rDerivativesOfResidual(c, i_node * residual_local_size + 1) * r[0]  - rDerivativesOfResidual(c, i_node * residual_local_size + 0) * r[1];
                    }

                    if constexpr(TDim == 3) {
                        moment_matrix(i_node * shape_local_size + 1, i_node * modified_residual_local_size + 0) -= reaction[2];
                        moment_matrix(i_node * shape_local_size + 2, i_node * modified_residual_local_size + 0) += reaction[1];

                        moment_matrix(i_node * shape_local_size + 2, i_node * modified_residual_local_size + 1) -= reaction[0];
                        moment_matrix(i_node * shape_local_size + 0, i_node * modified_residual_local_size + 1) += reaction[2];
                    }

                    moment_matrix(i_node * shape_local_size + 0, i_node * modified_residual_local_size + 2) -= reaction[1];
                    moment_matrix(i_node * shape_local_size + 1, i_node * modified_residual_local_size + 2) += reaction[0];

                    for (unsigned d = 0; d < 3; ++d) {
                        moment_flag_vector[i_node * modified_residual_local_size + d] = mMomentDirection[d];
                    }
                }
            }

            noalias(rDerivativesOfMoment) = prod(moment_matrix, moment_flag_vector);
        } else {
            rDerivativesOfMoment.clear();
        }
    }

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_MOMENT_RESPONSE_FUNCTION_H_INCLUDED defined */
