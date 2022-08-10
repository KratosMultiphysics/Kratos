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

#if !defined(KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED)
#define KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED

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
class DragResponseFunction : public AdjointResponseFunction
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DragResponseFunction);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DragResponseFunction(Parameters Settings, ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
        KRATOS_TRY;

        Parameters default_settings(R"(
        {
            "structure_model_part_name": "PLEASE_SPECIFY_STRUCTURE_MODEL_PART",
            "drag_direction": [1.0, 0.0, 0.0],
            "start_time"    : 0.0
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        mStructureModelPartName = Settings["structure_model_part_name"].GetString();
        mStartTime = Settings["start_time"].GetDouble();

        if (Settings["drag_direction"].IsArray() == false ||
            Settings["drag_direction"].size() != 3)
        {
            KRATOS_ERROR << "Invalid \"drag_direction\"." << std::endl;
        }

        for (unsigned int d = 0; d < TDim; ++d)
            mDragDirection[d] = Settings["drag_direction"][d].GetDouble();

        if (std::abs(norm_2(mDragDirection) - 1.0) > 1e-3)
        {
            const double magnitude = norm_2(mDragDirection);
            if (magnitude == 0.0)
                KRATOS_ERROR << "\"drag_direction\" is zero." << std::endl;

            KRATOS_WARNING("DragResponseFunction")
                << "Non unit magnitude in \"drag_direction\"." << std::endl;
            KRATOS_WARNING("DragResponseFunction") << "Normalizing ..." << std::endl;

            for (unsigned int d = 0; d < TDim; ++d)
                mDragDirection[d] /= magnitude;
        }

        KRATOS_CATCH("");
    }

    DragResponseFunction(
        ModelPart& rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    ~DragResponseFunction() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        Check();

        VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Nodes());
        VariableUtils().SetFlag(
            STRUCTURE, true,
            mrModelPart.GetSubModelPart(mStructureModelPartName).Nodes());

        KRATOS_CATCH("");
    }

    void CalculateGradient(const Element& rAdjointElement,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        CalculateDragContribution(
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

    void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override
    {
        CalculateDragContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient, rProcessInfo);
    }

    void CalculateFirstDerivativesGradient(const Condition& rAdjointCondition,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override
    {
        CalculateDragContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient, rProcessInfo);
    }

    void CalculateSecondDerivativesGradient(const Condition& rAdjointCondition,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override
    {
        rResponseGradient.clear();
    }

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        CalculateDragContribution(
            rSensitivityMatrix, rAdjointElement.GetGeometry().Points(), rSensitivityGradient, rProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculatePartialSensitivity(Condition& rAdjointCondition,
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
            const double local_drag = block_for_each<SumReduction<double>>(rModelPart.GetCommunicator().LocalMesh().Nodes(), [&](const ModelPart::NodeType &rNode) {
                const auto& reaction = rNode.FastGetSolutionStepValue(REACTION);
                const double coeff = static_cast<double>(rNode.Is(STRUCTURE));

                double nodal_drag_contribution = 0.0;
                for (unsigned int i = 0; i < TDim; ++i) {
                    nodal_drag_contribution += mDragDirection[i] * reaction[i] * coeff;
                }

                return nodal_drag_contribution; });
            return rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_drag);
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
    array_1d<double, TDim> mDragDirection;
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

    void CalculateDragContribution(const Matrix& rDerivativesOfResidual,
                                   const Element::NodesArrayType& rNodes,
                                   Vector& rDerivativesOfDrag,
                                   const ProcessInfo& rProcessInfo) const
    {
        if (rDerivativesOfDrag.size() != rDerivativesOfResidual.size1())
            rDerivativesOfDrag.resize(rDerivativesOfResidual.size1(), false);

        const double current_time = rProcessInfo[TIME];

        if (mStartTime <= current_time) {
            constexpr std::size_t max_size = 50;
            BoundedVector<double, max_size> drag_flag_vector(rDerivativesOfResidual.size2());
            drag_flag_vector.clear();

            const unsigned num_nodes = rNodes.size();
            const unsigned local_size = rDerivativesOfResidual.size2() / num_nodes;

            for (unsigned i_node = 0; i_node < num_nodes; ++i_node)
            {
                if (rNodes[i_node].Is(STRUCTURE))
                {
                    for (unsigned d = 0; d < TDim; ++d)
                        drag_flag_vector[i_node * local_size + d] = mDragDirection[d];
                }
            }

            noalias(rDerivativesOfDrag) = prod(rDerivativesOfResidual, drag_flag_vector);
        } else {
            rDerivativesOfDrag.clear();
        }
    }

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED defined */
