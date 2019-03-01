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
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"
#include "utilities/variable_utils.h"
#include "response_functions/adjoint_response_function.h"

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
            "drag_direction": [1.0, 0.0, 0.0]
        })");

        Settings.ValidateAndAssignDefaults(default_settings);

        mStructureModelPartName =
            Settings["structure_model_part_name"].GetString();

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

            KRATOS_WARNING("DragResponseFunction") << "Non unit magnitude in \"drag_direction\"." << std::endl;
            KRATOS_WARNING("DragResponseFunction") << "Normalizing ..." << std::endl;

            for (unsigned int d = 0; d < TDim; ++d)
                mDragDirection[d] /= magnitude;
        }

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~DragResponseFunction() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        KRATOS_TRY;

        Check();

        VariableUtils().SetFlag(STRUCTURE, false, mrModelPart.Nodes());
        VariableUtils().SetFlag(
            STRUCTURE, true, mrModelPart.GetSubModelPart(mStructureModelPartName).Nodes());

        KRATOS_CATCH("");
    }

    void CalculateGradient(const Element& rAdjointElement,
                           const Matrix& rResidualGradient,
                           Vector& rResponseGradient,
                           const ProcessInfo& rProcessInfo) override
    {
        CalculateDragContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient);
    }

    void CalculateFirstDerivativesGradient(const Element& rAdjointElement,
                                           const Matrix& rResidualGradient,
                                           Vector& rResponseGradient,
                                           const ProcessInfo& rProcessInfo) override
    {
        CalculateDragContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient);
    }

    void CalculateSecondDerivativesGradient(const Element& rAdjointElement,
                                            const Matrix& rResidualGradient,
                                            Vector& rResponseGradient,
                                            const ProcessInfo& rProcessInfo) override
    {
        CalculateDragContribution(
            rResidualGradient, rAdjointElement.GetGeometry().Points(), rResponseGradient);
    }

    void CalculatePartialSensitivity(Element& rAdjointElement,
                                     const Variable<array_1d<double, 3>>& rVariable,
                                     const Matrix& rSensitivityMatrix,
                                     Vector& rSensitivityGradient,
                                     const ProcessInfo& rProcessInfo) override
    {
        KRATOS_TRY;

        CalculateDragContribution(
            rSensitivityMatrix, rAdjointElement.GetGeometry().Points(), rSensitivityGradient);

        KRATOS_CATCH("");
    }

    double CalculateValue() override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "DragResponseFunction::CalculateValue() is not implemented!!!\n";
        KRATOS_CATCH("");
    }

    ///@}

protected:
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    std::string mStructureModelPartName;
    array_1d<double, TDim> mDragDirection;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
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
                                   Vector& rDerivativesOfDrag) const
    {
        constexpr std::size_t max_size = 50;
        BoundedVector<double, max_size> drag_flag_vector(rDerivativesOfResidual.size2());

        const unsigned num_nodes = rNodes.size();
        unsigned local_index = 0;
        for (unsigned i_node = 0; i_node < num_nodes; ++i_node)
        {
            if (rNodes[i_node].Is(STRUCTURE))
            {
                for (unsigned d = 0; d < TDim; ++d)
                    drag_flag_vector[local_index++] = mDragDirection[d];
            }
            else
            {
                for (unsigned int d = 0; d < TDim; ++d)
                    drag_flag_vector[local_index++] = 0.0;
            }

            drag_flag_vector[local_index++] = 0.0; // pressure dof
        }

        if (rDerivativesOfDrag.size() != rDerivativesOfResidual.size1())
            rDerivativesOfDrag.resize(rDerivativesOfResidual.size1(), false);

        noalias(rDerivativesOfDrag) = prod(rDerivativesOfResidual, drag_flag_vector);
    }

    ///@}
};

///@} // Kratos Classes

///@} // FluidDynamicsApplication group

} /* namespace Kratos.*/

#endif /* KRATOS_DRAG_RESPONSE_FUNCTION_H_INCLUDED defined */
