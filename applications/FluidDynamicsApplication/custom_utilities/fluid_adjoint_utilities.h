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

#if !defined(KRATOS_FLUID_ADJOINT_UTILITIES_H)
#define KRATOS_FLUID_ADJOINT_UTILITIES_H

// System includes
#include <array>

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"
#include "utilities/coordinate_transformation_utilities.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos classes
///@{

template<unsigned int TDim>
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) FluidAdjointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using IndexType = std::size_t;

    using GeometryType = Geometry<NodeType>;

    using CoordinateTransformationUtilities = CoordinateTransformationUtils<Matrix, Vector, double>;

    ///@}
    ///@name Static Operations
    ///@{

    template <class TDataType>
    static const Variable<double>& GetRelevantVariable(
        const IndexType DirectionIndex,
        const Variable<TDataType>& rVariable,
        const std::array<const Variable<double>*, 3>& rAllVariableComponents);

    template <class TDataType, std::size_t TGradientVariableTotalDimensionality>
    static std::array<const Variable<double>*, TDim> GetRelevantGradientVariableComponentList(
        const IndexType DirectionIndex,
        const Variable<TDataType>& rVariable,
        const std::array<const Variable<double>*, TGradientVariableTotalDimensionality>& rAllGradientVariableComponents);

    ///@}
    ///@name Classes
    ///@{

    template <unsigned int TBlockSize>
    class SlipUtilities
    {
    public:
        ///@name Static Operations
        ///@{

        /**
         * @brief Calculates rotated slip applied state derivatives
         *
         * This method calculates rotated slip applied state derivatives
         * when non-rotated residual derivatives are given. Nodal rotations
         * are determined by SLIP flag. If SLIP flag is true, then it is assumed
         * to be rotated as well as slip condition applied.
         *
         * This method assumes first dofs to be VELOCITY, and derivatives starting with
         * w.r.t. VELOCITY
         *
         * @param rOutput                   Rotated and slip applied state derivatives
         * @param rResidualDerivatives      Non-rotated and non-slip applied state derivatives
         * @param rGeometry                 Geometry of which rResidualDerivatives are computed on
         */
        static void CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const GeometryType& rGeometry);

        /**
         * @brief Calculates rotated slip applied state derivatives
         *
         * This method calculates rotated slip applied state derivatives
         * when non-rotated residual derivatives are given. Nodal rotations
         * are determined by SLIP flag. If SLIP flag is true, then it is assumed
         * to be rotated as well as slip condition applied.
         *
         * This method assumes first dofs to be VELOCITY, and derivatives not starting with
         * w.r.t. VELOCITY because SLIP is based on the VELOCITY variable.
         *
         * @param rOutput                   Rotated and slip applied state derivatives
         * @param rResidualDerivatives      Non-rotated and non-slip applied state derivatives
         * @param rGeometry                 Geometry of which rResidualDerivatives are computed on
         */
        static void CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const GeometryType& rGeometry);

        /**
         * @brief Rotates residual derivatives
         *
         * This method is used to rotate part of the non-rotated residual derivatives
         * which corresponds to a given node.
         *
         * @param rOutput                   Nodal rotated state derivatives
         * @param rResidualDerivatives      Non-rotated state derivatives
         * @param NodeStartIndex            Block starting column for that specific node
         * @param rNode                     Node which has rotated residual derivative contributions
         */
        static void AddNodalRotationDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const IndexType NodeStartIndex,
            const NodeType& rNode);

        /**
         * @brief Slip condition residual derivatives
         *
         * This method is used to apply slip condition state derivatives to part of the non-rotated residual derivatives
         * which corresponds to a given node.
         *
         * @param rOutput                   Nodal slip condition applied state derivatives
         * @param NodeStartIndex            Block starting column for that specific node
         * @param rNode                     Node which has rotated residual derivative contributions
         */
        static void AddNodalApplySlipConditionDerivatives(
            Matrix& rOutput,
            const IndexType NodeStartIndex,
            const NodeType& rNode);

        /**
         * @brief Adds the nodal state derivatives as it is
         *
         * This method is used to apply slip condition state derivatives to part of the non-rotated residual derivatives
         * which corresponds to a given node.
         *
         * @param rOutput                   Output state derivatives
         * @param rResidualDerivatives      Non-rotated state derivatives
         * @param NodeStartIndex            Block starting column for that specific node
         */
        static void AddNodalResidualDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const IndexType NodeStartIndex);

        /**
         * @brief Clears a residual derivative
         *
         * @param rOutput               Output matrix
         * @param ResidualIndex         Residual index which needs to be set to zero.
         */
        static void ClearNodalResidualDerivatives(
            Matrix& rOutput,
            const IndexType ResidualIndex);

        ///@}
    private:
        ///@name Private Static Members
        ///@{

        static const CoordinateTransformationUtilities mRotationTool;

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_UTILITIES_H
