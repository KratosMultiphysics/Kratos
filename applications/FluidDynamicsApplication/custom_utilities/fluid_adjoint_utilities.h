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

    class KRATOS_API(FLUID_DYNAMICS_APPLICATION) SlipUtilities
    {
    public:
        ///@name Life Cycle
        ///@{

        SlipUtilities(const IndexType BlockSize);

        ///@}
        ///@name Operations
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
        void CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const GeometryType& rGeometry) const;

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
        void CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const GeometryType& rGeometry) const;

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
        void AddNodalRotationDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const IndexType NodeStartIndex,
            const NodeType& rNode) const;

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
        void AddNodalApplySlipConditionDerivatives(
            Matrix& rOutput,
            const IndexType NodeStartIndex,
            const NodeType& rNode) const;

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
        void AddNodalResidualDerivatives(
            Matrix& rOutput,
            const Matrix& rResidualDerivatives,
            const IndexType NodeStartIndex) const;

        /**
         * @brief Clears a residual derivative
         *
         * @param rOutput               Output matrix
         * @param ResidualIndex         Residual index which needs to be set to zero.
         */
        void ClearNodalResidualDerivatives(
            Matrix& rOutput,
            const IndexType ResidualIndex) const;

        ///@}
    private:
        ///@name Private Members
        ///@{

        const IndexType mBlockSize;

        const CoordinateTransformationUtilities mRotationTool;

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_UTILITIES_H
