//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Jordi Cotela
//

#if !defined(KRATOS_FLUID_ELEMENT_UTILITIES_H )
#define  KRATOS_FLUID_ELEMENT_UTILITIES_H

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/define.h"
#include "includes/node.h"

// Application includes
#include "custom_utilities/fluid_element_data.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Auxiliary and specialized functions for elements derived from FluidElement.
template< std::size_t TNumNodes >
class FluidElementUtilities {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FluidElementUtilities
    KRATOS_CLASS_POINTER_DEFINITION(FluidElementUtilities);

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using ShapeDerivatives2DType = typename FluidElementData<2,TNumNodes,false>::ShapeDerivativesType;
    using ShapeDerivatives3DType = typename FluidElementData<3,TNumNodes,false>::ShapeDerivativesType;

    constexpr static std::size_t VoigtVector2DSize = FluidElementData<2,TNumNodes,false>::StrainSize;
    constexpr static std::size_t VoigtVector3DSize = FluidElementData<3,TNumNodes,false>::StrainSize;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    FluidElementUtilities() = delete;

    /// Deleted copy constructor.
    FluidElementUtilities(FluidElementUtilities const& rOther) = delete;

    /// Destructor.
    ~FluidElementUtilities();

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    FluidElementUtilities& operator=(FluidElementUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
    * Auxiliary function that writes the strain matrix (B) relating nodal degrees of freedom and the symmetric gradient of velocity.
    * Note that pressure Dofs are considered included in the array of nodal Dofs (and the corresponding rows set to zero).
    * 2D variant.
    * @param rDNDX Shape function gradients evaluated at the current integration point.
    * @param rStrainMatrix computed strain matrix for the current integration point (output).
    */
    static void GetStrainMatrix(
        const ShapeDerivatives2DType& rDNDX,
        BoundedMatrix<double, VoigtVector2DSize, 3 * TNumNodes>& rStrainMatrix);

    /**
    * Auxiliary function that writes the strain matrix (B) relating nodal degrees of freedom and the symmetric gradient of velocity.
    * Note that pressure Dofs are considered included in the array of nodal Dofs (and the corresponding rows set to zero).
    * 2D variant.
    * @param rDNDX Shape function gradients evaluated at the current integration point.
    * @param rStrainMatrix computed strain matrix for the current integration point (output).
    */
    static void GetStrainMatrix(
        const ShapeDerivatives3DType& rDNDX,
        BoundedMatrix<double, VoigtVector3DSize, 4 * TNumNodes>& rStrainMatrix);

    /**
    * Auxiliary function that writes the constitutive matrix (C) for a Newtonian fluid using the given dynamic viscosity (mu).
    * 2D variant.
    * @param DynamicViscosity Dynamic viscosity (mu) for the fluid.
    * @param rConstitutiveMatrix computed constitutive matrix for the fluid (output).
    */
    static void GetNewtonianConstitutiveMatrix(
        const double DynamicViscosity,
        BoundedMatrix<double, VoigtVector2DSize, VoigtVector2DSize>& rConstitutiveMatrix);

    /**
    * Auxiliary function that writes the constitutive matrix (C) for a Newtonian fluid using the given dynamic viscosity (mu).
    * 3D variant.
    * @param DynamicViscosity Dynamic viscosity (mu) for the fluid.
    * @param rConstitutiveMatrix computed constitutive matrix for the fluid (output).
    */
    static void GetNewtonianConstitutiveMatrix(
        const double DynamicViscosity,
        BoundedMatrix<double, VoigtVector3DSize, VoigtVector3DSize>& rConstitutiveMatrix);

    /**
     * This function transforms a vector n into a matrix P that can be used to compute
     * matrix-vector product A*n for a symmetric matrix A expressed in Voigt notation.
     * If a is the Voigt-notation representation of matrix A, A*n == P*a.
     * This is the 2D variant.
     * @param rVector the vector to transform.
     * @param rVoigtMatrix the transformed matrix (output).
     */
    static void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double, 2, VoigtVector2DSize>& rVoigtMatrix);

    /**
     * This function transforms a vector n into a matrix P that can be used to compute
     * matrix-vector product A*n for a symmetric matrix A expressed in Voigt notation.
     * If a is the Voigt-notation representation of matrix A, A*n == P*a.
     * This is the 3D variant.
     * @param rVector the vector to transform.
     * @param rVoigtMatrix the transformed matrix (output).
     */
    static void VoigtTransformForProduct(
        const array_1d<double,3>& rVector,
        BoundedMatrix<double, 3, VoigtVector3DSize>& rVoigtMatrix);

    /**
     * This function sets the normal projection matrix as the given unit normal outer product.
     * @param rUnitNormal reference to the unit normal vector
     * @param rNormProjMatrix reference to the normal projection matrix
     */
    static void SetNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, 2, 2>& rNormalProjMatrix);

    /**
     * This function sets the normal projection matrix as the given unit normal outer product.
     * @param rUnitNormal reference to the unit normal vector
     * @param rNormalProjMatrix reference to the normal projection matrix
     */
    static void SetNormalProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, 3, 3>& rNormalProjMatrix);

    /**
     * This function sets the tangential projection matrix as the identity matrix minus the given unit normal outer product.
     * @param rUnitNormal reference to the unit normal vector
     * @param rTangProjMatrix reference to the normal projection matrix
     */
    static void SetTangentialProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, 2, 2>& rTangProjMatrix);

    /**
     * This function sets the tangential projection matrix as the identity matrix minus the given unit normal outer product.
     * @param rUnitNormal reference to the unit normal vector
     * @param rTangProjMatrix reference to the normal projection matrix
     */
    static void SetTangentialProjectionMatrix(
        const array_1d<double, 3>& rUnitNormal,
        BoundedMatrix<double, 3, 3>& rTangProjMatrix);

    /**
     * Invert a system with 2 unknowns, defined using bounded matrix types.
     * @param[in] rA the 2x2 system matrix.
     * @param[in] rB the right hand side vector.
     * @param[out] rX the solution of the system.
     */
    static void DenseSystemSolve(
        const BoundedMatrix<double,2,2> &rA,
        const array_1d<double,2> &rB,
        array_1d<double,2> &rX);

    /**
     * Invert a system with 3 unknowns, defined using bounded matrix types.
     * @param[in] rA the 3x3 system matrix.
     * @param[in] rB the right hand side vector.
     * @param[out] rX the solution of the system.
     */
    static void DenseSystemSolve(
        const BoundedMatrix<double,3,3> &rA,
        const array_1d<double,3> &rB,
        array_1d<double,3> &rX);

    /**
     * @brief Evaluates given list of variable pairs at gauss point locations at step
     *
     * Example:
     *      double density;
     *      array_1d<double, 3> velocity
     *      EvaluateInPoint(rGeometry, rShapeFunction, Step,
     *                      std::tie(density, DENSITY),
     *                      std::tie(velocity, VELOCITY));
     *
     *      The above evaluation will fill density, and velocity variables with gauss point
     *      evaluated DENSITY and VELOCITY at gauss point with shape function values rShapeFunction
     *      in the geometry named rGeometry.
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunction            Shape function values at gauss point
     * @param[in] Step                      Step to be evaluated
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType, const Variable<TDataType>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void EvaluateInPoint(
        const GeometryType& rGeometry,
        const Vector& rShapeFunction,
        const int Step,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        KRATOS_TRY

        const auto& r_node = rGeometry[0];
        const double shape_function_value = rShapeFunction[0];

        int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
            std::get<0>(rValueVariablePairs) =
                r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value,
            0)...};

        // this can be removed with fold expressions in c++17
        *dummy = 0;

        for (unsigned int c = 1; c < TNumNodes; ++c) {
            const auto& r_node = rGeometry[c];
            const double shape_function_value = rShapeFunction[c];

            int dummy[sizeof...(TRefVariableValuePairArgs)] = {(
                UpdateValue<typename std::remove_reference<typename std::tuple_element<0, TRefVariableValuePairArgs>::type>::type>(
                    std::get<0>(rValueVariablePairs),
                    r_node.FastGetSolutionStepValue(std::get<1>(rValueVariablePairs), Step) * shape_function_value),
                0)...};

            // this can be removed with fold expressions in c++17
            *dummy = 0;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Evaluates given list of variable pairs at gauss point locations at current step
     *
     * Example:
     *      double density;
     *      array_1d<double, 3> velocity
     *      EvaluateInPoint(rGeometry, rShapeFunction,
     *                      std::tie(density, DENSITY),
     *                      std::tie(velocity, VELOCITY));
     *
     *      The above evaluation will fill density, and velocity variables with gauss point
     *      evaluated DENSITY and VELOCITY at gauss point with shape function values rShapeFunction
     *      in the geometry named rGeometry.
     *
     * @tparam TRefVariableValuePairArgs
     * @param[in] rGeometry                 Geometry to get nodes
     * @param[in] rShapeFunction            Shape function values at gauss point
     * @param[in/out] rValueVariablePairs   std::tuple<TDataType, const Variable<TDataType>> list of variables.
     */
    template <class... TRefVariableValuePairArgs>
    static void inline EvaluateInPoint(
        const GeometryType& rGeometry,
        const Vector& rShapeFunction,
        const TRefVariableValuePairArgs&... rValueVariablePairs)
    {
        EvaluateInPoint<TRefVariableValuePairArgs...>(rGeometry, rShapeFunction, 0, rValueVariablePairs...);
    }

    template<class TVectorType>
    static void Product(
        TVectorType& rOutput,
        const Matrix& rA,
        const Vector& rB)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rOutput.size() < rA.size1())
            << "rOutput.size() does not have enough indices to comply with "
               "matrix vector multiplication with rA. [ rOutput.size() = "
            << rOutput.size() << ", rA.size1() = " << rA.size1() << " ].\n";

        rOutput.clear();

        const IndexType cols = std::min(rA.size2(), rB.size());
        for (IndexType i = 0; i < rA.size1(); ++i) {
            for (IndexType j = 0; j < cols; ++j) {
                rOutput[i] += rA(i, j) * rB[j];
            }
        }

        KRATOS_CATCH("");
    }

    template<class TMatrixType>
    static void Product(
        TMatrixType& rOutput,
        const Matrix& rA,
        const Matrix& rB)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rOutput.size1() < rA.size1())
            << "rOutput.size1() does not have enough indices to comply with "
               "matrix matrix multiplication with rA and rB. [ rOutput.size1() = "
            << rOutput.size1() << ", rA.size1() = " << rA.size1() << " ].\n";

        KRATOS_DEBUG_ERROR_IF(rOutput.size2() < rB.size2())
            << "rOutput.size2() does not have enough indices to comply with "
               "matrix matrix multiplication with rA and rB. [ rOutput.size2() = "
            << rOutput.size2() << ", rB.size2() = " << rB.size2() << " ].\n";

        rOutput.clear();

        const IndexType cols = std::min(rA.size2(), rB.size1());
        for (IndexType i = 0; i < rA.size1(); ++i) {
            for (IndexType j = 0; j < rB.size2(); ++j) {
                for (IndexType k = 0; k < cols; ++k) {
                    rOutput(i, j) += rA(i, k) * rB(k, j);
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///@}

private:
    ///@name Private Operations
    ///@{

    template<class TDataType, class TDummy = void>
    static void UpdateValue(
        TDataType& rOutput,
        const TDataType& rInput)
    {
        noalias(rOutput) += rInput;
    }

    template<class TDummy>
    static void UpdateValue(
        double& rOutput,
        const double& rInput)
    {
        rOutput += rInput;
    }

    ///@}

};  // Class FluidElementUtilities

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FLUID_ELEMENT_UTILITIES_H  defined


