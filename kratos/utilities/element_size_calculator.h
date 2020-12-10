//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#if !defined(KRATOS_ELEMENT_SIZE_CALCULATOR_H )
#define  KRATOS_ELEMENT_SIZE_CALCULATOR_H

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"


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

/// A collection of functions to compute element size, to be used as the h parameter in stabilized CFD formulations.
template< std::size_t TDim, std::size_t TNumNodes >
class KRATOS_API(KRATOS_CORE) ElementSizeCalculator {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ElementSizeCalculator
    KRATOS_CLASS_POINTER_DEFINITION(ElementSizeCalculator);

    typedef BoundedMatrix<double,TNumNodes,TDim> ShapeDerivativesType;

    ///@}
    ///@name Life Cycle
    ///@{

    // Deleted default constructor
    ElementSizeCalculator() = delete;

    /// Deleted copy constructor.
    ElementSizeCalculator(ElementSizeCalculator const& rOther) = delete;

    /// Destructor.
    ~ElementSizeCalculator();

    ///@}
    ///@name Operators
    ///@{

    /// Deleted assignment operator.
    ElementSizeCalculator& operator=(ElementSizeCalculator const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /// Minimum element size based on the geometry.
    /** @param rGeometry The geometry of calling element.
     *  @return The computed size.
     */
    static double MinimumElementSize(const Geometry<Node<3> >& rGeometry);

    /// Minimum element size derivative based on the geometry
    /** @param DerivativeNodeIndex NodeIndex for which the derivative is obtained
     *  @param DerivativeDirectionIndex DirectionIndex for which the derivative is obtained
     *  @param rGeometry The geometry of calling element
     *  @return The derivative w.r.t. DerivativeNodeIndex and DerivativeDirectionIndex
     */
    static double MinimumElementSizeDerivative(
        const unsigned int DerivativeNodeIndex,
        const unsigned int DerivativeDirectionIndex,
        const Geometry<Node<3> >& rGeometry);

    /// Average element size based on the geometry.
    /** @param rGeometry The geometry of calling element.
     *  @return The computed size.
     */
    static double AverageElementSize(const Geometry<Node<3> >& rGeometry);

    /// Average element size derivative based on the geometry
    /** @param DerivativeNodeIndex NodeIndex for which the derivative is obtained
     *  @param DerivativeDirectionIndex DirectionIndex for which the derivative is obtained
     *  @param rGeometry The geometry of calling element
     *  @return The derivative w.r.t. DerivativeNodeIndex and DerivativeDirectionIndex
     */
    static double AverageElementSizeDerivative(
        const unsigned int DerivativeNodeIndex,
        const unsigned int DerivativeDirectionIndex,
        const Geometry<Node<3> >& rGeometry);

    /// Projected element size in the direction of the velocity vector.
    /** @param rGeometry The geometry of calling element.
     *  @param rVelocity The velocity defining the direction of projection.
     *  @return The computed size.
     */
    static double ProjectedElementSize(const Geometry<Node<3> >& rGeometry, const array_1d<double,3>& rVelocity);

    /// Element size based on the shape functions gradients. Triangle element version.
    /** @param rDN_DX The shape functions gradients.
     *  @return The computed size.
     */
    static double GradientsElementSize(const BoundedMatrix<double, 3, 2>& rDN_DX);


    /// Element size based on the shape functions gradients. Tetrahedral element version.
    /** @param rDN_DX The shape functions gradients.
     *  @return The computed size.
     */
    static double GradientsElementSize(const BoundedMatrix<double, 4, 3>& rDN_DX);

    ///@}

private:
    ///@name Private Operations
    ///@{

    static inline double HsqDerivative2D(
        const double Vx,
        const double VxDerivative,
        const double Nx,
        const double NxDerivative,
        const double Vy,
        const double VyDerivative,
        const double Ny,
        const double NyDerivative)
    {
        return 2 * (Vx * Nx + Vy * Ny) *
                   (VxDerivative * Nx + Vx * NxDerivative + VyDerivative * Ny +
                    Vy * NyDerivative) / (Nx * Nx + Ny * Ny) -
               std::pow((Vx * Nx + Vy * Ny) / (Nx * Nx + Ny * Ny), 2) *
                   (2 * Nx * NxDerivative + 2 * Ny * NyDerivative);
    }

    static inline double HsqDerivative3D(
        const double Vx,
        const double VxDerivative,
        const double Nx,
        const double NxDerivative,
        const double Vy,
        const double VyDerivative,
        const double Ny,
        const double NyDerivative,
        const double Vz,
        const double VzDerivative,
        const double Nz,
        const double NzDerivative)
    {
        return 2 * (Vx * Nx + Vy * Ny + Vz * Nz) *
                   (VxDerivative * Nx + Vx * NxDerivative + VyDerivative * Ny +
                    Vy * NyDerivative + VzDerivative * Nz + Vz * NzDerivative) / (Nx * Nx + Ny * Ny + Nz * Nz) -
               std::pow((Vx * Nx + Vy * Ny + Vz * Nz) / (Nx * Nx + Ny * Ny + Nz * Nz), 2) *
                   (2 * Nx * NxDerivative + 2 * Ny * NyDerivative + 2 * Nz * NzDerivative);
    }

    static inline double NormSquareDerivative(
        const array_1d<double, 3>& rValue,
        const array_1d<double, 3>& rValueDerivative)
    {
        return 2 * inner_prod(rValue, rValueDerivative);
    }

    static inline double EdgeLengthDerivative(
        const unsigned int DerivativeNodeIndex,
        const unsigned int DerivativeDirectionIndex,
        const unsigned int NodeIndexA,
        const unsigned int NodeIndexB,
        const unsigned int EdgeDirection)
    {
        return ((DerivativeNodeIndex == NodeIndexA) - (DerivativeNodeIndex == NodeIndexB)) * (DerivativeDirectionIndex == EdgeDirection);
    }

    ///@}

};  // Class ElementSizeCalculator

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ELEMENT_SIZE_CALCULATOR_H  defined
