// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES)
#define KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES

#include "utilities/math_utils.h"

namespace Kratos
{
///@}
///@name  Enum's
///@{

#if !defined(INITIAL_CURRENT)
#define INITIAL_CURRENT
    enum Configuration {Initial = 0, Current = 1};
#endif

class StructuralMechanicsMathUtilities
{
public:

    ///@name Type definitions
    ///@{

    typedef long double                                RealType;

    typedef Node<3>                                    NodeType;

    typedef Geometry<NodeType>                     GeometryType;

    ///@}
    ///@name Operations
    ///@{

    /**
     * It gives the orthonormal base based in a orthogonal base
     * @return tig: The orthonormal base
     * @param vxe and vye: X and Y direction vector
     */
    static inline void Comp_Orthonor_Vect(
            array_1d<double, 3 > & t1g,
            array_1d<double, 3 > & t2g,
            array_1d<double, 3 > & t3g,
            const array_1d<double, 3 > & vxe,
            const array_1d<double, 3 > & vye
            )
    {
        double n;

        MathUtils<double>::CrossProduct(t3g, vxe, vye);
        n = norm_2(t3g);
        t3g /= n;

        MathUtils<double>::CrossProduct(t2g, t3g, vxe);
        n = norm_2(t2g);
        t2g /= n;

        MathUtils<double>::CrossProduct(t1g, t2g, t3g);
        n = norm_2(t1g);
        t1g /= n;
    }

    /**
     * It gives the orthonormal base based in the cartesian derivatives
     * @return tig: The orthonormal base
     * @param vxe: X direction vector
     * @param Xdxi and Xeta: Cartesian derivatives in xi and eta
     */
    static inline void Comp_Orthonor_Base(
            array_1d<double, 3 > & t1g,
            array_1d<double, 3 > & t2g,
            array_1d<double, 3 > & t3g,
            const array_1d<double, 3 > & vxe,
            const array_1d<double, 3 > & Xdxi,
            const array_1d<double, 3 > & Xdeta
            )
    {
        double n;

        MathUtils<double>::CrossProduct(t3g, Xdxi, Xdeta);
        n = norm_2(t3g);
        t3g /= n;

        MathUtils<double>::CrossProduct(t2g, t3g, vxe);
        n = norm_2(t2g);
        t2g /= n;

        MathUtils<double>::CrossProduct(t1g, t2g, t3g);
        n = norm_2(t1g);
        t1g /= n;
    }

    static inline void Comp_Orthonor_Base(
        BoundedMatrix<double, 3, 3 > & t,
        const array_1d<double, 3 > & vxe,
        const array_1d<double, 3 > & Xdxi,
        const array_1d<double, 3 > & Xdeta
        )
    {
        double n;

        array_1d<double, 3 > t1g, t2g, t3g;

        MathUtils<double>::CrossProduct(t3g, Xdxi, Xdeta);

        n = norm_2(t3g);
        t3g /= n;

        MathUtils<double>::CrossProduct(t2g, t3g, vxe);
        n = norm_2(t2g);
        t2g /= n;

        MathUtils<double>::CrossProduct(t1g, t2g, t3g);
        n = norm_2(t1g);
        t1g /= n;

        for (std::size_t i = 0; i < 3; ++i) {
            t(0, i) = t1g[i];
            t(1, i) = t2g[i];
            t(2, i) = t3g[i];
        }
    }

    /**
     * Gives the interpolation of two Gauss Points (to introduce in GiD the result)
     * @param nG Number of Gauss Points
     * @return Matrix of interpolations
     */
    static inline Matrix InterpolPrismGiD(const int nG)
    {
        Matrix interpol;
        interpol.resize(nG, 6, false);

        /* Assigning values to the vectors */
        if (nG == 1)
        {
            for (unsigned int i = 0; i < 6; i++)
            {
                interpol(0, i) = 1.0;
            }
        }
        else if (nG == 2)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                interpol(0, i) = 1.0;
                interpol(1, i) = 0.0;
            }
            for (unsigned int i = 3; i < 6; i++)
            {
                interpol(0, i) = 0.0;
                interpol(1, i) = 1.0;
            }
        }
        else if (nG == 3)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                interpol(0, i) = 0.745326;
                interpol(1, i) = 0.254644;
                interpol(2, i) = 0.0;
            }
            for (unsigned int i = 3; i < 6; i++)
            {
                interpol(0, i) = 0.0;
                interpol(1, i) = 0.254644;
                interpol(2, i) = 0.745326;
            }
        }
        else if (nG == 4)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                interpol(0, i) = 0.455467382132614037538;
                interpol(1, i) = 0.544532617867385962462;
                interpol(2, i) = 0.0;
                interpol(3, i) = 0.0;
            }
            for (unsigned int i = 3; i < 6; i++)
            {
                interpol(0, i) = 0.0;
                interpol(1, i) = 0.0;
                interpol(2, i) = 0.544532617867385962462;
                interpol(3, i) = 0.455467382132614037538;
            }
        }
        else if (nG == 5)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                interpol(0, i) = 0.062831503549096234958;
                interpol(1, i) = 0.907868;
                interpol(2, i) = 0.0293;
                interpol(3, i) = 0.0;
                interpol(4, i) = 0.0;
            }
            for (unsigned int i = 3; i < 6; i++)
            {
                interpol(0, i) = 0.0;
                interpol(1, i) = 0.0;
                interpol(2, i) = 0.0293;
                interpol(3, i) = 0.907868;
                interpol(4, i) = 0.062831503549096234958;
            }
        }
        else if (nG == 7)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                interpol(0, i) = 0.0;
                interpol(1, i) = 0.51090930312223869755;
                interpol(2, i) = 0.48909069687776130245;
                interpol(3, i) = 0.0;
                interpol(4, i) = 0.0;
                interpol(5, i) = 0.0;
                interpol(6, i) = 0.0;
            }
            for (unsigned int i = 3; i < 6; i++)
            {
                interpol(0, i) = 0.0;
                interpol(1, i) = 0.0;
                interpol(2, i) = 0.0;
                interpol(3, i) = 0.0;
                interpol(4, i) = 0.48909069687776130245;
                interpol(5, i) = 0.51090930312223869755;
                interpol(6, i) = 0.0;
            }
        }
        else if (nG == 11)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                interpol( 0, i) = 0.0;
                interpol( 1, i) = 0.0;
                interpol( 2, i) = 0.27601287860590845062;
                interpol( 3, i) = 0.72398712139409154938;
                interpol( 4, i) = 0.0;
                interpol( 5, i) = 0.0;
                interpol( 6, i) = 0.0;
                interpol( 7, i) = 0.0;
                interpol( 8, i) = 0.0;
                interpol( 9, i) = 0.0;
                interpol(10, i) = 0.0;
            }
            for (unsigned int i = 3; i < 6; i++)
            {
                interpol( 0, i) = 0.0;
                interpol( 1, i) = 0.0;
                interpol( 2, i) = 0.0;
                interpol( 3, i) = 0.0;
                interpol( 4, i) = 0.0;
                interpol( 5, i) = 0.0;
                interpol( 6, i) = 0.0;
                interpol( 7, i) = 0.72398712139409154938;
                interpol( 8, i) = 0.27601287860590845062;
                interpol( 9, i) = 0.0;
                interpol(10, i) = 0.0;
            }
        }
        return interpol;
    }

    /**
     * Calculates the solutions for a given second order polynomial equation 0 = a*x^2 + b*x + c
     * @param a coefficient
     * @param b coefficient
     * @param c coefficient
     * @param ZeroTol number treated as zero
     * @return Vector of solutions
     */

    static inline bool SolveSecondOrderEquation(
            const RealType& a,
            const RealType& b,
            const RealType& c,
            std::vector<RealType>& solution
            )
    {
        const RealType disc = b*b - 4.00*a*c;
        RealType q = 0.0;

        solution.resize(2, false);

        if (b > 0.00)
        {
            q = -0.5 * (b + std::sqrt(disc));
        }
        else
        {
            q = -0.5 * (b - std::sqrt(disc));
        }

        solution[0] = q / a;
        solution[1] = c / q;

        return true;
    }

    /**
     * Calculates the radius of axisymmetry
     * @param N: The Gauss Point shape function
     * @param Geom: The geometry studied
     * @return Radius: The radius of axisymmetry
     */

    static inline double CalculateRadius(
        const Vector N,
        GeometryType& Geom,
        const Configuration ThisConfiguration = Current
        )
    {
        double Radius = 0.0;

        for (unsigned int iNode = 0; iNode < Geom.size(); iNode++)
        {
            // Displacement from the reference to the current configuration
            if (ThisConfiguration == Current)
            {
                const array_1d<double, 3 > CurrentPosition = Geom[iNode].Coordinates();
                Radius += CurrentPosition[0] * N[iNode];
            }
            else
            {
                const array_1d<double, 3 > DeltaDisplacement = Geom[iNode].FastGetSolutionStepValue(DISPLACEMENT) - Geom[iNode].FastGetSolutionStepValue(DISPLACEMENT,1);
                const array_1d<double, 3 > CurrentPosition = Geom[iNode].Coordinates();
                const array_1d<double, 3 > ReferencePosition = CurrentPosition - DeltaDisplacement;
                Radius += ReferencePosition[0] * N[iNode];
            }
        }

        return Radius;
    }

    /**
     * Calculates the radius of axisymmetry for a point
     * @param Geom: The geometry studied
     * @return The radius of axisymmetry
     */

    static inline double CalculateRadiusPoint(
        GeometryType& Geom,
        const Configuration ThisConfiguration = Current
        )
    {
        // Displacement from the reference to the current configuration
        if (ThisConfiguration == Current)
        {
            const array_1d<double, 3 > CurrentPosition = Geom[0].Coordinates();
            return CurrentPosition[0];
        }
        else
        {
            const array_1d<double, 3 > DeltaDisplacement = Geom[0].FastGetSolutionStepValue(DISPLACEMENT) - Geom[0].FastGetSolutionStepValue(DISPLACEMENT,1);
            const array_1d<double, 3 > CurrentPosition = Geom[0].Coordinates();
            const array_1d<double, 3 > ReferencePosition = CurrentPosition - DeltaDisplacement;
            return ReferencePosition[0];
        }
    }

    /**
     * Transforms a 2-point tensor from an origin system to a target system
     * M=M_ij origin_left x origin_right = M_lk target_left X target_right
     * @return rOriginLeft: matrix with the basis vectors of basis origin left as columns
     * @return rOriginRight: matrix with the basis vectors of basis origin right as columns
     * @return rTargetLeft: matrix with the basis vectors of basis target left as columns
     * @return rTsargetRight: matrix with the basis vectors of basis target right as columns
     * @return rTensor: the tensor to be tranformed
     */
    template<int TDim>
    static inline void TensorTransformation(
        BoundedMatrix<double,TDim,TDim>& rOriginLeft,
        BoundedMatrix<double,TDim,TDim>& rOriginRight,
        BoundedMatrix<double,TDim,TDim>& rTargetLeft,
        BoundedMatrix<double,TDim,TDim>& rTargetRight,
        BoundedMatrix<double,TDim,TDim>& rTensor)
    {
        // metric computation (of the target systems)
        BoundedMatrix<double,TDim,TDim> metric_left = ZeroMatrix(TDim,TDim);
        BoundedMatrix<double,TDim,TDim> metric_right = ZeroMatrix(TDim,TDim);
        for(int i=0;i<TDim;i++){
            for(int j=0;j<TDim;j++){
                metric_left(i,j) += inner_prod(column(rTargetLeft,i),column(rTargetLeft,j));
                metric_right(i,j) += inner_prod(column(rTargetRight,i),column(rTargetRight,j));
            }
        }

        // invert metric
        double det;
        Matrix inv_metric_left = Matrix(TDim,TDim);
        Matrix inv_metric_right = Matrix(TDim,TDim);
        MathUtils<double>::InvertMatrix(Matrix(metric_left),inv_metric_left,det);
        MathUtils<double>::InvertMatrix(metric_right,inv_metric_right,det);

        // Compute dual target base vectors
        BoundedMatrix<double,TDim,TDim> target_left_dual = ZeroMatrix(TDim,TDim); // Anna noalias?
        BoundedMatrix<double,TDim,TDim> target_right_dual = ZeroMatrix(TDim,TDim); // Anna noalias?
        for(int i=0;i<TDim;i++){
            for(int j=0;j<TDim;j++){
                column(target_left_dual,i) += inv_metric_left(i,j)*column(rTargetLeft,j);
                column(target_right_dual,i) += inv_metric_right(i,j)*column(rTargetRight,j);
            }
        }

        // Tensor transformation
        BoundedMatrix<double, TDim, TDim> transformed_tensor = ZeroMatrix(TDim, TDim); // Anna noalias?
        for(int k=0;k<TDim;k++){
            for(int l=0;l<TDim;l++){
                for(int i=0;i<TDim;i++){
                    for(int j=0;j<TDim;j++){
                        transformed_tensor(k,l) += rTensor(i,j)*inner_prod(column(target_left_dual,k),column(rOriginLeft,i))*inner_prod(column(target_right_dual,l),column(rOriginRight,j));
                    }
                }
            }
        }
        rTensor = transformed_tensor;
    }

private:
};// class StructuralMechanicsMathUtilities
}
#endif /* KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES defined */

