// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES)
#define KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES

#include "utilities/math_utils.h"

// TODO: Remove repeated functions

namespace Kratos
{
class StructuralMechanicsMathUtilities
{
public:

    ///@name Type definitions
    ///@{
     
    typedef long double RealType;
    
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

//        KRATOS_WATCH(t1g);
//        KRATOS_WATCH(t2g);
//        KRATOS_WATCH(t3g);
    }

    static inline void Comp_Orthonor_Base(
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & t,
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

        t(2, 0) = t3g[0];
        t(2, 1) = t3g[1];
        t(2, 2) = t3g[2];

        MathUtils<double>::CrossProduct(t2g, t3g, vxe);
        n = norm_2(t2g);
        t2g /= n;

        t(1, 0) = t2g[0];
        t(1, 1) = t2g[1];
        t(1, 2) = t2g[2];

        MathUtils<double>::CrossProduct(t1g, t2g, t3g);
        n = norm_2(t1g);
        t1g /= n;

        t(0, 0) = t1g[0];
        t(0, 1) = t1g[1];
        t(0, 2) = t1g[2];
    }

    /**
     * Gives the interpolation of two Gauss Points (to introduce in GiD the result)
     * @param nG Number of Gauss Points
     * @return Matrix of interpolations
     */
    static inline Matrix interpol_PrismGiD(const int nG)
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

    static inline bool Solve_Second_Order_Equation(
            const RealType& a,
            const RealType& b,
            const RealType& c,
            std::vector<RealType>& solution
            )
    {
        RealType disc = b*b - 4.00*a*c;
        RealType q = 0.00;

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
    
private:
};// class StructuralMechanicsMathUtilities
}
#endif /* KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES defined */
 
