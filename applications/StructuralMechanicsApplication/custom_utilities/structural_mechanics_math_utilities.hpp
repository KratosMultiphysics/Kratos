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
#define PI 3.1415926535898

#include "utilities/math_utils.h"
#include "geometries/point.h"
#include <algorithm> // To user std::sort()
#include <cmath>

namespace Kratos
{
class StructuralMechanicsMathUtilities
{
public:
    /**
     * @name Type definitions
     * @{
     */

    typedef boost::numeric::ublas::vector<Vector> Second_Order_Tensor;

    typedef boost::numeric::ublas::vector<Second_Order_Tensor> Third_Order_Tensor;

    typedef boost::numeric::ublas::vector<boost::numeric::ublas::vector<Matrix> > Fourth_Order_Tensor;

    typedef matrix<Second_Order_Tensor> Matrix_Second_Tensor;

    typedef long double RealType;

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

    /***********************************************************************************/
    /***********************************************************************************/

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

    /***********************************************************************************/
    /***********************************************************************************/

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

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives the x values of the Gauss-Legendre quadrature integral
     * @param nG Number of Gauss Points
     * @return Vector of coordinates x
     */
    static inline std::vector<double> GaussQuadrature_xcoor(const int nG)
    {
        // Note: To calculate more use http://keisan.casio.com/exec/system/1329114617
        std::vector<double> x_gauss(nG);
        /* Assigning values to the vectors */
        if (nG == 1)
        {
            x_gauss[0] = 0.0;
        }
        else if (nG == 2)
        {
            x_gauss[0] = -0.5773502691896257645091;
            x_gauss[1] = 0.57735026918962576450919;
        }
        else if (nG == 3)
        {
            x_gauss[0] = -0.7745966692414833770359;
            x_gauss[1] = 0.0;
            x_gauss[2] = 0.7745966692414833770359;
        }
        else if (nG == 4)
        {
            x_gauss[0] = -0.8611363115940525752239;
            x_gauss[1] = -0.3399810435848562648027;
            x_gauss[2] = 0.3399810435848562648027;
            x_gauss[3] = 0.8611363115940525752239;
        }
        else if (nG == 5)
        {
            x_gauss[0] = -0.9061798459386639927976;
            x_gauss[1] = -0.5384693101056830910363;
            x_gauss[2] = 0.0;
            x_gauss[3] = 0.5384693101056830910363;
            x_gauss[4] = 0.9061798459386639927976;
        }
        else if (nG == 6)
        {
            x_gauss[0] = -0.9324695142031520278123;
            x_gauss[1] = -0.6612093864662645136614;
            x_gauss[2] = -0.2386191860831969086305;
            x_gauss[3] = 0.2386191860831969086305;
            x_gauss[4] = 0.6612093864662645136614;
            x_gauss[5] = 0.9324695142031520278123;
        }
        else if (nG == 7)
        {
            x_gauss[0] = -0.9491079123427585245262;
            x_gauss[1] = -0.7415311855993944398639;
            x_gauss[2] = -0.4058451513773971669066;
            x_gauss[3] = 0.0;
            x_gauss[4] = 0.4058451513773971669066;
            x_gauss[5] = 0.7415311855993944398639;
            x_gauss[6] = 0.9491079123427585245262;
        }
        else if (nG == 8)
        {
            x_gauss[0] = -0.960289856497536231684;
            x_gauss[1] = -0.7966664774136267395916;
            x_gauss[2] = -0.5255324099163289858177;
            x_gauss[3] = -0.1834346424956498049395;
            x_gauss[4] = 0.1834346424956498049395;
            x_gauss[5] = 0.5255324099163289858177;
            x_gauss[6] = 0.7966664774136267395916;
            x_gauss[7] = 0.960289856497536231684;
        }
        else if (nG == 9)
        {
            x_gauss[0] = -0.968160239507626089836;
            x_gauss[1] = -0.8360311073266357942994;
            x_gauss[2] = -0.6133714327005903973087;
            x_gauss[3] = -0.3242534234038089290385;
            x_gauss[4] = 0.0;
            x_gauss[5] = 0.3242534234038089290385;
            x_gauss[6] = 0.6133714327005903973087;
            x_gauss[7] = 0.8360311073266357942994;
            x_gauss[8] = 0.968160239507626089836;
        }
        else if (nG == 10)
        {
            x_gauss[0] = -0.973906528517171720078;
            x_gauss[1] = -0.8650633666889845107321;
            x_gauss[2] = -0.6794095682990244062343;
            x_gauss[3] = -0.4333953941292471907993;
            x_gauss[4] = -0.148874338981631210885;
            x_gauss[5] = 0.148874338981631210885;
            x_gauss[6] = 0.4333953941292471907993;
            x_gauss[7] = 0.6794095682990244062343;
            x_gauss[8] = 0.8650633666889845107321;
            x_gauss[9] = 0.973906528517171720078;
        }
        else if (nG == 11)
        {
            x_gauss[0] = -0.9782286581460569928039;
            x_gauss[1] = -0.8870625997680952990752;
            x_gauss[2] = -0.7301520055740493240934;
            x_gauss[3] = -0.5190961292068118159257;
            x_gauss[4] = -0.2695431559523449723315;
            x_gauss[5] = 0.0000000000000000000000;
            x_gauss[6] = 0.2695431559523449723315;
            x_gauss[7] = 0.5190961292068118159257;
            x_gauss[8] = 0.7301520055740493240934;
            x_gauss[9] = 0.8870625997680952990752;
            x_gauss[10] = 0.9782286581460569928039;
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "The number of layers is bigger that the ones considered in the Gauss Quadrature", "");
        }

        return x_gauss;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Gives the weights of the Gauss-Legendre quadrature integral
     * @param nG Number of Gauss Points
     * @return Vector of weights
     */
    static inline std::vector<double> GaussQuadrature_weight(const int nG)
    {
        // To calculate more use http://keisan.casio.com/exec/system/1329114617
        std::vector<double> w_gauss(nG);
        /* Assigning values to the vectors */
        if (nG == 1)
        {
            w_gauss[0] = 2.0;
        }
        else if (nG == 2)
        {
            w_gauss[0] = 1.0;
            w_gauss[1] = 1.0;
        }
        else if (nG == 3)
        {
            w_gauss[0] = 0.5555555555555555555556;
            w_gauss[1] = 0.888888888888888888889;
            w_gauss[2] = 0.5555555555555555555556;
        }
        else if (nG == 4)
        {
            w_gauss[0] = 0.3478548451374538573731;
            w_gauss[1] = 0.6521451548625461426269;
            w_gauss[2] = 0.6521451548625461426269;
            w_gauss[3] = 0.3478548451374538573731;
        }
        else if (nG == 5)
        {
            w_gauss[0] = 0.2369268850561890875143;
            w_gauss[1] = 0.4786286704993664680413;
            w_gauss[2] = 0.568888888888888888889;
            w_gauss[3] = 0.4786286704993664680413;
            w_gauss[4] = 0.2369268850561890875143;
        }
        else if (nG == 6)
        {
            w_gauss[0] = 0.1713244923791703450403;
            w_gauss[1] = 0.36076157304813860757;
            w_gauss[2] = 0.46791393457269104739;
            w_gauss[3] = 0.46791393457269104739;
            w_gauss[4] = 0.36076157304813860757;
            w_gauss[5] = 0.1713244923791703450403;
        }
        else if (nG == 7)
        {
            w_gauss[0] = 0.1294849661688696932706;
            w_gauss[1] = 0.279705391489276667901;
            w_gauss[2] = 0.38183005050511894495;
            w_gauss[3] = 0.4179591836734693877552;
            w_gauss[4] = 0.38183005050511894495;
            w_gauss[5] = 0.279705391489276667901;
            w_gauss[6] = 0.1294849661688696932706;
        }
        else if (nG == 8)
        {
            w_gauss[0] = 0.1012285362903762591525;
            w_gauss[1] = 0.222381034453374470544;
            w_gauss[2] = 0.3137066458778872873389;
            w_gauss[3] = 0.3626837833783619829652;
            w_gauss[4] = 0.3626837833783619829652;
            w_gauss[5] = 0.3137066458778872873389;
            w_gauss[6] = 0.222381034453374470544;
            w_gauss[7] = 0.1012285362903762591525;
        }
        else if (nG == 9)
        {
            w_gauss[0] = 0.0812743883615744119723;
            w_gauss[1] = 0.180648160694857404058;
            w_gauss[2] = 0.2606106964029354623187;
            w_gauss[3] = 0.312347077040002840069;
            w_gauss[4] = 0.3302393550012597631654;
            w_gauss[5] = 0.312347077040002840069;
            w_gauss[6] = 0.2606106964029354623187;
            w_gauss[7] = 0.180648160694857404058;
            w_gauss[8] = 0.0812743883615744119723;
        }
        else if (nG == 10)
        {
            w_gauss[0] = 0.0666713443086881375936;
            w_gauss[1] = 0.1494513491505805931458;
            w_gauss[2] = 0.2190863625159820439955;
            w_gauss[3] = 0.2692667193099963550912;
            w_gauss[4] = 0.295524224714752870174;
            w_gauss[5] = 0.295524224714752870174;
            w_gauss[6] = 0.2692667193099963550912;
            w_gauss[7] = 0.2190863625159820439955;
            w_gauss[8] = 0.1494513491505805931458;
            w_gauss[9] = 0.0666713443086881375936;
        }
        else if (nG == 11)
        {
            w_gauss[0] = 0.0556685671161736664828;
            w_gauss[1] = 0.125580369464904624635;
            w_gauss[2] = 0.1862902109277342514261;
            w_gauss[3] = 0.233193764591990479919;
            w_gauss[4] = 0.262804544510246662181;
            w_gauss[5] = 0.2729250867779006307145;
            w_gauss[6] = 0.262804544510246662181;
            w_gauss[7] = 0.233193764591990479919;
            w_gauss[8] = 0.1862902109277342514261;
            w_gauss[9] = 0.125580369464904624635;
            w_gauss[10] = 0.0556685671161736664828;
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "The number of layers is bigger that the ones considered in the Gauss Quadrature implementation", "");
        }

        return w_gauss;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives the inverse of a 2x2 matrix
     * @param rMatrixOrig, matrix to invert
     * @return InvMatrix, inverted matrix
     */
    static inline void InvMat2x2(
            const boost::numeric::ublas::bounded_matrix<double, 2, 2 > & rMatrixOrig,
            boost::numeric::ublas::bounded_matrix<double, 2, 2 > & InvMatrix
            )
    {
        double det = rMatrixOrig(0, 0) * rMatrixOrig(1, 1) - rMatrixOrig(0, 1) * rMatrixOrig(1, 0);

//        KRATOS_WATCH(rMatrixOrig);
//        KRATOS_WATCH(det);

        if (det < 1.0e-18)
        {
            KRATOS_THROW_ERROR( std::invalid_argument," Determinant of the matrix is 0 or negative!!!, det = ", det);
        }

        /* Compute inverse of the Matrix */
        InvMatrix(0, 0) =   rMatrixOrig(1, 1) / det;
        InvMatrix(0, 1) = - rMatrixOrig(0, 1) / det;
        InvMatrix(1, 0) = - rMatrixOrig(1, 0) / det;
        InvMatrix(1, 1) =   rMatrixOrig(0, 0) / det;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives the determinant of a 3x3 matrix
     * @param rMatrixOrig, matrix to calculate determinant
     * @return det, the matrix determinant
     */
    static inline void DetMat3x3(
            const Matrix& rMatrixOrig,
            double& det
            )
    {
        /* Compute determinant of the Jacobian */
        det = rMatrixOrig(0, 0) * rMatrixOrig(1, 1) * rMatrixOrig(2, 2) \
                + rMatrixOrig(1, 0) * rMatrixOrig(2, 1) * rMatrixOrig(0, 2)\
                + rMatrixOrig(0, 1) * rMatrixOrig(1, 2) * rMatrixOrig(2, 0)\
                - rMatrixOrig(2, 0) * rMatrixOrig(1, 1) * rMatrixOrig(0, 2)\
                - rMatrixOrig(2, 1) * rMatrixOrig(1, 2) * rMatrixOrig(0, 0)\
                - rMatrixOrig(1, 0) * rMatrixOrig(0, 1) * rMatrixOrig(2,2);

        if (det < 1.0e-18)
        {
            KRATOS_THROW_ERROR( std::invalid_argument," Determinant of the matrix is 0 or negative!!!, det = ", det);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives the inverse of a 3x3 matrix
     * @param rMatrixOrig, matrix to invert
     * @return InvMatrix, inverted matrix
     */
    static inline void InvMat3x3(
            const Matrix& rMatrixOrig,
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InvMatrix
            )
    {
        /* Compute determinant of the Jacobian */
        double det;
        DetMat3x3(rMatrixOrig, det);

        /* Compute inverse of the Matrix */
        InvMatrix(0, 0) =   (rMatrixOrig(1, 1) * rMatrixOrig(2, 2) - rMatrixOrig(1, 2) * rMatrixOrig(2, 1)) / det;
        InvMatrix(1, 0) = - (rMatrixOrig(1, 0) * rMatrixOrig(2, 2) - rMatrixOrig(2, 0) * rMatrixOrig(1, 2)) / det;
        InvMatrix(2, 0) =   (rMatrixOrig(1, 0) * rMatrixOrig(2, 1) - rMatrixOrig(1, 1) * rMatrixOrig(2, 0)) / det;

        InvMatrix(0, 1) = - (rMatrixOrig(0, 1) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 1)) / det;
        InvMatrix(1, 1) =   (rMatrixOrig(0, 0) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 0)) / det;
        InvMatrix(2, 1) = - (rMatrixOrig(0, 0) * rMatrixOrig(2, 1) - rMatrixOrig(0, 1) * rMatrixOrig(2, 0)) / det;

        InvMatrix(0, 2) =   (rMatrixOrig(0, 1) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 1)) / det;
        InvMatrix(1, 2) = - (rMatrixOrig(0, 0) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 0)) / det;
        InvMatrix(2, 2) =   (rMatrixOrig(0, 0) * rMatrixOrig(1, 1) - rMatrixOrig(1, 0) * rMatrixOrig(0, 1)) / det;

//        boost::numeric::ublas::bounded_matrix<double, 3, 3 > aux;
//        noalias(aux) = prod(InvMatrix, rMatrixOrig);
//        KRATOS_WATCH(aux);

    }

    /***********************************************************************************/
    /***********************************************************************************/

    static inline void InvMat3x3(
            const Matrix& rMatrixOrig,
            double& det,
            boost::numeric::ublas::bounded_matrix<double, 3, 3 > & InvMatrix
            )
    {
        /* Compute determinant of the Jacobian */
        DetMat3x3(rMatrixOrig, det);

        /* Compute inverse of the Matrix */
        InvMatrix(0, 0) = (rMatrixOrig(1, 1) * rMatrixOrig(2, 2) - rMatrixOrig(1, 2) * rMatrixOrig(2, 1)) / det;
        InvMatrix(1, 0) = - (rMatrixOrig(1, 0) * rMatrixOrig(2, 2) - rMatrixOrig(2, 0) * rMatrixOrig(1, 2)) / det;
        InvMatrix(2, 0) = (rMatrixOrig(1, 0) * rMatrixOrig(2, 1) - rMatrixOrig(1, 1) * rMatrixOrig(2, 0)) / det;
        InvMatrix(0, 1) = - (rMatrixOrig(0, 1) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 1)) / det;
        InvMatrix(1, 1) = (rMatrixOrig(0, 0) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 0)) / det;
        InvMatrix(2, 1) = - (rMatrixOrig(0, 0) * rMatrixOrig(2, 1) - rMatrixOrig(0, 1) * rMatrixOrig(2, 0)) / det;
        InvMatrix(0, 2) = (rMatrixOrig(0, 1) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 1)) / det;
        InvMatrix(1, 2) = - (rMatrixOrig(0, 0) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 0)) / det;
        InvMatrix(2, 2) = (rMatrixOrig(0, 0) * rMatrixOrig(1, 1) - rMatrixOrig(1, 0) * rMatrixOrig(0, 1)) / det;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives one column of the inverse of a 3x3 matrix
     * @param rMatrixOrig, matrix to invert
     * @param column, index of the column
     * @return InvMatrix, inverted matrix
     */
    static inline void InvMat3x3Col(
            const Matrix& rMatrixOrig,
            const int& column,
            boost::numeric::ublas::bounded_matrix<double, 3, 1 > & InvMatrix
            )
    {
        /* Compute determinant of the Jacobian */
        double det = rMatrixOrig(0, 0) * rMatrixOrig(1, 1) * rMatrixOrig(2, 2) \
                + rMatrixOrig(1, 0) * rMatrixOrig(2, 1) * rMatrixOrig(0, 2)\
                + rMatrixOrig(0, 1) * rMatrixOrig(1, 2) * rMatrixOrig(2, 0)\
                - rMatrixOrig(2, 0) * rMatrixOrig(1, 1) * rMatrixOrig(0, 2)\
                - rMatrixOrig(2, 1) * rMatrixOrig(1, 2) * rMatrixOrig(0, 0)\
                - rMatrixOrig(1, 0) * rMatrixOrig(0, 1) * rMatrixOrig(2,2);

        if (det < 1.0e-18)
        {
            KRATOS_THROW_ERROR( std::invalid_argument," Determinant of the matrix is 0 or negative!!!, det = ", det);
        }

        /* Compute inverse of the Matrix (just one column) */
        if (column == 1)
        {
            InvMatrix(0, 0) = (rMatrixOrig(1, 1) * rMatrixOrig(2, 2) - rMatrixOrig(1, 2) * rMatrixOrig(2, 1)) / det;
            InvMatrix(1, 0) = - (rMatrixOrig(1, 0) * rMatrixOrig(2, 2) - rMatrixOrig(2, 0) * rMatrixOrig(1, 2)) / det;
            InvMatrix(2, 0) = (rMatrixOrig(1, 0) * rMatrixOrig(2, 1) - rMatrixOrig(1, 1) * rMatrixOrig(2, 0)) / det;
        }
        else if (column == 2)
        {
            InvMatrix(0, 0) = - (rMatrixOrig(0, 1) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 1)) / det;
            InvMatrix(1, 0) = (rMatrixOrig(0, 0) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 0)) / det;
            InvMatrix(2, 0) = - (rMatrixOrig(0, 0) * rMatrixOrig(2, 1) - rMatrixOrig(0, 1) * rMatrixOrig(2, 0)) / det;
        }
        else if (column == 3)
        {
            InvMatrix(0, 0) = (rMatrixOrig(0, 1) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 1)) / det;
            InvMatrix(1, 0) = - (rMatrixOrig(0, 0) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 0)) / det;
            InvMatrix(2, 0) = (rMatrixOrig(0, 0) * rMatrixOrig(1, 1) - rMatrixOrig(1, 0) * rMatrixOrig(0, 1)) / det;
        }
        else
        {
            KRATOS_THROW_ERROR( std::invalid_argument," The column index must be between 1 and 3, column = ", column);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    static inline void InvMat3x3Col(
            const Matrix& rMatrixOrig,
            const int& column,
            const double& det,
            boost::numeric::ublas::bounded_matrix<double, 3, 1 > & InvMatrix
            )
    {
        if (det < 1.0e-18)
        {
            KRATOS_THROW_ERROR( std::invalid_argument," Determinant of the matrix is 0 or negative!!!, det = ", det);
        }

        /* Compute inverse of the Matrix (just one column) */
        if (column == 1)
        {
            InvMatrix(0, 0) =   (rMatrixOrig(1, 1) * rMatrixOrig(2, 2) - rMatrixOrig(1, 2) * rMatrixOrig(2, 1)) / det;
            InvMatrix(1, 0) = - (rMatrixOrig(1, 0) * rMatrixOrig(2, 2) - rMatrixOrig(2, 0) * rMatrixOrig(1, 2)) / det;
            InvMatrix(2, 0) =   (rMatrixOrig(1, 0) * rMatrixOrig(2, 1) - rMatrixOrig(1, 1) * rMatrixOrig(2, 0)) / det;
        }
        else if (column == 2)
        {
            InvMatrix(0, 0) = - (rMatrixOrig(0, 1) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 1)) / det;
            InvMatrix(1, 0) =   (rMatrixOrig(0, 0) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 0)) / det;
            InvMatrix(2, 0) = - (rMatrixOrig(0, 0) * rMatrixOrig(2, 1) - rMatrixOrig(0, 1) * rMatrixOrig(2, 0)) / det;
        }
        else if (column == 3)
        {
            InvMatrix(0, 0) =   (rMatrixOrig(0, 1) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 1)) / det;
            InvMatrix(1, 0) = - (rMatrixOrig(0, 0) * rMatrixOrig(1, 2) - rMatrixOrig(0, 2) * rMatrixOrig(1, 0)) / det;
            InvMatrix(2, 0) =   (rMatrixOrig(0, 0) * rMatrixOrig(1, 1) - rMatrixOrig(1, 0) * rMatrixOrig(0, 1)) / det;
        }
        else
        {
            KRATOS_THROW_ERROR( std::invalid_argument," The column index must be between 1 and 3, column = ", column);
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives the in plane components of the inverse of a 3x3 matrix
     * @param rMatrixOrig, matrix to invert
     * @return InvMatrix, inverted matrix
     */
    static inline void InvMat3x3Plane(
            const Matrix& rMatrixOrig,
            boost::numeric::ublas::bounded_matrix<double, 2, 2 > & InvMatrix
            )
    {
        /* Compute determinant of the Jacobian */
        double det;
        DetMat3x3(rMatrixOrig, det);

        /* Compute inverse of the Matrix */
        InvMatrix(0, 0) = (rMatrixOrig(1, 1) * rMatrixOrig(2, 2) - rMatrixOrig(1, 2) * rMatrixOrig(2, 1)) / det;
        InvMatrix(1, 0) = - (rMatrixOrig(1, 0) * rMatrixOrig(2, 2) - rMatrixOrig(2, 0) * rMatrixOrig(1, 2)) / det;
        InvMatrix(0, 1) = - (rMatrixOrig(0, 1) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 1)) / det;
        InvMatrix(1, 1) = (rMatrixOrig(0, 0) * rMatrixOrig(2, 2) - rMatrixOrig(0, 2) * rMatrixOrig(2, 0)) / det;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Gives the eigenvalues of a symetric 3x3 matrix
     * @param rMatrixOrig , the matrix to decompose
     * @return EigenValuesVector, the vector with the eigenvalues
     */
    static inline void Eigenvalues3x3(
            const Matrix& rMatrixOrig,
            Vector& EigenValuesVector
            )
    {
        // Given a real symmetric 3x3 matrix A, compute the eigenvalues
        // Look: http://dl.acm.org/citation.cfm?doid=355578.366316

        double p1 = pow(rMatrixOrig(0, 1), 2.0) + pow(rMatrixOrig(0, 2), 2.0)+ pow(rMatrixOrig(1, 2), 2.0);
        if (p1 == 0.0)
        {
            EigenValuesVector[0] = rMatrixOrig(0, 0);
            EigenValuesVector[1] = rMatrixOrig(1, 1);
            EigenValuesVector[2] = rMatrixOrig(2, 2);
        }
        else
        {
            // Calculate eigenvalues
            double q = (rMatrixOrig(0, 0) + rMatrixOrig(1, 1) + rMatrixOrig(2, 2)) / 3.0;
            double p2 = pow((rMatrixOrig(0, 0) - q), 2.0) + pow((rMatrixOrig(1, 1) - q), 2.0) + pow((rMatrixOrig(2, 2) - q), 2.0) + 2 *p1;
            double p = sqrt(p2 / 6.0);

            boost::numeric::ublas::bounded_matrix<double, 3, 3 > aux_Matrix = ZeroMatrix(3, 3);

            noalias(aux_Matrix) = rMatrixOrig;
            aux_Matrix(0, 0) -= q;
            aux_Matrix(1, 1) -= q;
            aux_Matrix(2, 2) -= q;

            double r = aux_Matrix(0, 0) * aux_Matrix(1, 1) * aux_Matrix(2, 2) + aux_Matrix(1, 0) * aux_Matrix(2, 1) * aux_Matrix(0, 2) + aux_Matrix(0, 1) * aux_Matrix(1, 2) * aux_Matrix(2, 0)\
                    - aux_Matrix(2, 0) * aux_Matrix(1, 1) * aux_Matrix(0, 2) - aux_Matrix(2, 1) * aux_Matrix(1, 2) * aux_Matrix(0, 0) - aux_Matrix(1, 0) * aux_Matrix(0, 1) * aux_Matrix(2,2);

            r /= 2.0;

            double phi;
            if (r < 1.0)
            {
                phi = PI / 3.0;
            }
            else if (r > 1.0)
            {
                phi = 0.0;
            }
            else
            {
                phi = acos(r) / 3.0;
            }

            EigenValuesVector[0] = q + 2 * p * cos(phi);
            EigenValuesVector[2] = q + 2 * p * cos(phi + 2 * PI / 3.0);
            EigenValuesVector[1] = 3 * q - EigenValuesVector[0] - EigenValuesVector[2];
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

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
            q = -0.5 * (b + sqrt(disc));
        }
        else
        {
            q = -0.5 * (b - sqrt(disc));
        }

        solution[0] = q / a;
        solution[1] = c / q;

        return true;
    }
    
    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Calculates the solutions for a given cubic polynomial equation 0 = a*x^3 + b*x^2 + c*x + d
     * @param a coefficient
     * @param b coefficient
     * @param c coefficient
     * @param d coefficient
     * @param ZeroTol number treated as zero
     * @return Vector of solutions
     * WARNING only valid cubic (not quadratic, not linear, not constant) equations with three real (not complex) solutions
     */

    static inline bool CardanoFormula(
            const double& a,
            const double& b,
            const double& c,
            const double& d,
            Vector& solution
            )
    {
        solution.resize(3, false);
        noalias(solution)= ZeroVector(3);

        if(a == 0)
        {
            std::cout << "This is not a cubic equation: CardanoFormula" << std::endl;
            return false;
        }

        double p = (3.0 * a * c - b * b) / (3.0 * a * a);
        double q = 2.0 * b * b * b / (27.0 * a * a * a)-b * c/(3.0 * a * a)+d / a;
        double disc = p * p * p / 27.0 + q * q / 4.0;

        if(disc > 0)
        {
            return false;
        }

        if(disc == 0)
        {
            if( a == 0 )
            {
                return false;
            }

            solution(0) = pow(q / 2.0, 1.0 / 3.0) - b / (3 * a);
            solution(1) = pow(q / 2.0, 1.0 / 3.0) - b / (3 * a);
            solution(2) = pow(-4.0 * q, 1.0 / 3.0) - b / (3 * a);

            return true;
        }

        solution(0) = -sqrt(-4.0 / 3.0 * p) * cos(1.0 / 3.0 * acos(-q / 2.0 * sqrt(-27.0 / (p * p * p))) + PI / 3.0) - b / ( 3 * a);
        solution(1) =  sqrt(-4.0 / 3.0 * p) * cos(1.0 / 3.0 * acos(-q / 2.0 * sqrt(-27.0 / (p * p * p)))) - b / (3 * a);
        solution(2) = -sqrt(-4.0 / 3.0 * p) * cos(1.0 / 3.0 * acos(-q/2.0 * sqrt(-27.0 / (p * p * p))) - PI / 3.0) - b / (3 * a);

        return true;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Calculates Eigenvalues of given square matrix A.
     * The QR Algorithm with shifts is used
     * @param A the given square matrix the eigenvalues are to be calculated.
     * @param crit convergence criteria
     * @param zero number treated as zero
     * @return Vector of eigenvalues
     * WARNING only valid for 2*2 and 3*3 Matrices yet
     */

    static inline Vector EigenValues(
                const Matrix& A,
                double crit,
                double zero
                )
    {
        int dim= A.size1();

        Matrix Convergence(2, dim);

        double delta;
        double abs;
        Vector Result = ZeroVector(dim);
        Matrix HelpA = ZeroMatrix(dim, dim);
        Matrix HelpQ = ZeroMatrix(dim, dim);
        Matrix HelpR = ZeroMatrix(dim, dim);
        HelpA = A;

        bool is_converged = false;
        while(!(is_converged))
        {
            double shift= HelpA((dim-1),(dim-1));
            for(int i = 0; i<dim; i++)
            {
                HelpA(i, i) = HelpA(i, i)- shift;
            }

            StructuralMechanicsMathUtilities::QRFactorization(HelpA, HelpQ, HelpR);

            HelpA = ZeroMatrix(dim, dim);

            for(int i = 0; i < dim; i++)
            {
                HelpA(i, i) += shift;
                for(int j = 0; j < dim; j++)
                {
                    for(int k = 0; k < dim; k++)
                    {
                        HelpA(i, j) += HelpR(i, k) * HelpQ(k, j);
                    }
                }
            }

            delta = 0.0;
            abs = 0.0;

            for(int i = 0; i < dim; i++)
            {
                Convergence(0, i) = Convergence(1, i);
                Convergence(1, i) = HelpA(i, i);
                delta += (Convergence(1, i)-Convergence(0, i))*(Convergence(1, i)-Convergence(0, i));
                abs += (Convergence(1, i))*(Convergence(1, i));
            }

            delta = sqrt(delta);
            abs = sqrt(abs);

            if(abs < zero)
            {
                abs = 1.0;
            }

            if(delta < zero || (delta / abs) < crit)
            {
                is_converged = true;
            }

        }

        for(int i = 0; i < dim; i++)
        {
            Result(i)= HelpA(i, i);

            if(fabs(Result(i)) < zero)
            {
                Result(i) = 0.0;
            }
        }

        return Result;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Calculates the QR Factorization of given square matrix A=QR.
     * The Factorization is performed using the householder algorithm
     * @param A the given square matrix the factorization is to be calculated.
     * @param Q the result matrix Q
     * @param R the result matrix R
     */

    static inline void QRFactorization(
                const Matrix& A,
                Matrix& Q,
                Matrix& R
                )
    {
        // QR Factorization with Householder-Algo
        int dim = A.size1();

        Vector y(dim);
        Vector w(dim);
        R.resize(dim, dim, false);
        R = ZeroMatrix(dim, dim);
        Q.resize(dim, dim,false);
        Q = ZeroMatrix(dim, dim);
        Matrix Help= A;
        Matrix unity= ZeroMatrix(dim, dim);

        for(int j = 0; j < dim; j++)
        {
            unity(j, j) = 1.0;
        }

        std::vector<Matrix> HelpQ(dim-1);
        std::vector<Matrix> HelpR(dim-1);

        for(int i = 0; i < dim-1; i++)
        {
            HelpQ[i].resize(dim,dim,false);
            HelpR[i].resize(dim,dim,false);
            noalias(HelpQ[i])= unity;
            noalias(HelpR[i])= ZeroMatrix(dim, dim);
        }

        for(int iteration = 0; iteration < dim-1; iteration++)
        {
            // Vector y
            for(int i = iteration; i < dim; i++)
            {
                y(i) = Help(i, iteration);
            }

            // Helpvalue l
            double normy = 0.0;

            for(int i = iteration; i < dim; i++)
            {
                normy += y(i) * y(i);
            }

            normy = sqrt(normy);

            double l = sqrt((normy * (normy + fabs(y(iteration)))) / 2.0);
            double k = 0.0;

            if(y[iteration] != 0)
            {
                k = - y(iteration) / fabs(y(iteration)) * normy;
            }
            else
            {
                k = -normy;
            }

            for(int i=iteration; i<dim; i++)
            {
                double e = 0.0;

                if(i==iteration)
                {
                    e = 1.0;
                }

                w(i) = 1.0 / (2.0 * l) * (y(i) - k * e);
            }

            for(int i = iteration; i < dim; i++)
            {
                for(int j = iteration; j < dim; j++)
                {
                    HelpQ[iteration](i, j)= unity(i, j)- 2.0 * w(i) * w(j);
                }
            }


            for(int i = iteration; i < dim; i++)
            {
                for(int j = iteration; j < dim; j++)
                {
                    for(int k = iteration; k < dim; k++)
                    {
                        HelpR[iteration](i, j) += HelpQ[iteration](i, k) * Help(k, j);
                    }
                }
            }

            Help= HelpR[iteration];

        }

        // Assembling R
        for(int k = 0; k < dim-1; k++)
        {
            for(int i = k; i < dim; i++)
            {
                for(int j = k; j < dim; j++)
                {
                    R(i, j) = HelpR[k](i, j);
                }
            }

        }


        for(int k = 1; k < dim-1; k++)
        {
            for(int i = 0; i < dim; i++)
            {
                for(int j = 0; j < dim; j++)
                {
                    for(int l = 0; l < dim; l++)
                    {
                        Q(i,j) += HelpQ[(k-1)](i, l) * HelpQ[k](l, j);
                    }
                }
            }
            noalias(HelpQ[k]) = Q;
        }
        if(dim-1 == 1)
        {
            noalias(Q) = HelpQ[0];
        }

    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Calculates the eigenvectors and eigenvalues of given symmetric matrix A.
     * The eigenvectors and eigenvalues are calculated using the iterative Gauss-Seidel-method
     * @param A: The given symmetric matrix the eigenvectors are to be calculated.
     * @return vectors: The result matrix (will be overwritten with the eigenvectors)
     * @return lambda: The result diagonal matrix with the eigenvalues
     * @param zero_tolerance: The largest value considered to be zero
     * @param max_iterations: Maximum number of iterations
     * WARNING: Matrix A will be overwritten and has to be symmetric
     */

    static inline void EigenVectors(
            const Matrix& A,
            Matrix& vectors,
            Matrix& lambda,
            double zero_tolerance = 1e-9,
            int max_iterations = 10
            )
    {

        Matrix Help = A;

        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                Help(i, j) = Help(i, j);
            }
        }

        vectors.resize(Help.size1(), Help.size2(), false);
        lambda.resize(Help.size1(), Help.size2(), false);
        lambda = ZeroMatrix(Help.size1());
        Matrix HelpDummy(Help.size1(),Help.size2());
        bool is_converged = false;
        Matrix unity = ZeroMatrix(Help.size1(), Help.size2());

        for(unsigned int i = 0; i < Help.size1(); i++)
        {
            unity(i,i) = 1.0;
        }

        Matrix V = unity;
        Matrix VDummy(Help.size1(),Help.size2());
        Matrix Rotation(Help.size1(),Help.size2());

        for(int iterations = 0; iterations < max_iterations; iterations++)
        {
            is_converged = true;

            double a = 0.0;
            unsigned int index1 = 0;
            unsigned int index2 = 1;

            for(unsigned int i = 0; i < Help.size1(); i++)
            {
                for(unsigned int j = (i + 1); j < Help.size2(); j++)
                {
                    if((fabs(Help(i, j)) > a ) && (fabs(Help(i, j)) > zero_tolerance))
                    {
                        a = fabs(Help(i,j));
                        index1 = i;
                        index2 = j;
                        is_converged = false;
                    }
                }
            }

            if(is_converged)
            {
                break;
            }

            // Calculation of Rotation angle
            double gamma = (Help(index2, index2)-Help(index1, index1)) / (2 * Help(index1, index2));
            double u = 1.0;

            if(fabs(gamma) > zero_tolerance && fabs(gamma)< (1/zero_tolerance))
            {
                u = gamma / fabs(gamma) * 1.0 / (fabs(gamma) + sqrt(1.0 + gamma * gamma));
            }
            else
            {
                if  (fabs(gamma) >= (1.0/zero_tolerance))
                {
                    u = 0.5 / gamma;
                }
            }

            double c = 1.0 / (sqrt(1.0 + u * u));
            double s = c * u;
            double teta = s / (1.0 + c);

            // Rotation of the Matrix
            HelpDummy = Help;
            HelpDummy(index2, index2) = Help(index2,index2) + u * Help(index1, index2);
            HelpDummy(index1, index1) = Help(index1,index1) - u * Help(index1, index2);
            HelpDummy(index1, index2) = 0.0;
            HelpDummy(index2, index1) = 0.0;

            for(unsigned int i=0; i<Help.size1(); i++)
            {
                if((i!= index1) && (i!= index2))
                {
                    HelpDummy(index2, i) = Help(index2, i) + s * (Help(index1, i)- teta * Help(index2, i));
                    HelpDummy(i, index2) = Help(index2, i) + s * (Help(index1, i)- teta * Help(index2, i));
                    HelpDummy(index1, i) = Help(index1, i) - s * (Help(index2, i) + teta * Help(index1, i));
                    HelpDummy(i, index1) = Help(index1, i) - s * (Help(index2, i) + teta * Help(index1, i));
                }
            }

            Help = HelpDummy;

            // Calculation of the eigenvectors V
            Rotation = unity;
            Rotation(index2, index1) = -s;
            Rotation(index1, index2) = s;
            Rotation(index1, index1) = c;
            Rotation(index2, index2) = c;

            VDummy = ZeroMatrix(Help.size1(), Help.size2());

            for(unsigned int i = 0; i < Help.size1(); i++)
            {
                for(unsigned int j = 0; j < Help.size1(); j++)
                {
                    for(unsigned int k = 0; k < Help.size1(); k++)
                    {
                        VDummy(i, j) += V(i, k) * Rotation(k, j);
                    }
                }
            }
            V = VDummy;
        }

        if(!(is_converged))
        {
            std::cout<<"########################################################"<<std::endl;
            std::cout<<"Max_Iterations exceed in Jacobi-Seidel-Iteration (eigenvectors)"<<std::endl;
            std::cout<<"########################################################"<<std::endl;
        }

        for(unsigned int i = 0; i < Help.size1(); i++)
        {
            for(unsigned int j = 0; j < Help.size1(); j++)
            {
                vectors(i, j) = V(j, i);
            }
        }

        for(unsigned int i = 0; i < Help.size1(); i++)
        {
            lambda(i, i) = Help(i, i);
        }

        return;
    }

    /**
     * Calculates the eigenvectors and eigenvalues of given matrix A.
     * The eigenvectors and eigenvalues are calculated using the iterative JACOBI-method
     * @param A the given matrix the eigenvectors are to be calculated.
     * @return V the result matrix (will be overwritten with the eigenvectors)
     * @param error_tolerance the desired accuracy for the convergence check
     * @param zero_tolerance the largest value considered to be zero
     * WARNING: Matrix A will be overwritten
     */
    static inline void EigenVectors(
            Matrix& A,
            Matrix& V,
            double& error_tolerance,
            double& zero_tolerance
            )
    {
        // Initial error
        double error = 1.0;
        int n = A.size2();
        // Setting V to identity matrix
        V = IdentityMatrix( V.size1() );
        // Calculation loop (as long as there is no convergence)
        // WARNING: Iteration never exceeds
        while( error > error_tolerance )
        {
            for( int i = 0; i < n; i++ )
            {
                for( int j = i + 1; j<n; j++ )
                {
                    double theta = 0.0;
                    if( MathUtils<double>::Abs( A(i,j) ) >= zero_tolerance )
                    {
                        if( MathUtils<double>::Abs( A(i, i)-A(j, j) ) > 0.0 )
                        {
                            theta = 0.5 * atan(2 * A(i, j)/(A(i, i)-A(j, j)));
                        }
                        else theta = 0.25 * PI;
                    }
                    Matrix T = IdentityMatrix(n);

                    T(i, i) = cos(theta);
                    T(i, j) = -sin(theta);
                    T(j, i) = -T(i, j);
                    T(j, j) = T(i, i);

                    noalias(A) = prod(A, T);
                    Matrix TT = trans(T);
                    noalias(A) = prod(TT, A);
                    noalias(V) = prod(V, T);
                }
            }
            double sTot = 0.0;
            double sDiag = 0.0;
            for( unsigned int i = 0; i < A.size1(); i++ )
            {
                for( unsigned int j = 0; j < A.size2(); j++ )
                {
                    sTot += MathUtils<double>::Abs(A(i, j));
                }
                sDiag += MathUtils<double>::Abs(A(i, i));
            }
            error = (sTot - sDiag) / sDiag;
        }
        // Sorting eigenvalues
        int maxIndex = 0;
        double maxEv = A(0, 0);
        for( unsigned int i = 0; i < A.size1(); i++ )
        {
            for( unsigned int j = i; j < A.size1(); j++ )
            {
                // Searching current maximum
                if( A(j,j) > maxEv )
                {
                    maxIndex = j;
                    maxEv = A(j, j);
                }
                // Swapping eigenvalue matrix
                double dummy = A(i, i);
                A(i, i) = A(maxIndex, maxIndex);
                A(maxIndex, maxIndex) = dummy;
                // Swapping eigenvector matrix
                for( unsigned int k = 0; k < A.size2(); k++ )
                {
                    dummy = V(k, i);
                    V(k, i) = V(k, maxIndex);
                    V(k, maxIndex) = dummy;
                }
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    
    /**
     * Converts a strain vector into a matrix. Strains are assumed to be stored in the following way:
     * \f$ [ e11, e22, e33, 2*e12, 2*e23, 2*e13 ] \f$ for 3D case and
     * \f$ [ e11, e22, 2*e12 ] \f$ for 2D case.
     * Hence the deviatoric components of the strain vector are divided by 2 while they are stored into the matrix
     * @param Strains the given strain vector
     * @return the corresponding strain tensor in matrix form
     */
    static inline Matrix StrainVectorToTensor( const Vector& Strains )
    {
        KRATOS_TRY;

        Matrix StrainTensor;

        if (Strains.size() == 3)
        {
            StrainTensor.resize(2, 2, false);
            StrainTensor(0, 0) = Strains[0];
            StrainTensor(0, 1) = 0.5 * Strains[2];
            StrainTensor(1, 0) = 0.5 * Strains[2];
            StrainTensor(1, 1) = Strains[1];
        }
        else if (Strains.size() == 6)
        {
            StrainTensor.resize(3, 3, false);
            StrainTensor(0, 0) = Strains[0];
            StrainTensor(0, 1) = 0.5 * Strains[3];
            StrainTensor(0, 2) = 0.5 * Strains[5];
            StrainTensor(1, 0) = 0.5 * Strains[3];
            StrainTensor(1, 1) = Strains[1];
            StrainTensor(1, 2) = 0.5 * Strains[4];
            StrainTensor(2, 0) = 0.5 * Strains[5];
            StrainTensor(2, 1) = 0.5 * Strains[4];
            StrainTensor(2, 2) = Strains[2];
        }

        return StrainTensor;
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    static inline Vector TensorToStrainVector( const Matrix& Tensor )
    {
        KRATOS_TRY;
        Vector StrainVector;

        if (Tensor.size1()==2)
        {
            StrainVector.resize(3, false);
            noalias(StrainVector) = ZeroVector(3);
            StrainVector(0) = Tensor(0, 0);
            StrainVector(1) = Tensor(1, 1);
            StrainVector(2) = 2.0 * Tensor(0, 1);
        }
        else if (Tensor.size1()==3)
        {
            StrainVector.resize(6, false);
            noalias(StrainVector) = ZeroVector(6);
            StrainVector(0) = Tensor(0, 0);
            StrainVector(1) = Tensor(1, 1);
            StrainVector(2) = Tensor(2, 2);
            StrainVector(3) = 2.0 * Tensor(0, 1);
            StrainVector(4) = 2.0 * Tensor(1, 2);
            StrainVector(5) = 2.0 * Tensor(0, 2);
        }

        return StrainVector;
        KRATOS_CATCH("");
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Builds the Inverse of Matrix input
    * @param input the given Matrix
    * @param inverse inverse of the given Matrix
    */
    static int InvertMatrix(
            const Matrix& input,
            Matrix& inverse
            )
    {
        int singular = 0;
        using namespace boost::numeric::ublas;
        typedef permutation_matrix<std::size_t> pmatrix;
        Matrix A(input);
        pmatrix pm(A.size1());
        singular = lu_factorize(A,pm);
        inverse.assign( IdentityMatrix(A.size1()));
        lu_substitute(A, pm, inverse);
        return singular;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Builds the norm of a gibven second order tensor
    * @param Tensor the given second order tensor
    * @return the norm of the given tensor
    */
    static double normTensor(Matrix& Tensor)
    {
        double result = 0.0;
        for(unsigned int i = 0; i < Tensor.size1(); i++)
        {
            for(unsigned int j = 0; j < Tensor.size2(); j++)
            {
                result += Tensor(i,j)*Tensor(i,j);
            }
        }

        result = sqrt(result);

        return result;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Transforms a given 6*1 Vector to a corresponing symmetric Tensor of second order (3*3)
    * @param Vector the given vector
    * @param Tensor the symmetric second order tensor
    */
    static inline void VectorToTensor(
            const Vector& Stress,
            Matrix& Tensor
            )
    {
        if(Stress.size() == 6)
        {
            Tensor.resize(3, 3, false);
            Tensor(0, 0) = Stress(0);
            Tensor(0, 1) = Stress(3);
            Tensor(0, 2) = Stress(5);
            Tensor(1, 0) = Stress(3);
            Tensor(1, 1) = Stress(1);
            Tensor(1, 2) = Stress(4);
            Tensor(2, 0) = Stress(5);
            Tensor(2, 1) = Stress(4);
            Tensor(2, 2) = Stress(2);
        }
        if(Stress.size() == 3)
        {
            Tensor.resize(2, 2, false);
            Tensor(0, 0)= Stress(0);
            Tensor(0, 1)= Stress(2);
            Tensor(1, 0)= Stress(2);
            Tensor(1, 1)= Stress(1);
        }
        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Transforms a given symmetric Tensor of second order (3*3) to a corresponing 6*1 Vector
    * @param Tensor the given symmetric second order tensor
    * @param Vector the vector
    */
    static void TensorToVector(
            const Matrix& Tensor,
            Vector& Vector
            )
    {
        unsigned int  dim  =  Tensor.size1();
        if (dim==3)
        {
            Vector.resize(6, false);
            Vector(0)= Tensor(0, 0);
            Vector(1)= Tensor(1, 1);
            Vector(2)= Tensor(2, 2);
            Vector(3)= Tensor(0, 1);
            Vector(4)= Tensor(1, 2);
            Vector(5)= Tensor(2, 0);
        }
        else if(dim == 2)
        {
            Vector.resize(3,false);
            Vector(0) = Tensor(0, 0);
            Vector(1) = Tensor(1, 1);
            Vector(2) = Tensor(0, 1);
        }
        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    static inline void TensorToMatrix(
            Fourth_Order_Tensor& Tensor,
            Matrix& Matrix
            )
    {
        if (Tensor[0].size()== 3)
        {
            // Tensor de cuarto orden cuyos componentes correspondes a una matriz de 3x3
            if(Matrix.size1()!=6 || Matrix.size2()!=6)
            {
                Matrix.resize(6,6,false);
            }

            Matrix(0, 0) = Tensor[0][0](0, 0);
            Matrix(0, 1) = Tensor[0][0](1, 1);
            Matrix(0, 2) = Tensor[0][0](2, 2);
            Matrix(0, 3) = Tensor[0][0](0, 1);
            Matrix(0, 4) = Tensor[0][0](0, 2);
            Matrix(0, 5) = Tensor[0][0](1, 2);

            Matrix(1, 0) = Tensor[1][1](0, 0);
            Matrix(1, 1) = Tensor[1][1](1, 1);
            Matrix(1, 2) = Tensor[1][1](2, 2);
            Matrix(1, 3) = Tensor[1][1](0, 1);
            Matrix(1, 4) = Tensor[1][1](0, 2);
            Matrix(1, 5) = Tensor[1][1](1, 2);

            Matrix(2, 0) = Tensor[2][2](0, 0);
            Matrix(2, 1) = Tensor[2][2](1, 1);
            Matrix(2, 2) = Tensor[2][2](2, 2);
            Matrix(2, 3) = Tensor[2][2](0, 1);
            Matrix(2, 4) = Tensor[2][2](0, 2);
            Matrix(2, 5) = Tensor[2][2](1, 2);

            Matrix(3, 0) = Tensor[0][1](0, 0);
            Matrix(3, 1) = Tensor[0][1](1, 1);
            Matrix(3, 2) = Tensor[0][1](2, 2);
            Matrix(3, 3) = Tensor[0][1](0, 1);
            Matrix(3, 4) = Tensor[0][1](0, 2);
            Matrix(3, 5) = Tensor[0][1](1, 2);

            Matrix(4, 0) = Tensor[0][2](0, 0);
            Matrix(4, 1) = Tensor[0][2](1, 1);
            Matrix(4, 2) = Tensor[0][2](2, 2);
            Matrix(4, 3) = Tensor[0][2](0, 1);
            Matrix(4, 4) = Tensor[0][2](0, 2);
            Matrix(4, 5) = Tensor[0][2](1, 2);

            Matrix(5, 0) = Tensor[1][2](0, 0);
            Matrix(5, 1) = Tensor[1][2](1, 1);
            Matrix(5, 2) = Tensor[1][2](2, 2);
            Matrix(5, 3) = Tensor[1][2](0, 1);
            Matrix(5, 4) = Tensor[1][2](0, 2);
            Matrix(5, 5) = Tensor[1][2](1, 2);
        }
        else
        {
            // Tensor de cuarto orden cuyos componentes correspondes a una matriz de 2x2
            if(Matrix.size1()!=3 || Matrix.size2()!=3)
            {
                Matrix.resize(3,3,false);
            }

            Matrix(0, 0) = Tensor[0][0](0, 0);
            Matrix(0, 1) = Tensor[0][0](1, 1);
            Matrix(0, 2) = Tensor[0][0](0, 1);
            Matrix(1, 0) = Tensor[1][1](0, 0);
            Matrix(1, 1) = Tensor[1][1](1, 1);
            Matrix(1, 2) = Tensor[1][1](0, 1);
            Matrix(2, 0) = Tensor[0][1](0, 0);
            Matrix(2, 1) = Tensor[0][1](1, 1);
            Matrix(2, 2) = Tensor[0][1](0, 1);
        }
        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Transforms a given 6*6 Matrix to a corresponing 4th order tensor
    * @param Tensor the given Matrix
    * @param Vector the Tensor
    */
    static void MatrixToTensor(
            Matrix& A,
            std::vector<std::vector<Matrix> >& Tensor)
    {
        int help1 = 0;
        int help2 = 0;
        double coeff = 1.0;

        Tensor.resize(3);

        for(unsigned int i = 0; i < 3; i++)
        {
            Tensor[i].resize(3);
            for(unsigned int j = 0; j < 3; j++)
            {
                Tensor[i][j].resize(3, 3, false);
                noalias(Tensor[i][j])= ZeroMatrix(3, 3);
                for(unsigned int k = 0; k < 3; k++)
                    for(unsigned int l = 0; l < 3; l++)
                    {
                        if(i == j) help1 = i;
                        else
                        {
                            if((i == 0 && j == 1) || (i == 1 && j == 0)) help1 = 3;
                            if((i == 1 && j == 2) || (i == 2 && j == 1)) help1 = 4;
                            if((i == 2 && j == 0) || (i == 0 && j == 2)) help1 = 5;
                        }
                        if(k == l)
                        {
                            help2 = k;
                            coeff = 1.0;
                        }
                        else
                        {
                            coeff = 0.5;
                            if((k == 0 && l == 1) || (k == 1 && l == 0)) help2 = 3;
                            if((k == 1 && l == 2) || (k == 2 && l == 1)) help2 = 4;
                            if((k == 2 && l == 0) || (k == 0 && l == 2)) help2 = 5;
                        }

                        Tensor[i][j](k,l)= A(help1, help2)*coeff;
                    }
            }
        }

        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Transforms a given 6*6 Matrix to a corresponing 4th order tensor
    * @param Tensor the given Matrix
    * @param Vector the Tensor
    */
    static void MatrixToTensor(
            Matrix& A,
            array_1d<double,81>& Tensor
            )
    {
        int help1 = 0;
        int help2 = 0;
        double coeff = 1.0;
        std::fill(Tensor.begin(), Tensor.end(), 0.0);
        for(unsigned int i = 0; i < 3; i++)
        {
            for(unsigned int j = 0; j<3; j++)
            {
                for(unsigned int k = 0; k < 3; k++)
                    for(unsigned int l = 0; l < 3; l++)
                    {
                        if(i == j) help1 = i;
                        else
                        {
                            if((i == 0 && j == 1) || (i == 1 && j == 0)) help1 = 3;
                            if((i == 1 && j == 2) || (i == 2 && j == 1)) help1 = 4;
                            if((i == 2 && j == 0) || (i == 0 && j == 2)) help1 = 5;
                        }
                        if(k == l)
                        {
                            help2 = k;
                            coeff = 1.0;
                        }
                        else
                        {
                            coeff = 0.5;
                            if((k == 0 && l == 1) || (k == 1 && l == 0)) help2 = 3;
                            if((k == 1 && l == 2) || (k == 2 && l == 1)) help2 = 4;
                            if((k == 2 && l == 0) || (k == 0 && l == 2)) help2 = 5;
                        }

                        Tensor[i * 27 + j * 9 + k * 3 + l]= A(help1, help2) * coeff;
                    }
            }
        }

        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
    * Transforms a given 4th order tensor to a corresponing 6*6 Matrix
    * @param Tensor the given Tensor
    * @param Vector the Matrix
    */
    static void TensorToMatrix(
            std::vector<std::vector<Matrix> >& Tensor,
            Matrix& Matrix
            )
    {
        int help1 = 0;
        int help2 = 0;
        int help3 = 0;
        int help4 = 0;
        double coeff = 1.0;

        if(Matrix.size1()!=6 || Matrix.size2()!=6)
        {
            Matrix.resize(6,6,false);
        }

        for(unsigned int i = 0; i < 6; i++)
            for(unsigned int j = 0; j < 6; j++)
            {
                if(i < 3)
                {
                    help1 = i;
                    help2 = i;
                }
                else
                {
                    if(i == 3)
                    {
                        help1 = 0;
                        help2 = 1;
                    }
                    if(i == 4)
                    {
                        help1 = 1;
                        help2 = 2;
                    }
                    if(i == 5)
                    {
                        help1 = 2;
                        help2 = 0;
                    }
                }

                if(j < 3)
                {
                    help3 = j;
                    help4 = j;
                    coeff = 1.0;
                }
                else
                {
                    if(j == 3)
                    {
                        help3 = 0;
                        help4 = 1;
                    }
                    if(j == 4)
                    {
                        help3 = 1;
                        help4 = 2;
                    }
                    if(j == 5)
                    {
                        help3 = 2;
                        help4 = 0;
                    }
                    coeff = 2.0;
                }
                Matrix(i,j) = Tensor[help1][help2](help3,help4)*coeff;
            }
        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Transforms a given 4th order tensor to a corresponing 6*6 Matrix
     * @param Tensor the given Tensor
     * @param Vector the Matrix
     */
    static void TensorToMatrix(
            const array_1d<double,
            81>& Tensor,
            Matrix& Matrix
            )
    {
        if(Matrix.size1()!=6 || Matrix.size2()!=6)
        {
            Matrix.resize(6,6, false);
        }

        Matrix(0, 0) = Tensor[0];
        Matrix(0, 1) = Tensor[4];
        Matrix(0, 2) = Tensor[8];
        Matrix(0, 3) = 2.0 * Tensor[1];
        Matrix(0, 4) = 2.0 * Tensor[5];
        Matrix(0, 5) = 2.0 * Tensor[6];

        Matrix(1, 0) = Tensor[36];
        Matrix(1, 1) = Tensor[40];
        Matrix(1, 2) = Tensor[44];
        Matrix(1, 3) = 2.0 * Tensor[37];
        Matrix(1, 4) = 0.0 * Tensor[41];
        Matrix(1, 5) = 0.0 * Tensor[42];

        Matrix(2, 0) = Tensor[72];
        Matrix(2, 1) = Tensor[76];
        Matrix(2, 2) = Tensor[80];
        Matrix(2, 3) = 2.0 * Tensor[73];
        Matrix(2, 4) = 2.0 * Tensor[77];
        Matrix(2, 5) = 2.0 * Tensor[78];

        Matrix(3, 0) = Tensor[9];
        Matrix(3, 1) = Tensor[13];
        Matrix(3, 2) = Tensor[18];
        Matrix(3, 3) = 2.0 * Tensor[10];
        Matrix(3, 4) = 2.0 * Tensor[14];
        Matrix(3, 5) = 2.0 * Tensor[15];

        Matrix(4, 0) = Tensor[45];
        Matrix(4, 1) = Tensor[49];
        Matrix(4, 2) = Tensor[53];
        Matrix(4, 3) = 2.0 * Tensor[46];
        Matrix(4, 4) = 0.0 * Tensor[50];
        Matrix(4, 5) = 0.0 * Tensor[51];

        Matrix(5, 0) = Tensor[54];
        Matrix(5, 1) = Tensor[58];
        Matrix(5, 2) = Tensor[62];
        Matrix(5, 3) = 2.0 * Tensor[55];
        Matrix(5, 4) = 2.0 * Tensor[59];
        Matrix(5, 5) = 2.0 * Tensor[60];

        return;
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Scales a given 4th order tensor by a scalar (C = C*a)
     * @param C the given Tensor
     * @param alpha
     */
    static void ScaleFourthOrderTensor(
            Fourth_Order_Tensor& C,
            double alpha
            )
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                        C[i][j](k,l) *= alpha;
                }
            }
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/

    /**
     * Computes outer product of two 2nd order tensors (Matrix) and add to a given 4th order tensor (Result += alpha * (A \odot B))
     * @param C the given Tensor
     * @param alpha
     */
    static void OuterProductFourthOrderTensor(
            const double alpha,
            const Matrix& A,
            const Matrix& B,
            Fourth_Order_Tensor& Result
            )
    {
        for(unsigned int i = 0; i < 3; ++i)
        {
            for(unsigned int j = 0; j < 3; ++j)
            {
                for(unsigned int k = 0; k < 3; ++k)
                {
                    for(unsigned int l = 0; l < 3; ++l)
                    {
                        Result[i][j](k, l) += alpha * A(i, j) * B(k, l);
                    }
                }
            }
        }
    }

private:
};// class StructuralMechanicsMathUtilities
}
#endif /* KRATOS_STRUCTURAL_MECHANICS_MATH_UTILITIES defined */
 
