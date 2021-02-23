//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi (porting from utility made by Fabio Petracca and Peter Wilson)
//					 Philipp Bucher  (https://github.com/philbucher)
//

#if !defined(KRATOS_SHELL_UTILITIES_H_INCLUDED )
#define  KRATOS_SHELL_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "containers/array_1d.h"

namespace Kratos {

class ShellQ4_LocalCoordinateSystem; // forward declare to avoid circular includes

namespace ShellUtilities {

using SizeType = std::size_t;
using IndexType = std::size_t;

/** \brief JacobianOperator
 *
 * This class is a utility to compute at a given integration point,
 * the Jacobian, its inverse, its determinant
 * and the derivatives of the shape functions in the local
 * cartesian coordinate system.
 */
class JacobianOperator
{
public:

    JacobianOperator();

    template<class TLocalCoordinateSystem>
    void Calculate(const TLocalCoordinateSystem& CS, const Matrix& dN)
    {
        KRATOS_ERROR << "This function should not be called, this type of coordinate transformation is unknown!" << std::endl;
    }

    void Calculate(const ShellQ4_LocalCoordinateSystem& CS, const Matrix& dN);

    inline const Matrix& Jacobian()const
    {
        return mJac;
    }

    inline const Matrix& Inverse()const
    {
        return mInv;
    }

    inline const Matrix& XYDerivatives()const
    {
        return mXYDeriv;
    }

    inline double Determinant()const
    {
        return mDet;
    }

private:

    Matrix mJac;     /*!< Jacobian matrix */
    Matrix mInv;     /*!< Inverse of the Jacobian matrix */
    Matrix mXYDeriv; /*!< Shape function derivatives in cartesian coordinates */
    double mDet;     /*!< Determinant of the Jacobian matrix */
};

template<class TVec>
inline void ShapeFunc(double Xi, double Eta, TVec& rN)
{
    rN(0) = 0.25 * (1.0 - Xi) * (1.0 - Eta); // node 1
    rN(1) = 0.25 * (1.0 + Xi) * (1.0 - Eta); // node 2
    rN(2) = 0.25 * (1.0 + Xi) * (1.0 + Eta); // node 3
    rN(3) = 0.25 * (1.0 - Xi) * (1.0 + Eta); // node 4
}

template<class TMat>
inline void ShapeFunc_NaturalDerivatives(double Xi, const double Eta, TMat& rDN)
{
    rDN(0, 0) = -(1.0 - Eta) * 0.25;
    rDN(1, 0) = (1.0 - Eta) * 0.25;
    rDN(2, 0) = (1.0 + Eta) * 0.25;
    rDN(3, 0) = -(1.0 + Eta) * 0.25;

    rDN(0, 1) = -(1.0 - Xi)  * 0.25;
    rDN(1, 1) = -(1.0 + Xi)  * 0.25;
    rDN(2, 1) = (1.0 + Xi)  * 0.25;
    rDN(3, 1) = (1.0 - Xi)  * 0.25;
}

double dN_seren_dxi(const int nNode, const double Xi, const double Eta);

double dN_seren_deta(const int nNode, const double Xi, const double Eta);

void InterpToStandardGaussPoints(double& rV1, double& rV2, double& rV3);

void InterpToStandardGaussPoints(std::vector< double >& rV);

void InterpToStandardGaussPoints(std::vector< array_1d<double, 3> >& rV);

void InterpToStandardGaussPoints(std::vector< array_1d<double, 6> >& rV);

void InterpToStandardGaussPoints(std::vector< Vector >& rV);

void InterpToStandardGaussPoints(std::vector< Matrix >& rV);

bool IsOrthotropic(const Properties& rProps);

double GetThickness(const Properties& rProps);

double GetThickness(const Properties& rProps, const IndexType Index);

double GetDensity(const Properties& rProps, const IndexType Index);

double GetOrientationAngle(const Properties& rProps, const IndexType Index);

double GetOffset(const Properties& rProps);

}  // namespace Shell Utilities.

}  // namespace Kratos.

#endif // KRATOS_SHELL_UTILITIES_H_INCLUDED  defined


