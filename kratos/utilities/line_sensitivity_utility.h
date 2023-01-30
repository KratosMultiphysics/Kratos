//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#ifndef KRATOS_LINE_SENSITIVITY_UTILITY_H_INCLUDED
#define KRATOS_LINE_SENSITIVITY_UTILITY_H_INCLUDED

#include "includes/define.h"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(KRATOS_CORE) LineSensitivityUtility
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(LineSensitivityUtility);

    typedef DenseMatrix<double> MatrixType;

    typedef MatrixType JacobianType;

    typedef MatrixType ShapeFunctionsLocalGradientType;

    ///@}
    ///@name Life Cycle
    ///@{

    LineSensitivityUtility(const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De);

    ~LineSensitivityUtility() = default;

    ///@}
    ///@name Operations
    ///@{

    void CalculateSensitivity(const ShapeParameter Deriv, double& rIntegrationWeightSensitivity) const;

    ///@}

private:
    ///@name Member Variables
    ///@{

    const JacobianType& mrJ;

    const ShapeFunctionsLocalGradientType& mrDN_De;

    MatrixType mJtJ;

    MatrixType mCofactorJtJ;

    double mDetJtJ;

    ///@}
    ///@name Private Operations
    ///@{

    void Initialize();

    ///@}

}; // class LineSensitivityUtility

///@}

} // namespace Kratos

#endif // KRATOS_LINE_SENSITIVITY_UTILITY_H_INCLUDED