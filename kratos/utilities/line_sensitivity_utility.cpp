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

#include "line_sensitivity_utility.h"
#include "utilities/math_utils.h"

namespace Kratos
{

LineSensitivityUtility::LineSensitivityUtility(
    const JacobianType& rJ, const ShapeFunctionsLocalGradientType& rDN_De)
    : mrJ(rJ), mrDN_De(rDN_De)
{
    this->Initialize();
}

void LineSensitivityUtility::CalculateSensitivity(
    const ShapeParameter Deriv, double& rIntegrationWeightSensitivity) const
{
    const unsigned int physical_dimension = mrJ.size1();
    const unsigned int parameter_dimension = mrJ.size2();
    MatrixType derivatives_of_jacobian(physical_dimension, parameter_dimension);
    noalias(derivatives_of_jacobian) = ZeroMatrix(physical_dimension, parameter_dimension);

    for (unsigned int i = 0; i < parameter_dimension; i++)
    {
        derivatives_of_jacobian(Deriv.Direction,i) = mrDN_De(Deriv.NodeIndex, i);
    }

    MatrixType jt_j_partial = prod(trans(derivatives_of_jacobian), mrJ);
    jt_j_partial += prod(trans(mrJ), derivatives_of_jacobian);

    double det_jt_j_partial = 0.0;
    for (unsigned int i = 0; i < parameter_dimension; i++)
    {
        for (unsigned int j = 0; j < parameter_dimension; j++)
        {
            det_jt_j_partial += mCofactorJtJ(i,j) * jt_j_partial(i,j);
        }
    }

    rIntegrationWeightSensitivity = det_jt_j_partial / (2.0 * std::sqrt(mDetJtJ));
}

void LineSensitivityUtility::Initialize()
{
    const unsigned int parameter_dimension = mrJ.size2();

    mJtJ.resize(parameter_dimension, parameter_dimension,false);
    noalias(mJtJ) = prod(trans(mrJ), mrJ);

    mCofactorJtJ.resize(parameter_dimension, parameter_dimension, false);
    noalias(mCofactorJtJ) = MathUtils<double>::CofactorMatrix(mJtJ);

    mDetJtJ = MathUtils<double>::Det(mJtJ);
}

}