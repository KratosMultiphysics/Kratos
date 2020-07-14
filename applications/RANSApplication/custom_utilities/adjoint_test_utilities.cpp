//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//
//

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"

// Application includes
#include "rans_calculation_utilities.h"
#include "test_utilities.h"

// Include base h
#include "adjoint_test_utilities.h"

namespace Kratos
{
namespace RansApplicationTestUtilities
{
/**
 * @brief based on python's math.IsNear()
 */
bool IsNear(const double ValueA, const double ValueB, const double RelTol, const double AbsTol)
{
    if (std::isnan(ValueA) || std::isnan(ValueB))
    {
        return false;
    }
    // Special cases inf, -inf or exact:
    if (ValueA == ValueB)
    {
        return true;
    }
    // Regular floating point number:
    return std::abs(ValueA - ValueB) <=
           std::max(RelTol * std::max(std::abs(ValueA), std::abs(ValueB)), AbsTol);
}

void CheckNear(const double ValueA, const double ValueB, const double RelTol, const double AbsTol)
{
    if (IsNear(ValueA, ValueB, RelTol, AbsTol))
    {
        if (IsNear(ValueA, 0.0, 0.0, AbsTol) &&
            IsNear(ValueB, 0.0, 0.0, AbsTol) && ValueA != 0.0 && ValueB != 0.0)
        {
            // Warn if ValueA and ValueB are non-zero but below the absolute tolerance threshold.
            KRATOS_WARNING("CheckNear")
                << "Comparing values smaller than Tolerance. ValueA / ValueB < "
                   "Tolerance [ "
                << ValueA << " / " << ValueB << " < " << AbsTol << " ]\n";
        }
    }
    else
    {
        // Currently KRATOS_ERROR doesn't handle I/O formatting so stringstream
        // is used here to create the message.
        std::stringstream msg;
        msg << std::fixed << std::scientific << std::setprecision(20);
        msg << "Check failed because\n"
            << "\t ValueA = " << ValueA << " is not close to\n"
            << "\t ValueB = " << ValueB << std::setprecision(2) << '\n'
            << "\t with rel. tol. = " << RelTol << " and abs. tol. = " << AbsTol;
        KRATOS_ERROR << msg.str();
    }
}

void CheckNear(const Matrix& rA, const Matrix& rB, const double RelTol, const double AbsTol)
{
    KRATOS_ERROR_IF(rA.size1() != rB.size1())
        << "rA.size1() = " << rA.size1()
        << " is not equal to rB.size1() = " << rB.size1();

    KRATOS_ERROR_IF(rA.size2() != rB.size2())
        << "rA.size2() = " << rA.size2()
        << " is not equal to rB.size2() = " << rB.size2();

    for (std::size_t i = 0; i < rA.size1(); ++i)
        for (std::size_t j = 0; j < rA.size2(); ++j)
            CheckNear(rA(i, j), rB(i, j), RelTol, AbsTol);
}

void CheckNear(const Vector& rA, const Vector& rB, const double RelTol, const double AbsTol)
{
    KRATOS_ERROR_IF(rA.size() != rB.size())
        << "rA.size() = " << rA.size() << " is not equal to rB.size() = " << rB.size();

    for (std::size_t i = 0; i < rA.size(); ++i)
        CheckNear(rA(i), rB(i), RelTol, AbsTol);
}

std::function<double&(NodeType&)> GetPerturbationMethod(const Variable<double>& rPerturbationVariable)
{
    std::function<double&(NodeType&)> perturbation_method =
        [&rPerturbationVariable](NodeType& rNode) -> double& {
        return rNode.FastGetSolutionStepValue(rPerturbationVariable);
    };

    return perturbation_method;
}

std::function<double&(NodeType&, const int)> GetPerturbationMethod(
    const Variable<array_1d<double, 3>>& rPerturbationVariable)
{
    if (rPerturbationVariable == SHAPE_SENSITIVITY)
    {
        std::function<double&(NodeType&, const int)> perturbation_method =
            [](NodeType& rNode, const int iDim) -> double& {
            array_1d<double, 3>& r_coordinates = rNode.Coordinates();
            return r_coordinates[iDim];
        };
        return perturbation_method;
    }
    else
    {
        std::function<double&(NodeType&, const int)> perturbation_method =
            [&rPerturbationVariable](NodeType& rNode, const int iDim) -> double& {
            array_1d<double, 3>& r_vector =
                rNode.FastGetSolutionStepValue(rPerturbationVariable);
            return r_vector[iDim];
        };
        return perturbation_method;
    }
}

} // namespace RansApplicationTestUtilities
} // namespace Kratos
