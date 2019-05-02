//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

#if !defined(IGA_GEOMETRY_UTILITIES_H_INCLUDED)
#define IGA_GEOMETRY_UTILITIES_H_INCLUDED

#include "includes/model_part.h"

namespace Kratos
{

namespace IgaGeometryUtilities
{
    static void CalculateJacobian(
        const Element::GeometryType& rGeometry,
        const Matrix& rDN_De,
        const unsigned int rWorkingSpaceDimension,
        const unsigned int rLocalSpaceDimension,
        Matrix& rJacobian)
    {
        const unsigned int number_of_nodes = rGeometry.size();

        if ((rJacobian.size1() != rWorkingSpaceDimension) && (rJacobian.size2() != rLocalSpaceDimension))
            rJacobian.resize(rWorkingSpaceDimension, rLocalSpaceDimension);
        noalias(rJacobian) = ZeroMatrix(rWorkingSpaceDimension, rLocalSpaceDimension);

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            for (unsigned int k = 0; k<rWorkingSpaceDimension; k++)
            {
                for (unsigned int m = 0; m<rLocalSpaceDimension; m++)
                {
                    rJacobian(k, m) += (rGeometry[i]).Coordinates()[k] * rDN_De(i, m);
                }
            }
        }
    }

    static void CalculateInitialJacobian(
        const Element::GeometryType& rGeometry,
        const Matrix& rDN_De,
        const unsigned int rWorkingSpaceDimension,
        const unsigned int rLocalSpaceDimension,
        Matrix& rJacobian)
    {
        const unsigned int number_of_nodes = rGeometry.size();

        if ((rJacobian.size1() != rWorkingSpaceDimension) && (rJacobian.size2() != rLocalSpaceDimension))
            rJacobian.resize(rWorkingSpaceDimension, rLocalSpaceDimension);
        noalias(rJacobian) = ZeroMatrix(rWorkingSpaceDimension, rLocalSpaceDimension);

        for (unsigned int i = 0; i < number_of_nodes; i++)
        {
            for (unsigned int k = 0; k<rWorkingSpaceDimension; k++)
            {
                for (unsigned int m = 0; m<rLocalSpaceDimension; m++)
                {
                    rJacobian(k, m) += (rGeometry[i]).GetInitialPosition()[k] * rDN_De(i, m);
                }
            }
        }
    }

    static void CalculateHessian(
        const Element::GeometryType& rGeometry,
        const Matrix& rDDN_DDe,
        const unsigned int rWorkingSpaceDimension,
        Matrix& rHessian)
    {
        const unsigned int number_of_nodes = rGeometry.size();

        if ((rHessian.size1() != rWorkingSpaceDimension) && (rHessian.size2() != rWorkingSpaceDimension))
            rHessian.resize(rWorkingSpaceDimension, rWorkingSpaceDimension);
        rHessian = ZeroMatrix(rWorkingSpaceDimension, rWorkingSpaceDimension);

        for (int k = 0; k<number_of_nodes; k++)
        {
            const array_1d<double, 3> coords = rGeometry[k].Coordinates();

            for (int i = 0; i < rWorkingSpaceDimension; ++i)
            {
                for (int j = 0; j < rWorkingSpaceDimension; ++j)
                {
                    rHessian(i, j) += rDDN_DDe(k, j)*coords[i];
                }
            }
        }
    }


    static void CalculateCoordinates(
        const Element::GeometryType& rGeometry,
        const Vector& rN,
        const unsigned int rWorkingSpaceDimension,
        array_1d<double, 3>& rCoordinates)
    {
        const int& number_of_control_points = rGeometry.size();
        rCoordinates = ZeroVector(rWorkingSpaceDimension);
        for (int i = 0; i < rN.size(); i++)
        {
            const Node<3> & node = rGeometry[i];
            const array_1d<double, 3>& coords = node.Coordinates();

            for (int dim = 0; dim < rWorkingSpaceDimension; ++dim)
            {
                rCoordinates[dim] += rN[i] * coords[dim];
            }
        }
    }

    static void CalculateSolutionStepValue(
        const Variable<array_1d<double, 3>>& rVariable,
        const Element::GeometryType& rGeometry,
        const Vector& rN,
        const unsigned int rWorkingSpaceDimension,
        array_1d<double, 3>& rOutput,
        const int rSolutionStepIndex = 0)
    {
        rOutput = ZeroVector(3);
        for (int i = 0; i < rN.size(); i++)
        {
            const Node<3> & node = rGeometry[i];
            const array_1d<double, 3>& solution_step_value = node.FastGetSolutionStepValue(rVariable, rSolutionStepIndex);

            for (int dim = 0; dim < rWorkingSpaceDimension; ++dim)
            {
                rOutput[dim] += rN[i] * solution_step_value[dim];
            }
        }
    }

    static void CalculateValue(
        const Variable<Vector>& rVariable,
        const Element::GeometryType& rGeometry,
        const Vector& rN,
        const unsigned int rWorkingSpaceDimension,
        Vector& rOutput)
    {
        if (rOutput.size() != rWorkingSpaceDimension)
            rOutput.resize(rWorkingSpaceDimension);
        rOutput = ZeroVector(rWorkingSpaceDimension);

        for (int i = 0; i < rN.size(); i++)
        {
            const Node<3> & node = rGeometry[i];
            const array_1d<double, 3>& variable = node.GetValue(rVariable);

            for (int dim = 0; dim < rWorkingSpaceDimension; ++dim)
            {
                rOutput[dim] += rN[i] * variable[dim];
            }
        }
    }



    static double DistanceNodeToElement(
        const Element::GeometryType& rGeometry,
        const Vector& rN,
        const Node<3>::Pointer pNode,
        const unsigned int rWorkingSpaceDimension
    )
    {
        array_1d<double, 3> coordinates = ZeroVector(3);
        CalculateCoordinates(
            rGeometry,
            rN,
            rWorkingSpaceDimension,
            coordinates);

        double pythagoras = 0.0;
        for (int dim = 0; dim < rWorkingSpaceDimension; ++dim)
        {
            pythagoras += std::pow((coordinates[dim] - pNode->Coordinates()[dim]), 2);
        }

        return std::sqrt(pythagoras);
    }
}

} // namespace Kratos

#endif // IGA_GEOMETRY_UTILITIES_H_INCLUDED
