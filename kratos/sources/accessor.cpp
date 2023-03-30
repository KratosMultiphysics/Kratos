//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Riccardo Rossi
//                   Carlos Roig
//
//

// System includes

// External includes

// Project includes
#include "includes/accessor.h"

namespace Kratos
{

double Accessor::GetProperty(
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

Vector Accessor::GetProperty(
    const Variable<Vector>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    Vector aux(0);
    return aux;
}

/***********************************************************************************/
/***********************************************************************************/

bool Accessor::GetProperty(
    const Variable<bool>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

int Accessor::GetProperty(
    const Variable<int>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix Accessor::GetProperty(
    const Variable<Matrix>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    Matrix mat(0, 0);
    return mat;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> Accessor::GetProperty(
    const Variable<array_1d<double, 3>>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    array_1d<double, 3> mat;
    return mat;
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6> Accessor::GetProperty(
    const Variable<array_1d<double, 6>>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    array_1d<double, 6> mat;
    return mat;
}

/***********************************************************************************/
/***********************************************************************************/

std::string Accessor::GetProperty(
    const Variable<std::string>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
    return "string";
}

/***********************************************************************************/
/***********************************************************************************/

Accessor::Pointer Accessor::Clone() const
{
    return Kratos::make_shared<Accessor>(*this);
}
} // namespace Kratos
