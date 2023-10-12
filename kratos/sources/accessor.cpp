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

double Accessor::GetValue(
    const Variable<double>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

Vector Accessor::GetValue(
    const Variable<Vector>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return Vector();
}

/***********************************************************************************/
/***********************************************************************************/

bool Accessor::GetValue(
    const Variable<bool>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return false;
}

/***********************************************************************************/
/***********************************************************************************/

int Accessor::GetValue(
    const Variable<int>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix Accessor::GetValue(
    const Variable<Matrix>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return Matrix();
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> Accessor::GetValue(
    const Variable<array_1d<double, 3>>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return array_1d<double, 3>();
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6> Accessor::GetValue(
    const Variable<array_1d<double, 6>>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return array_1d<double, 6>();
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 4> Accessor::GetValue(
    const Variable<array_1d<double, 4>>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return array_1d<double, 4>();
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 9> Accessor::GetValue(
    const Variable<array_1d<double, 9>>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return array_1d<double, 9>();
}

/***********************************************************************************/
/***********************************************************************************/

std::string Accessor::GetValue(
    const Variable<std::string>& rVariable,
    const Properties& rProperties,
    const GeometryType& rGeometry,
    const Vector& rShapeFunctionVector,
    const ProcessInfo& rProcessInfo
    ) const
{
    KRATOS_ERROR << "You are calling the GetValue of the virtual Accessor class..." << std::endl;
    return "string";
}

/***********************************************************************************/
/***********************************************************************************/

Accessor::UniquePointer Accessor::Clone() const
{
    return Kratos::make_unique<Accessor>(*this);
}

} // namespace Kratos
