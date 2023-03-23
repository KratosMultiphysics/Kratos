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
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

    double Accessor::GetProperty(
        const Variable<double> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling to the virtual Accessor class..." << std::endl;
        return 0.0;
    }


} // namespace Kratos
