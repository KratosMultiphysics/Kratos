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

    /**
     * @brief This method implements the way to retrieve a double type variable
     */
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

    /**
     * @brief This method returns a pointer of the class
     */
    Accessor::Pointer Accessor::Clone() const
    {
        return Kratos::make_shared<Accessor>(*this);
    }

    /**
     * @brief Copy method
     */
    // Accessor::Accessor(const Accessor& rOther)
    // {
    // }


} // namespace Kratos
