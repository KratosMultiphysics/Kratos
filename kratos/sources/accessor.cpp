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
    double const& Accessor::GetProperty(
        const Variable<double> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        double aux;
        return aux;
    }

    /**
     * @brief This method implements the way to retrieve a Vector type variable
     */
    Vector const& Accessor::GetProperty(
        const Variable<Vector> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Vector aux(0);
        return aux;
    }

    /**
     * @brief This method implements the way to retrieve a bool type variable
     */
    bool const& Accessor::GetProperty(
        const Variable<bool> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        bool aux;
        return aux;
    }

    /**
     * @brief This method implements the way to retrieve a int type variable
     */
    int const& Accessor::GetProperty(
        const Variable<int> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        int aux;
        return aux;
    }

    /**
     * @brief This method implements the way to retrieve a Matrix type variable
     */
    Matrix const& Accessor::GetProperty(
        const Variable<Matrix> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Matrix mat(0, 0);
        return mat;
    }

    /**
     * @brief This method implements the way to retrieve a array_1d<double, 3 > type variable
     */
    array_1d<double, 3 > const& Accessor::GetProperty(
        const Variable<array_1d<double, 3 >> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 3 > mat;
        return mat;
    }

    /**
     * @brief This method implements the way to retrieve a array_1d<double, 6 > type variable
     */
    array_1d<double, 6 > const& Accessor::GetProperty(
        const Variable<array_1d<double, 6 >> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 6 > mat;
        return mat;
    }

    /**
     * @brief This method implements the way to retrieve a std::string type variable
     */
    std::string const& Accessor::GetProperty(
        const Variable<std::string> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        ) const
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        std::string str;
        return str;
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
    Accessor::Accessor(const Accessor& rOther)
    {
    }


} // namespace Kratos
