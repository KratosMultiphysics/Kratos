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
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0.0;
    }

    double Accessor::GetProperty(
        const Variable<double> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0.0;
    }

    double Accessor::GetProperty(
        const Variable<double> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0.0;
    }

    double Accessor::GetProperty(
        const Variable<double> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0.0;
    }

    /**
     * @brief This method implements the way to retrieve a Vector type variable
     */
    Vector Accessor::GetProperty(
        const Variable<Vector> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Vector aux(0);
        return aux;
    }

    Vector Accessor::GetProperty(
        const Variable<Vector> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Vector aux(0);
        return aux;
    }

    Vector Accessor::GetProperty(
        const Variable<Vector> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Vector aux(0);
        return aux;
    }

    Vector Accessor::GetProperty(
        const Variable<Vector> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Vector aux(0);
        return aux;
    }

    /**
     * @brief This method implements the way to retrieve a bool type variable
     */
    bool Accessor::GetProperty(
        const Variable<bool> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return false;
    }

    bool Accessor::GetProperty(
        const Variable<bool> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return false;
    }

    bool Accessor::GetProperty(
        const Variable<bool> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return false;
    }

    bool Accessor::GetProperty(
        const Variable<bool> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return false;
    }

    /**
     * @brief This method implements the way to retrieve a int type variable
     */
    int Accessor::GetProperty(
        const Variable<int> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0;
    }

    int Accessor::GetProperty(
        const Variable<int> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0;
    }

    int Accessor::GetProperty(
        const Variable<int> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0;
    }

    int Accessor::GetProperty(
        const Variable<int> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        return 0;
    }

    /**
     * @brief This method implements the way to retrieve a Matrix type variable
     */
    Matrix Accessor::GetProperty(
        const Variable<Matrix> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Matrix mat(0, 0);
        return mat;
    }

    Matrix Accessor::GetProperty(
        const Variable<Matrix> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Matrix mat(0, 0);
        return mat;
    }

    Matrix Accessor::GetProperty(
        const Variable<Matrix> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Matrix mat(0, 0);
        return mat;
    }

    Matrix Accessor::GetProperty(
        const Variable<Matrix> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        Matrix mat(0, 0);
        return mat;
    }

    /**
     * @brief This method implements the way to retrieve a array_1d<double, 3 > type variable
     */
    array_1d<double, 3 > Accessor::GetProperty(
        const Variable<array_1d<double, 3 >> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 3 > mat;
        return mat;
    }

    array_1d<double, 3 > Accessor::GetProperty(
        const Variable<array_1d<double, 3 >> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 3 > mat;
        return mat;
    }

    array_1d<double, 3 > Accessor::GetProperty(
        const Variable<array_1d<double, 3 >> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 3 > mat;
        return mat;
    }

    array_1d<double, 3 > Accessor::GetProperty(
        const Variable<array_1d<double, 3 >> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 3 > mat;
        return mat;
    }

    /**
     * @brief This method implements the way to retrieve a array_1d<double, 6 > type variable
     */
    array_1d<double, 6 > Accessor::GetProperty(
        const Variable<array_1d<double, 6 >> &rVariable,
        const Properties &rProperties,
        const GeometryType &rGeometry,
        const Vector &rShapeFunctionVector,
        const ProcessInfo &rProcessInfo
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 6 > mat;
        return mat;
    }

    array_1d<double, 6 > Accessor::GetProperty(
        const Variable<array_1d<double, 6 >> &rVariable,
        const Properties &rProperties
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 6 > mat;
        return mat;
    }

    array_1d<double, 6 > Accessor::GetProperty(
        const Variable<array_1d<double, 6 >> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 6 > mat;
        return mat;
    }

    array_1d<double, 6 > Accessor::GetProperty(
        const Variable<array_1d<double, 6 >> &rVariable,
        const Properties &rProperties,
        NodeType& rThisNode,
        IndexType SolutionStepIndex
        )
    {
        KRATOS_ERROR << "You are calling the GetProperty of the virtual Accessor class..." << std::endl;
        array_1d<double, 6 > mat;
        return mat;
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
