//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_SHALLOW_ELEMENT_DATA_H_INCLUDED
#define KRATOS_SHALLOW_ELEMENT_DATA_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "custom_friction_laws/friction_law.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @brief Auxiliary class to store the nodal and Gauss point data, as well as the linearization matrices.
 * @ingroup ShallowWaterApplication
 * @author Miguel Maso Sotomayor
 * @tparam TNumNodes 
 */
template<std::size_t TNumNodes>
struct ShallowElementData
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    typedef Geometry<NodeType> GeometryType;

    static constexpr std::size_t NumNodes = TNumNodes;

    static constexpr std::size_t TLocalSize = 3 * TNumNodes;

    typedef array_1d<double, TLocalSize> LocalVectorType;

    typedef BoundedMatrix<double, TLocalSize, TLocalSize> LocalMatrixType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    ///@}
    ///@name Public member Variables
    ///@{

    bool integrate_by_parts;
    double stab_factor;
    double relative_dry_height;
    double gravity;
    double length;

    double height;
    array_1d<double,3> velocity;

    BoundedMatrix<double,3,3> A1;
    BoundedMatrix<double,3,3> A2;
    array_1d<double,3> b1;
    array_1d<double,3> b2;

    FrictionLaw::Pointer p_bottom_friction;

    array_1d<double, TNumNodes> nodal_z;

    ///@}
    ///@name Public member Operations
    ///@{

    void InitializeData(
        const GeometryType& rGeometry,
        const Properties& rProperties,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateGeometryData(const GeometryType& rGeometry);

    void SetNodalData(const GeometryType& rGeometry, int Step = 0);

    void UpdateGaussPointData(const array_1d<double,TNumNodes>& rN);

    ///@}
    ///@name Static Operations
    ///@{

    static double ShapeFunctionProduct(
        const array_1d<double,TNumNodes>& rN,
        const std::size_t I,
        const std::size_t J);

    static inline array_1d<double,3> VectorInnerProduct(
        const array_1d<array_1d<double,3>,TNumNodes>& rV,
        const array_1d<double,TNumNodes>& rN);

    static const Variable<double>& UnknownComponent(int Index);

    ///@}
    ///@name Access
    ///@{

    std::size_t NumberOfGaussPoints();

    const double GetWeight(int GaussPointIndex);

    const array_1d<double,TNumNodes> GetShapeFunctions(int GaussPointIndex);

    const BoundedMatrix<double,TNumNodes,2> GetShapeFunctionDerivatives(int GaussPointIndex);

    LocalVectorType& GetUnknownVector() const;

    ///@}

protected:

    ///@name Protected member Variables
    ///@{

    Vector mWeights;
    Matrix mN_Container;
    ShapeFunctionsGradientsType mDN_DX_Container;

    ///@}
};

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_SHALLOW_ELEMENT_DATA_H_INCLUDED  defined
