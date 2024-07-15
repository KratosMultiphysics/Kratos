// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter,
//                   Vahid Galavi
//

// System includes

// External includes

// Project includes
#include "custom_elements/geo_linear_truss_element_base.hpp"
#include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

namespace Kratos
{
//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementLinearBase<TDim, TNumNodes>::GeoTrussElementLinearBase(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoTrussElementBase<TDim, TNumNodes>(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementLinearBase<TDim, TNumNodes>::GeoTrussElementLinearBase(IndexType NewId,
                                                                      GeometryType::Pointer pGeometry,
                                                                      PropertiesType::Pointer pProperties)
    : GeoTrussElementBase<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElementLinearBase<TDim, TNumNodes>::Create(IndexType NewId,
                                                                    NodesArrayType const& rThisNodes,
                                                                    PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_intrusive<GeoTrussElementLinearBase>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElementLinearBase<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                    GeometryType::Pointer pGeom,
                                                                    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoTrussElementLinearBase>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElementLinearBase<TDim, TNumNodes>::~GeoTrussElementLinearBase()
{
}

//----------------------------------------------------------------------------------------
template class GeoTrussElementLinearBase<2, 2>;
template class GeoTrussElementLinearBase<3, 2>;

} // namespace Kratos.
