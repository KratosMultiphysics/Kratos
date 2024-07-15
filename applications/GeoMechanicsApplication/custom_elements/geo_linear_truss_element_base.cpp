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
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CreateElementStiffnessMatrix(MatrixType& rLocalStiffnessMatrix,
                                                                              const ProcessInfo& rCurrentProcessInfo)

{
    KRATOS_TRY

    this->CalculateElasticStiffnessMatrix(rLocalStiffnessMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                                       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);
    // creating LHS
    this->CreateElementStiffnessMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
double GeoTrussElementLinearBase<TDim, TNumNodes>::CalculateLinearStrain()
{
    Vector current_disp = ZeroVector(TDim * TNumNodes);
    this->GetValuesVector(current_disp);
    FullDofMatrixType transformation_matrix;
    this->CreateTransformationMatrix(transformation_matrix);

    current_disp = prod(Matrix(trans(transformation_matrix)), current_disp);

    double length_0;
    if constexpr (TDim == 2) {
        length_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength2D2N(*this);
    } else if constexpr (TDim == 3) {
        length_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    } else {
        KRATOS_ERROR << "Dimension of truss element should be either 2D or 3D" << std::endl;
    }

    const double e = (current_disp[TDim] - current_disp[0]) / length_0;

    return e;
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElementLinearBase<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0]     = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template class GeoTrussElementLinearBase<2, 2>;
template class GeoTrussElementLinearBase<3, 2>;

} // namespace Kratos.
