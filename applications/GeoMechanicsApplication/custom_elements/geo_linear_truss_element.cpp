// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes

// External includes

// Project includes
#include "custom_elements/geo_linear_truss_element.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoLinearTrussElement<TDim, TNumNodes>::GeoLinearTrussElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoTrussElementLinearBase<TDim, TNumNodes>(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoLinearTrussElement<TDim, TNumNodes>::GeoLinearTrussElement(IndexType               NewId,
                                                              GeometryType::Pointer   pGeometry,
                                                              PropertiesType::Pointer pProperties)
    : GeoTrussElementLinearBase<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoLinearTrussElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                NodesArrayType const& rThisNodes,
                                                                PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_intrusive<GeoLinearTrussElement>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoLinearTrussElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                GeometryType::Pointer pGeom,
                                                                PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoLinearTrussElement>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoLinearTrussElement<TDim, TNumNodes>::~GeoLinearTrussElement()
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoLinearTrussElement<TDim, TNumNodes>::ResetConstitutiveLaw()
{
    KRATOS_TRY

    mInternalStresses                  = ZeroVector(mStressVectorSize);
    mInternalStressesFinalized         = ZeroVector(mStressVectorSize);
    mInternalStressesFinalizedPrevious = ZeroVector(mStressVectorSize);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoLinearTrussElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeoTrussElementLinearBase<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

    if (rCurrentProcessInfo.Has(RESET_DISPLACEMENTS)) {
        bool ResetDisplacement = rCurrentProcessInfo[RESET_DISPLACEMENTS];
        if (ResetDisplacement) {
            mInternalStressesFinalizedPrevious = mInternalStressesFinalized;
        } else {
            mInternalStressesFinalized = mInternalStressesFinalizedPrevious;
        }
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoLinearTrussElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points =
        this->GetGeometry().IntegrationPoints();
    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == FORCE) {
        BoundedVector<double, 3> truss_forces = ZeroVector(3);
        const double             A            = this->GetProperties()[CROSS_AREA];

        double prestress = 0.00;
        if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
        Vector temp_strain = ZeroVector(mStressVectorSize);
        Vector temp_stress = ZeroVector(mStressVectorSize);
        temp_strain[0]     = this->CalculateLinearStrain();
        Values.SetStrainVector(temp_strain);
        Values.SetStressVector(temp_stress);
        mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

        temp_stress += mInternalStressesFinalizedPrevious;

        truss_forces[0] = (temp_stress[0] + prestress) * A;

        rOutput[0] = truss_forces;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoLinearTrussElement<TDim, TNumNodes>::UpdateInternalForces(FullDofVectorType& rInternalForces,
                                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Vector temp_strain = ZeroVector(mStressVectorSize);
    Vector temp_stress = ZeroVector(mStressVectorSize);
    temp_strain[0]     = this->CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    mInternalStresses = temp_stress;

    temp_stress += mInternalStressesFinalizedPrevious;

    Vector temp_internal_stresses = ZeroVector(TDim * TNumNodes);
    temp_internal_stresses[0]     = -1.0 * temp_stress[0];
    temp_internal_stresses[TDim]  = temp_stress[0];

    rInternalForces = temp_internal_stresses * this->GetProperties()[CROSS_AREA];

    FullDofMatrixType transformation_matrix;
    this->CreateTransformationMatrix(transformation_matrix);

    rInternalForces = prod(transformation_matrix, rInternalForces);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoLinearTrussElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    GeoTrussElementLinearBase<TDim, TNumNodes>::FinalizeSolutionStep(rCurrentProcessInfo);
    mInternalStressesFinalized = mInternalStresses + mInternalStressesFinalizedPrevious;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template class GeoLinearTrussElement<2, 2>;
template class GeoLinearTrussElement<3, 2>;

} // namespace Kratos.
