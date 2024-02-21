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
#include "custom_elements/geo_truss_element.hpp"
#include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/define.h"

namespace Kratos
{
//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElement<TDim, TNumNodes>::GeoTrussElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : GeoTrussElementBase<TDim, TNumNodes>(NewId, pGeometry)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElement<TDim, TNumNodes>::GeoTrussElement(IndexType               NewId,
                                                  GeometryType::Pointer   pGeometry,
                                                  PropertiesType::Pointer pProperties)
    : GeoTrussElementBase<TDim, TNumNodes>(NewId, pGeometry, pProperties)
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElement<TDim, TNumNodes>::Create(IndexType               NewId,
                                                          NodesArrayType const&   rThisNodes,
                                                          PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = this->GetGeometry();
    return Kratos::make_intrusive<GeoTrussElement>(NewId, rGeom.Create(rThisNodes), pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoTrussElement<TDim, TNumNodes>::Create(IndexType               NewId,
                                                          GeometryType::Pointer   pGeom,
                                                          PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoTrussElement>(NewId, pGeom, pProperties);
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
GeoTrussElement<TDim, TNumNodes>::~GeoTrussElement()
{
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElement<TDim, TNumNodes>::ResetConstitutiveLaw()
{
    KRATOS_TRY

    mInternalStresses                  = ZeroVector(mStressVectorSize);
    mInternalStressesFinalized         = ZeroVector(mStressVectorSize);
    mInternalStressesFinalizedPrevious = ZeroVector(mStressVectorSize);

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    GeoTrussElementBase<TDim, TNumNodes>::Initialize(rCurrentProcessInfo);

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
void GeoTrussElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                    std::vector<Vector>& rOutput,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points =
        this->GetGeometry().IntegrationPoints();

    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        Vector strain = ZeroVector(TDim);
        strain[0]     = this->CalculateGreenLagrangeStrain();
        rOutput[0]    = strain;
    }
    if (rVariable == PK2_STRESS_VECTOR) {
        ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
        Vector temp_strain = ZeroVector(1);
        temp_strain[0]     = this->CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);

        array_1d<double, 3> temp_internal_stresses = ZeroVector(3);
        mpConstitutiveLaw->CalculateValue(Values, FORCE, temp_internal_stresses);

        for (unsigned int i = 0; i < TDim; ++i)
            temp_internal_stresses[i] += mInternalStressesFinalizedPrevious[i];

        rOutput[0] = temp_internal_stresses;
    }
    if (rVariable == CAUCHY_STRESS_VECTOR) {
        ProcessInfo temp_process_information;

        ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), temp_process_information);
        Vector temp_strain = ZeroVector(1);
        temp_strain[0]     = this->CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);
        array_1d<double, 3> temp_internal_stresses = ZeroVector(3);
        mpConstitutiveLaw->CalculateValue(Values, FORCE, temp_internal_stresses);

        for (unsigned int i = 0; i < TDim; ++i)
            temp_internal_stresses[i] += mInternalStressesFinalizedPrevious[i];

        const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
        const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

        rOutput[0] = temp_internal_stresses * l / L0;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                                                    std::vector<array_1d<double, 3>>& rOutput,
                                                                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points =
        this->GetGeometry().IntegrationPoints();
    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == FORCE) {
        BoundedVector<double, TDim> truss_forces = ZeroVector(TDim);
        const double                A            = this->GetProperties()[CROSS_AREA];

        double prestress = 0.00;
        if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);

        ProcessInfo temp_process_information;
        ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), temp_process_information);

        Vector temp_strain = ZeroVector(1);
        temp_strain[0]     = this->CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);
        array_1d<double, 3> temp_internal_stresses = ZeroVector(3);
        mpConstitutiveLaw->CalculateValue(Values, FORCE, temp_internal_stresses);

        for (unsigned int i = 0; i < TDim; ++i)
            temp_internal_stresses[i] += mInternalStressesFinalizedPrevious[i];

        truss_forces[0] = ((temp_internal_stresses[0] + prestress) * l * A) / L0;

        rOutput[0] = truss_forces;
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElement<TDim, TNumNodes>::UpdateInternalForces(BoundedVector<double, TDim * TNumNodes>& rInternalForces,
                                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    BoundedMatrix<double, TDim * TNumNodes, TDim * TNumNodes> transformation_matrix =
        ZeroMatrix(TDim * TNumNodes, TDim * TNumNodes);

    this->CreateTransformationMatrix(transformation_matrix);

    const double l  = StructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double L0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double A  = this->GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (this->GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = this->GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    ConstitutiveLaw::Parameters Values(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0]     = this->CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);

    mInternalStresses = temp_stress;

    temp_stress += mInternalStressesFinalizedPrevious;

    const double normal_force = ((temp_stress[0] + prestress) * l * A) / L0;

    // internal force vectors
    BoundedVector<double, TDim * TNumNodes> f_local = ZeroVector(TDim * TNumNodes);
    f_local[0]                                      = -1.00 * normal_force;
    f_local[TDim]                                   = 1.00 * normal_force;
    rInternalForces                                 = ZeroVector(TDim * TNumNodes);
    noalias(rInternalForces)                        = prod(transformation_matrix, f_local);
    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void GeoTrussElement<TDim, TNumNodes>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    GeoTrussElementBase<TDim, TNumNodes>::FinalizeSolutionStep(rCurrentProcessInfo);
    mInternalStressesFinalized = mInternalStresses + mInternalStressesFinalizedPrevious;

    KRATOS_CATCH("");
}

//--------------------------------------------------------------------------------------------
template class GeoTrussElement<2, 2>;
template class GeoTrussElement<3, 2>;

} // namespace Kratos.
