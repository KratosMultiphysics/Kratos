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

// Application includes
#include "custom_elements/drained_U_Pw_small_strain_element.hpp"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer DrainedUPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                       NodesArrayType const& ThisNodes,
                                                                       PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DrainedUPwSmallStrainElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer DrainedUPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                       GeometryType::Pointer pGeom,
                                                                       PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DrainedUPwSmallStrainElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int DrainedUPwSmallStrainElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify generic variables
    int ierr = UPwBaseElement::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geometry   = this->GetGeometry();

    CheckUtilities::CheckDomainSize(r_geometry.DomainSize(), this->Id());

    const CheckProperties check_properties(r_properties, "property", CheckProperties::Bounds::AllExclusive);
    check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW);
    ierr = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, r_geometry, rCurrentProcessInfo);
    const auto expected_sizes = (TDim == 2 ? std::vector<std::size_t>{4} : std::vector<std::size_t>{6});
    ConstitutiveLawUtilities::CheckStrainSize(r_properties, expected_sizes, TDim, this->Id());

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void DrainedUPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                       ElementVariables& rVariables)
{
    KRATOS_TRY

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void DrainedUPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables,
                                                                       unsigned int      GPoint)
{
    KRATOS_TRY

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(rRightHandSideVector,
                                                                          rVariables, GPoint);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template class DrainedUPwSmallStrainElement<2, 3>;
template class DrainedUPwSmallStrainElement<2, 4>;

template class DrainedUPwSmallStrainElement<3, 4>;
template class DrainedUPwSmallStrainElement<3, 8>;

} // Namespace Kratos
