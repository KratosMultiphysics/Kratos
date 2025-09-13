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
#include "custom_elements/undrained_U_Pw_small_strain_element.hpp"
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UndrainedUPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         NodesArrayType const& ThisNodes,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UndrainedUPwSmallStrainElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties,
        this->GetStressStatePolicy().Clone(), this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UndrainedUPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                         GeometryType::Pointer pGeom,
                                                                         PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UndrainedUPwSmallStrainElement(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier()));
}

template <unsigned int TDim, unsigned int TNumNodes>
int UndrainedUPwSmallStrainElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const PropertiesType& r_properties = this->GetProperties();
    const GeometryType&   r_geometry   = this->GetGeometry();

    // Base class checks for positive area and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const auto element_Id = this->Id();
    CheckUtilities::CheckDomainSize(r_geometry.DomainSize(), element_Id);

    // Verify generic variables
    ierr = UPwBaseElement::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const CheckProperties check_properties(r_properties, "material properties at element",
                                           element_Id, CheckProperties::Bounds::AllExclusive);
    check_properties.Check(BULK_MODULUS_FLUID);

    check_properties.CheckAvailabilityAndSpecified(CONSTITUTIVE_LAW);
    ierr = r_properties[CONSTITUTIVE_LAW]->Check(r_properties, r_geometry, rCurrentProcessInfo);

    const auto expected_sizes = (TDim == 2 ? std::vector<std::size_t>{4} : std::vector<std::size_t>{6});
    ConstitutiveLawUtilities::CheckStrainSize(r_properties, expected_sizes, TDim, this->Id());

    return ierr;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UndrainedUPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                         ElementVariables& rVariables)
{
    KRATOS_TRY

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingMatrix(rLeftHandSideMatrix, rVariables);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCompressibilityMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void UndrainedUPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                         ElementVariables& rVariables,
                                                                         unsigned int GPoint)
{
    KRATOS_TRY

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(rRightHandSideVector,
                                                                          rVariables, GPoint);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCouplingTerms(rRightHandSideVector, rVariables);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddCompressibilityFlow(rRightHandSideVector, rVariables);

    KRATOS_CATCH("")
}

template class UndrainedUPwSmallStrainElement<2, 3>;
template class UndrainedUPwSmallStrainElement<2, 4>;
template class UndrainedUPwSmallStrainElement<3, 4>;
template class UndrainedUPwSmallStrainElement<3, 8>;

} // Namespace Kratos
