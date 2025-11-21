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

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_Pw_diff_order_element.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

UpdatedLagrangianUPwDiffOrderElement::UpdatedLagrangianUPwDiffOrderElement(
    IndexType                                       NewId,
    GeometryType::Pointer                           pGeometry,
    std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
    std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : SmallStrainUPwDiffOrderElement(NewId, pGeometry, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

UpdatedLagrangianUPwDiffOrderElement::UpdatedLagrangianUPwDiffOrderElement(
    IndexType                                       NewId,
    GeometryType::Pointer                           pGeometry,
    PropertiesType::Pointer                         pProperties,
    std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
    std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
    : SmallStrainUPwDiffOrderElement(
          NewId, pGeometry, pProperties, std::move(pStressStatePolicy), std::move(pCoefficientModifier))
{
}

Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              const NodesArrayType& rNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    return Create(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

Element::Pointer UpdatedLagrangianUPwDiffOrderElement::Create(IndexType             NewId,
                                                              GeometryType::Pointer pGeom,
                                                              PropertiesType::Pointer pProperties) const
{
    return make_intrusive<UpdatedLagrangianUPwDiffOrderElement>(
        NewId, pGeom, pProperties, this->GetStressStatePolicy().Clone(),
        this->CloneIntegrationCoefficientModifier());
}


void UpdatedLagrangianUPwDiffOrderElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)
}

void UpdatedLagrangianUPwDiffOrderElement::load(Serializer& rSerializer){
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallStrainUPwDiffOrderElement)}

std::string UpdatedLagrangianUPwDiffOrderElement::Info() const
{
    const std::string constitutive_info =
        !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
    std::ostringstream oss;
    oss << "Updated Lagrangian U-Pw different order Element #" << this->Id()
        << "\nConstitutive law: " << constitutive_info;

    return oss.str();
}

void UpdatedLagrangianUPwDiffOrderElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void UpdatedLagrangianUPwDiffOrderElement::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}
} // Namespace Kratos
