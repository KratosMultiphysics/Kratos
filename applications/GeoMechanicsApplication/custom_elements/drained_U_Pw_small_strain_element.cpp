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

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer DrainedUPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                       NodesArrayType const& ThisNodes,
                                                                       PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DrainedUPwSmallStrainElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties, this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer DrainedUPwSmallStrainElement<TDim, TNumNodes>::Create(IndexType             NewId,
                                                                       GeometryType::Pointer pGeom,
                                                                       PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new DrainedUPwSmallStrainElement(NewId, pGeom, pProperties,
                                                             this->GetStressStatePolicy().Clone()));
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void DrainedUPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim, TNumNodes>::InitializeSolutionStep(rCurrentProcessInfo);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
int DrainedUPwSmallStrainElement<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Verify generic variables
    int ierr = UPwBaseElement<TDim, TNumNodes>::Check(rCurrentProcessInfo);
    if (ierr != 0) return ierr;

    const PropertiesType& Prop = this->GetProperties();
    const GeometryType&   Geom = this->GetGeometry();

    if (Geom.DomainSize() < 1.0e-15)
        KRATOS_ERROR << "DomainSize < 1.0e-15 for the element " << this->Id() << std::endl;

    // Verify specific properties

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has(CONSTITUTIVE_LAW))
        << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = this->GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    if (TDim == 2) {
        KRATOS_ERROR_IF(strain_size < 3 || strain_size > 4)
            << "Wrong constitutive law used. This is a 2D element! expected "
               "strain size is 3 or 4 (el id = ) "
            << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6)
            << "Wrong constitutive law used. This is a 3D element! expected "
               "strain size is 6 (el id = ) "
            << this->Id() << std::endl;
    }

    // Check constitutive law
    if (mConstitutiveLawVector.size() > 0) {
        return mConstitutiveLawVector[0]->Check(Prop, Geom, rCurrentProcessInfo);
    }

    return ierr;

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void DrainedUPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddLHS(MatrixType& rLeftHandSideMatrix,
                                                                       ElementVariables& rVariables)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessMatrix(rLeftHandSideMatrix, rVariables);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
void DrainedUPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddRHS(VectorType& rRightHandSideVector,
                                                                       ElementVariables& rVariables,
                                                                       unsigned int      GPoint)
{
    KRATOS_TRY;

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddStiffnessForce(rRightHandSideVector,
                                                                          rVariables, GPoint);

    UPwSmallStrainElement<TDim, TNumNodes>::CalculateAndAddMixBodyForce(rRightHandSideVector, rVariables);

    KRATOS_CATCH("");
}

//----------------------------------------------------------------------------------------------------

template class DrainedUPwSmallStrainElement<2, 3>;
template class DrainedUPwSmallStrainElement<2, 4>;

template class DrainedUPwSmallStrainElement<3, 4>;
template class DrainedUPwSmallStrainElement<3, 8>;

} // Namespace Kratos
