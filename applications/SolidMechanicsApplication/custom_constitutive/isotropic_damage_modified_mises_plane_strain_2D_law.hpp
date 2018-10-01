//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_ISOTROPIC_DAMAGE_MODIFIED_MISES_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_DAMAGE_MODIFIED_MISES_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plastic_plane_strain_2D_law.hpp"


namespace Kratos
{


class KRATOS_API(SOLID_MECHANICS_APPLICATION) IsotropicDamageModifiedMisesPlaneStrain2DLaw : public LinearElasticPlasticPlaneStrain2DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    typedef FlowRule::Pointer                FlowRulePointer;
    typedef YieldCriterion::Pointer    YieldCriterionPointer;
    typedef HardeningLaw::Pointer        HardeningLawPointer;
    typedef Properties::Pointer            PropertiesPointer;

    /**
     * Counted pointer of IsotropicDamageModifiedMisesPlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( IsotropicDamageModifiedMisesPlaneStrain2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    IsotropicDamageModifiedMisesPlaneStrain2DLaw();


    IsotropicDamageModifiedMisesPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    IsotropicDamageModifiedMisesPlaneStrain2DLaw (const IsotropicDamageModifiedMisesPlaneStrain2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //IsotropicDamageModifiedMisesPlaneStrain2DLaw& operator=(const IsotropicDamageModifiedMisesPlaneStrain2DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~IsotropicDamageModifiedMisesPlaneStrain2DLaw() override;

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;



    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //String Info() const override;
    /**
     * Print information about this object.
     */
    //void PrintInfo(std::ostream& rOStream) const override;
    /**
     * Print object's data.
     */
    //void PrintData(std::ostream& rOStream) const override;

protected:

    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{
    ///@}


    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElasticPlasticPlaneStrain2DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElasticPlasticPlaneStrain2DLaw )
    }



}; // Class IsotropicDamageModifiedMisesPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_DAMAGE_MODIFIED_MISES_PLANE_STRAIN_2D_LAW_H_INCLUDED defined
