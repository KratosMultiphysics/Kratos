//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:              IPouplana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_ISOTROPIC_DAMAGE_SIMO_JU_3D_LAW_H_INCLUDED)
#define  KRATOS_ISOTROPIC_DAMAGE_SIMO_JU_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plastic_3D_law.hpp"


namespace Kratos
{


class KRATOS_API(SOLID_MECHANICS_APPLICATION) IsotropicDamageSimoJu3DLaw : public LinearElasticPlastic3DLaw
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
     * Counted pointer of IsotropicDamageSimoJu3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( IsotropicDamageSimoJu3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    IsotropicDamageSimoJu3DLaw();


    IsotropicDamageSimoJu3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    IsotropicDamageSimoJu3DLaw (const IsotropicDamageSimoJu3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //IsotropicDamageSimoJu3DLaw& operator=(const IsotropicDamageSimoJu3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    ~IsotropicDamageSimoJu3DLaw() override;

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

    /**
     * Calculates the characteristic size of the element.
     * It is the diameter of a sphere with the same volume as the element
     * @param rCharacteristicSize, the diameter of the sphere
     * @param DomainGeometry geometric information of the element
     */

    void CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry ) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LinearElasticPlastic3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LinearElasticPlastic3DLaw )
    }



}; // Class IsotropicDamageSimoJu3DLaw
}  // namespace Kratos.
#endif // KRATOS_ISOTROPIC_DAMAGE_SIMO_JU_3D_LAW_H_INCLUDED defined
