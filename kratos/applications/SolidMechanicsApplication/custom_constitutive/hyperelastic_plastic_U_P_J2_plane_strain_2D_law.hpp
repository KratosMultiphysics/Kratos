//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_PLASTIC_U_P_J2_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_PLASTIC_U_P_J2_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_plastic_U_P_plane_strain_2D_law.hpp"



namespace Kratos
{
/**
 * Defines a hyperelastic-plastic isotropic constitutive law J2 in plane strain 2D 
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements 
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) HyperElasticPlasticUPJ2PlaneStrain2DLaw : public HyperElasticPlasticUPPlaneStrain2DLaw
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
     * Counted pointer of HyperElasticPlasticUPJ2PlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticPlasticUPJ2PlaneStrain2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticPlasticUPJ2PlaneStrain2DLaw();


    HyperElasticPlasticUPJ2PlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    HyperElasticPlasticUPJ2PlaneStrain2DLaw (const HyperElasticPlasticUPJ2PlaneStrain2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticPlasticUPJ2PlaneStrain2DLaw& operator=(const HyperElasticPlasticUPJ2PlaneStrain2DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~HyperElasticPlasticUPJ2PlaneStrain2DLaw();

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
    //int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);



    /**
     * Input and output
     */
    /**
     * Turn back information as a string.
     */
    //virtual String Info() const;
    /**
     * Print information about this object.
     */
    //virtual void PrintInfo(std::ostream& rOStream) const;
    /**
     * Print object's data.
     */
    //virtual void PrintData(std::ostream& rOStream) const;

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

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticPlasticUPPlaneStrain2DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticPlasticUPPlaneStrain2DLaw )
    }



}; // Class HyperElasticPlasticUPJ2PlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_PLASTIC_U_P_J2_PLANE_STRAIN_2D_LAW_H_INCLUDED defined
