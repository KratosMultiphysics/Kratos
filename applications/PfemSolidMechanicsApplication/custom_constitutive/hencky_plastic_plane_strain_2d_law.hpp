//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_HENCKY_PLASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_HENCKY_PLASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hencky_plastic_3d_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic-plastic isotropic constitutive law in plane strain 2D 
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements 
 */

class HenckyElasticPlasticPlaneStrain2DLaw : public HenckyElasticPlastic3DLaw
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
     * Counted pointer of HyperElasticPlasticPlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HenckyElasticPlasticPlaneStrain2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HenckyElasticPlasticPlaneStrain2DLaw();


    HenckyElasticPlasticPlaneStrain2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    HenckyElasticPlasticPlaneStrain2DLaw (const HenckyElasticPlasticPlaneStrain2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticPlasticPlaneStrain2DLaw& operator=(const HyperElasticPlasticPlaneStrain2DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~HenckyElasticPlasticPlaneStrain2DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension()
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()
    {
        return 3;
    };



    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
    //int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);



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

    /**
     * Calculates the GreenLagrange strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    virtual void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
            Vector& rStrainVector );


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    virtual void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector );

							      
    virtual void ConvertConstitutiveMatrixToAppropiateDimension(Matrix& rConstitutiveMatrix);


    virtual void CalculateHenckyMainStrain(const Matrix& rCauchyGreeMatrix, 
					   FlowRule::RadialReturnVariables& rReturnMappingVariables, 
					   Vector& rMainStrain);
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HenckyElasticPlastic3DLaw )
    }



}; // Class HyperElasticPlasticPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_PLASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED defined

