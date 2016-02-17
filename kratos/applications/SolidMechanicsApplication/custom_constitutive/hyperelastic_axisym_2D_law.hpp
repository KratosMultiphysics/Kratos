//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_AXISYM_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_AXISYM_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic isotropic constitutive law in 2D Neohookean Model (Axisymmetric)
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 * As there are no further parameters the functionality is limited
 * to large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) HyperElasticAxisym2DLaw : public HyperElastic3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticAxisym2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticAxisym2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticAxisym2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    HyperElasticAxisym2DLaw (const HyperElasticAxisym2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticAxisym2DLaw& operator=(const HyperElasticAxisym2DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~HyperElasticAxisym2DLaw();

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
        return 4;
    };

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);


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


    /**
     * Calculates the constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    void CalculateConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
                                      Matrix& rConstitutiveMatrix);


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw);
    }



}; // Class HyperElasticAxisym2DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_AXISYM_2D_LAW_H_INCLUDED  defined 
