//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_U_P_AXISYM_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_U_P_AXISYM_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic isotropic constitutive law in 2D Neohookean Model (Axisymmetric)
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 * As there are no further parameters the functionality is limited
 * to large displacements elasticity.
 */

class KRATOS_API(SOLID_MECHANICS_APPLICATION) HyperElasticUPAxisym2DLaw : public HyperElasticUP3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticUPAxisym2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticUPAxisym2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticUPAxisym2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    HyperElasticUPAxisym2DLaw (const HyperElasticUPAxisym2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticUPAxisym2DLaw& operator=(const HyperElasticUPAxisym2DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~HyperElasticUPAxisym2DLaw();

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
    void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
                                       Vector& rStrainVector );


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                 Vector& rStrainVector );




    /**
     * Calculates the isochoric constitutive matrix
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const Matrix & rIsoStressMatrix,
            Matrix& rConstitutiveMatrix);


    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticUP3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticUP3DLaw )
    }



}; // Class HyperElasticUPAxisym2DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_U_P_AXISYM_2D_LAW_H_INCLUDED  defined 
