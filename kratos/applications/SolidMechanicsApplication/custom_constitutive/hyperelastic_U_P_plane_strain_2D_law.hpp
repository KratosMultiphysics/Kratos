//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_U_P_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_U_P_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic isotropic constitutive law in 2D Neohookean Model (Plane Strain)
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters:
 * 1) YOUNG MODULUS
 * 2) POISSON RATIO
 * As there are no further parameters the functionality is limited
 * to large displacements elasticity.
 */

class HyperElasticUPPlaneStrain2DLaw : public HyperElasticUP3DLaw
{
public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticUPPlaneStrain2DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticUPPlaneStrain2DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticUPPlaneStrain2DLaw();

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Copy constructor.
     */
    HyperElasticUPPlaneStrain2DLaw (const HyperElasticUPPlaneStrain2DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticUPPlaneStrain2DLaw& operator=(const HyperElasticUPPlaneStrain2DLaw& rOther);


    /**
     * Destructor.
     */
    virtual ~HyperElasticUPPlaneStrain2DLaw();

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
            const Vector & rIsoStressVector,
            Matrix& rConstitutiveMatrix);


    /**
     * Calculates the isochoric constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * @param rInverseDeformationGradientF
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const Vector & rIsoStressVector,
            const Matrix & rInverseDeformationGradientF,
            Matrix& rConstitutiveMatrix);



    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const GeometryType& rElementGeometry,
            const Vector & rShapeFunctions,
            Matrix& rConstitutiveMatrix);


    /**
     * Calculates the volumetric constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rInverseDeformationGradientF
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
            const GeometryType& rElementGeometry,
            const Vector & rShapeFunctions,
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



}; // Class HyperElasticUPPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_U_P_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined 
