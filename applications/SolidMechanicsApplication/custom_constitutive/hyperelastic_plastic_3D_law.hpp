//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"


namespace Kratos
{
/**
 * Defines a hyperelastic-plastic isotropic constitutive law in 3D 
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements 
 */

class HyperElasticPlastic3DLaw : public HyperElastic3DLaw
{
protected:

    /**
     * Parameters to be used in the volumetric and deviatoric split
     */
    struct VectorSplit
    {
        Vector  Isochoric;
        Vector  Volumetric;
    };

    struct MatrixSplit
    {
        Matrix  Isochoric;
        Matrix  Volumetric;
	Matrix  Plastic
    };


public:
    /**
     * Type Definitions
     */
    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(HyperElasticPlastic3DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    HyperElasticPlastic3DLaw();


    /**
     * Copy constructor.
     */
    HyperElasticPlastic3DLaw (const HyperElasticPlastic3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //HyperElasticPlastic3DLaw& operator=(const HyperElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~HyperElasticPlastic3DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues);

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

    YieldCriterion mpYieldCriterion;
	
    HardeningLaw   mpHardeningLaw;
	
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


    /**
     * Calculates the Pressure of the domain (element)
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rPressure the calculated pressure to be returned
     */
    double& CalculateDomainPressure (const GeometryType& rDomainGeometry,
                                     const Vector & rShapeFunctions,
                                     double & rPressure);


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
     * Constitutive isochoric component
     */

    double& IsochoricConstitutiveComponent( double & rCabcd,
                                            const MaterialResponseVariables& rElasticVariables,
                                            const Matrix & rIsoStressMatrix,
                                            const unsigned int& a, const unsigned int& b,
                                            const unsigned int& c, const unsigned int& d);

    /**
     * Constitutive isochoric component pull-back
     */

    double& IsochoricConstitutiveComponent( double & rCabcd,
                                            const MaterialResponseVariables& rElasticVariables,
                                            const Matrix & rIsoStressMatrix,
                                            const Matrix & rInverseDeformationGradientF,
                                            const unsigned int& a, const unsigned int& b,
                                            const unsigned int& c, const unsigned int& d);



    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rResult Matrix the result (Constitutive Matrix) will be stored in
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const GeometryType& rDomainGeometry,
            const Vector & rShapeFunctions,
            Matrix& rConstitutiveMatrix);


    /**
     * Calculates the volumetric constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rInverseDeformationGradientF
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
            const GeometryType& rDomainGeometry,
            const Vector & rShapeFunctions,
            Matrix& rConstitutiveMatrix);


    /**
     * Constitutive volumetric component
     */

    double& VolumetricConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const double & rPressure,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);

    /**
     * Constitutive volumetric component pull-back
     */

    double& VolumetricConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
            const double & rPressure,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);



    /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculateIsochoricStress( const MaterialResponseVariables & rElasticVariables,
					   FlowRule::RadialReturnVariables & rReturnMappingVariables,
                                           StressMeasure rStressMeasure,
                                           Vector& rIsoStressVector);

    /**
     * Calculates the volumetric stress vector
     * @param rElasticVariables
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rVolStressVector vector where the stress result is stored
     */
    virtual void CalculateVolumetricStress( const MaterialResponseVariables & rElasticVariables,
                                            const GeometryType& rDomainGeometry,
                                            const Vector & rShapeFunctions,
                                            Vector& rVolStressVector );




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



}; // Class HyperElasticPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED defined
