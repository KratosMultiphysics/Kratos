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
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"
#include "includes/constitutive_law.h"


namespace Kratos
{
/**
 * Defines a hyperelastic-plastic isotropic constitutive law in 3D 
 * With stress split in an isochoric and volumetric parts
 * This material law is defined by the parameters needed by the yield criterion:

 * The functionality is limited to large displacements 
 */

class HyperElasticPlastic3DLaw : public ConstitutiveLaw
{
protected:


    struct MaterialResponseVariables
    {

        //general material properties
        double LameMu;
        double LameLambda;

        //kinematic properties
        double J_pow13;
        double DeterminantF0;
        double traceCG;           //LeftCauchyGreen or RightCauchyGreen
        Matrix CauchyGreenMatrix; //LeftCauchyGreen or InverseRightCauchyGreen
        Matrix IdentityMatrix;

        //element properties
        const Vector*        mpShapeFunctionsValues;
        const GeometryType*  mpElementGeometry;

    public: 
      void SetShapeFunctionsValues (const Vector& rShapeFunctionsValues)      {mpShapeFunctionsValues=&rShapeFunctionsValues;};
      void SetElementGeometry      (const GeometryType& rElementGeometry)     {mpElementGeometry =&rElementGeometry;};
      const Vector& GetShapeFunctionsValues      () const {return *mpShapeFunctionsValues;};
      const GeometryType& GetElementGeometry     () const {return *mpElementGeometry;};
      
    };

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
        Matrix  Plastic;
    };


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


    HyperElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

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
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension()
    {
        return 3;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize()
    {
        return 6;
    };


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);


    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );


    void SetValue( const Variable<double>& rVariable,
                   const double& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial( const Properties& rMaterialProperties,
                             const GeometryType& rElementGeometry,
                             const Vector& rShapeFunctionsValues );


    void InitializeSolutionStep( const Properties& rMaterialProperties,
                                 const GeometryType& rElementGeometry, //this is just to give the array of nodes
                                 const Vector& rShapeFunctionsValues ,
                                 const ProcessInfo& rCurrentProcessInfo);

    void FinalizeSolutionStep( const Properties& rMaterialProperties,
                               const GeometryType& rElementGeometry, //this is just to give the array of nodes
                               const Vector& rShapeFunctionsValues ,
                               const ProcessInfo& rCurrentProcessInfo);

   /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK1 (Parameters & rValues);

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponsePK2 (Parameters & rValues);

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (Parameters & rValues); //Ll:

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseCauchy (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK2 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    virtual void FinalizeMaterialResponseKirchhoff (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy (Parameters & rValues);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);



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
  
    Matrix mElasticLeftCauchyGreen;
    
    FlowRulePointer mpFlowRule;

    YieldCriterionPointer mpYieldCriterion;
	
    HardeningLawPointer   mpHardeningLaw;
	
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
     */
    Matrix& DeformationGradient3D (Matrix & Matrix2D);

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
     * Calculates the isochoric constitutive matrix
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
            const Matrix & rIsoStressMatrix,
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
						       const Matrix & rInverseDeformationGradientF,
						       const Matrix & rIsoStressMatrix,
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
                                            const Matrix & rInverseDeformationGradientF,
					    const Matrix & rIsoStressMatrix,
					    const unsigned int& a, const unsigned int& b,
                                            const unsigned int& c, const unsigned int& d);



    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
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
							Matrix& rConstitutiveMatrix);


    /**
     * Constitutive volumetric component
     */

    double& VolumetricConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Vector& rFactors,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);

    /**
     * Constitutive volumetric component pull-back
     */

    double& VolumetricConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
            const Vector & rFactors,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);


    /**
     * Calculates the plastic constitutive matrix
     * @param rElasticVariables
     * @param rReturnMappingVariables, plastic variables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						     FlowRule::RadialReturnVariables & rReturnMappingVariables,		     
						     Matrix& rConstitutiveMatrix);


    /**
     * Calculates the plastic constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rReturnMappingVariables, plastic variables
     * @param rInverseDeformationGradientF
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
    virtual void CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						     const Matrix & rInverseDeformationGradientF,
						     FlowRule::RadialReturnVariables & rReturnMappingVariables,
						     Matrix& rConstitutiveMatrix);


    /**
     * Constitutive volumetric component
     */

    double& PlasticConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rIsoStressMatrix,
            const FlowRule::PlasticFactors & rScalingFactors,			 
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);

    /**
     * Constitutive volumetric component pull-back
     */

    double& PlasticConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
  	    const Matrix & rIsoStressMatrix,
	    const FlowRule::PlasticFactors & rScalingFactors,
	    const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);


    /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculateIsochoricStress( MaterialResponseVariables & rElasticVariables,
					   FlowRule::RadialReturnVariables & rReturnMappingVariables,
                                           StressMeasure rStressMeasure,
					   Matrix& rIsoStressMatrix,
                                           Vector& rIsoStressVector);

    /**
     * Calculates the volumetric stress vector
     * @param rElasticVariables
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rVolStressVector vector where the stress result is stored
     */
    virtual void CalculateVolumetricStress( const MaterialResponseVariables & rElasticVariables,
                                            Vector& rVolStressVector );


    /**
     * Calculates the Temperature of the domain (element)
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rTemperature the calculated temperature to be returned
     */
    virtual double& CalculateDomainTemperature (const MaterialResponseVariables & rElasticVariables,
                                     double & rTemperature);
    /**
     * Calculates the Pressure of the domain (element)
     * @param rElementGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rPressure the calculated pressure to be returned
     */
    virtual double& CalculateDomainPressure (const MaterialResponseVariables & rElasticVariables,
                                     double & rPressure);

    /**
     * Calculates the Volumetric part factors
     * @param rElasticResponseVariables the material variables
     * @param rFactors Volumetric stress factors
     */
    virtual Vector&  CalculateDomainPressureFactors (const MaterialResponseVariables & rElasticVariables,
					     Vector & rFactors);


    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
    virtual bool CheckParameters(Parameters& rValues);

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
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    }



}; // Class HyperElasticPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED defined
