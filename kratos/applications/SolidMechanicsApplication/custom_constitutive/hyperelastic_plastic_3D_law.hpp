//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  1.0 $
//
//

#if !defined (KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/custom_flow_rules/flow_rule.hpp"
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

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticPlastic3DLaw );

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
    virtual void CalculateMaterialResponseKirchhoff (Parameters & rValues);


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
    
    FlowRulePointer       mpFlowRule;

    YieldCriterionPointer mpYieldCriterion;
	
    HardeningLawPointer   mpHardeningLaw;
	
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{


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
     * Constitutive volumetric component
     */

    double& PlasticConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rIsoStressMatrix,
            const FlowRule::PlasticFactors & rScalingFactors,			 
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);


    /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressMatrix matrix where the stress result is stored
     * @param rIsoStressVector vector where the stress result is stored
     */
    virtual void CalculatePlasticIsochoricStress( MaterialResponseVariables & rElasticVariables,
						  FlowRule::RadialReturnVariables & rReturnMappingVariables,
						  StressMeasure rStressMeasure,
						  Matrix& rIsoStressMatrix,
						  Vector& rIsoStressVector);



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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElastic3DLaw )
	
	rSerializer.save("mElasticLeftCauchyGreen",mElasticLeftCauchyGreen);
	rSerializer.save("mpFlowRule",mpFlowRule);
	rSerializer.save("mpYieldCriterion",mpYieldCriterion);
	rSerializer.save("mpHardeningLaw",mpHardeningLaw);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElastic3DLaw )
	  
	rSerializer.load("mElasticLeftCauchyGreen",mElasticLeftCauchyGreen);
	rSerializer.load("mpFlowRule",mpFlowRule);
	rSerializer.load("mpYieldCriterion",mpYieldCriterion);
	rSerializer.load("mpHardeningLaw",mpHardeningLaw);
    }



}; // Class HyperElasticPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_HYPERELASTIC_PLASTIC_3D_LAW_H_INCLUDED defined
