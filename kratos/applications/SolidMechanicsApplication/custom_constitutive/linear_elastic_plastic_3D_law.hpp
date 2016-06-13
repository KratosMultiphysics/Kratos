//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_LINEAR_ELASTIC_PLASTIC_3D_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_PLASTIC_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"



namespace Kratos
{


class KRATOS_API(SOLID_MECHANICS_APPLICATION) LinearElasticPlastic3DLaw : public HyperElasticPlastic3DLaw
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
     * Counted pointer of LinearElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION( LinearElasticPlastic3DLaw );

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    LinearElasticPlastic3DLaw();


    LinearElasticPlastic3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

    /**
     * Copy constructor.
     */
    LinearElasticPlastic3DLaw (const LinearElasticPlastic3DLaw& rOther);


    /**
     * Assignment operator.
     */

    //LinearElasticPlastic3DLaw& operator=(const LinearElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const;

    /**
     * Destructor.
     */
    virtual ~LinearElasticPlastic3DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures);


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
     * Calculates the characteristic size of the element
     * @param rCharacteristicSize
     * @param DomainGeometry geometric information of the element
     */

    virtual void CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry );
    
    
    /**
     * Calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * Poisson ratio
     * @param E the Young's modulus
     * @param NU the Poisson ratio
     * @return the linear elastic constitutive matrix
     */

    virtual void CalculateLinearElasticMatrix( Matrix& rLinearElasticMatrix,const double& YoungModulus,const double& PoissonCoefficient );


    /**
     * Calculates the internal state variables and the stress vector
     * @param rReturnMappingVariables, plastic variables
     * @param rStressMatrix
     * @param rStressVector (same but in Voigt notation)
     * @param LinearElasticMatrix
     * @param StrainVector
     */
    
    virtual void CalculateReturnMapping( FlowRule::RadialReturnVariables& rReturnMappingVariables, Matrix& rStressMatrix,
                                            Vector& rStressVector, const Matrix& LinearElasticMatrix, const Vector& StrainVector );


    /**
     * Calculates the constitutive tensor: the secant or the tangent
     * @param rConstitutiveMatrix
     * @param rReturnMappingVariables, plastic variables
     * @param LinearElasticMatrix
     */

    virtual void CalculateConstitutiveTensor( Matrix& rConstitutiveMatrix, FlowRule::RadialReturnVariables& rReturnMappingVariables, const Matrix& LinearElasticMatrix );


    /**
     * Updates the internal state variables (to finalize the step)
     * @param rReturnMappingVariables, plastic variables
     * @param rStressVector
     * @param LinearElasticMatrix
     * @param StrainVector     
     */

    virtual void UpdateInternalStateVariables( FlowRule::RadialReturnVariables& rReturnMappingVariables, Vector& rStressVector,
                                            const Matrix& LinearElasticMatrix, const Vector& StrainVector );
                                  
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, HyperElasticPlastic3DLaw )
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, HyperElasticPlastic3DLaw )
    }



}; // Class LinearElasticPlastic3DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_PLASTIC_3D_LAW_H_INCLUDED defined
