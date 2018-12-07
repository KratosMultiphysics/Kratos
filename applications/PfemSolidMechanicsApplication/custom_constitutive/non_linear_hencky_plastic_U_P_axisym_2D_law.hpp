//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_NON_LINEAR_HENCKY_PLASTIC_AXISYM_U_P_2D_LAW_H_INCLUDED)
#define       KRATOS_NON_LINEAR_HENCKY_PLASTIC_AXISYM_U_P_2D_LAW_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/non_linear_hencky_plastic_U_P_3D_law.hpp"

namespace Kratos
{
/**
  *

 */


class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) NonLinearHenckyElasticPlasticUPAxisym2DLaw: public NonLinearHenckyElasticPlasticUP3DLaw
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
     * Counted pointer of HyperElasticPlastic3DLaw
     */

    KRATOS_CLASS_POINTER_DEFINITION(NonLinearHenckyElasticPlasticUPAxisym2DLaw);

    /**
     * Life Cycle
     */

    /**
     * Default constructor.
     */
    NonLinearHenckyElasticPlasticUPAxisym2DLaw();


    NonLinearHenckyElasticPlasticUPAxisym2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw);

    /**
     * Copy constructor.
     */
    NonLinearHenckyElasticPlasticUPAxisym2DLaw (const NonLinearHenckyElasticPlasticUPAxisym2DLaw& rOther)  ;


    /**
     * Assignment operator.
     */

    //HyperElasticPlastic3DLaw& operator=(const HyperElasticPlastic3DLaw& rOther);

    /**
     * Clone function (has to be implemented by any derived class)
     * @return a pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Destructor.
     */
    virtual ~NonLinearHenckyElasticPlasticUPAxisym2DLaw();

    /**
     * Operators
     */

    /**
     * Operations needed by the base class:
     */

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 4;
    };


    void GetLawFeatures(Features& rFeatures) override;

/*    bool Has( const Variable<double>& rThisVariable );
    bool Has( const Variable<Vector>& rThisVariable );
    bool Has( const Variable<Matrix>& rThisVariable );

    double& GetValue( const Variable<double>& rThisVariable, double& rValue );
    Vector& GetValue( const Variable<Vector>& rThisVariable, Vector& rValue );
    Matrix& GetValue( const Variable<Matrix>& rThisVariable, Matrix& rValue );


    void SetValue( const Variable<double>& rVariable,
                   const double& Value,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Vector>& rThisVariable,
                   const Vector& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
    void SetValue( const Variable<Matrix>& rThisVariable,
                   const Matrix& rValue,
                   const ProcessInfo& rCurrentProcessInfo );
*/    /**
     * Material parameters are inizialized
     */
//    void InitializeMaterial( const Properties& rProps,
//                             const GeometryType& rGeom,
//                             const Vector& rShapeFunctionsValues );


/*    void InitializeSolutionStep( const Properties& props,
                                 const GeometryType& geom, //this is just to give the array of nodes
                                 const Vector& ShapeFunctionsValues ,
                                 const ProcessInfo& CurrentProcessInfo);

    void FinalizeSolutionStep( const Properties& props,
                               const GeometryType& geom, //this is just to give the array of nodes
                               const Vector& ShapeFunctionsValues ,
                               const ProcessInfo& CurrentProcessInfo);
*/
   /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
//    void CalculateMaterialResponsePK1 (Parameters & rValues);
// Esta función sirve como está definida en la clase madre

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
//    void CalculateMaterialResponsePK2 (Parameters & rValues);
// Esta función sirve como está definida en la clase madre

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
//    void CalculateMaterialResponseKirchhoff (Parameters & rValues);


    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
//    void CalculateMaterialResponseCauchy (Parameters & rValues);
// Esta función sirve como está definida

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
//    void FinalizeMaterialResponsePK1 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
//    void FinalizeMaterialResponsePK2 (Parameters & rValues);

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues
      * @see   Parameters
      */
//    void FinalizeMaterialResponseCauchy (Parameters & rValues);

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param props
     * @param geom
     * @param CurrentProcessInfo
     * @return
     */
//    int Check(const Properties& rProperties, const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo);



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

//    Matrix mElasticLeftCauchyGreen;

//    FlowRulePointer mpFlowRule;

//    YieldCriterionPointer mpYieldCriterion;

//    HardeningLawPointer   mpHardeningLaw;

    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    Matrix SetConstitutiveMatrixToAppropiateDimension(const Matrix& rElastoPlasticTangentMatrix) override;


    /**
     * Takes a matrix 2x2 and transforms it to a 3x3 adding a 3rd row and a 3rd column with a 1 in the diagonal
     */
//    Matrix& DeformationGradient3D (Matrix & Matrix2D);

    /**
     * Calculates the GreenLagrange strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateGreenLagrangeStrain( const Matrix & rRightCauchyGreen,
            Vector& rStrainVector ) override;


    /**
     * Calculates the Almansi strains
     * @param rRightCauchyGreen
     * @param rStrainVector
     */
    void CalculateAlmansiStrain( const Matrix & rLeftCauchyGreen,
                                         Vector& rStrainVector ) override;
    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS
        in a pullback fashion
    */



    /** First and secod term of the CONSISTENT ELASTOPLASTIC MATRIX FOR LARGE DEFORMATIONS
        in the actual configuration
    */

    /**
     * Calculates the isochoric constitutive matrix
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */

    /**
     * Calculates the isochoric constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rIsoStressVector the isochoric stress vector
     * @param rInverseDeformationGradientF
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
/*    virtual void CalculateIsochoricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						       const Matrix & rInverseDeformationGradientF,
						       const Vector & rIsoStressVector,
						       Matrix& rConstitutiveMatrix);
*/

    /**
     * Constitutive isochoric component
     */

/*    double& IsochoricConstitutiveComponent( double & rCabcd,
                                            const MaterialResponseVariables& rElasticVariables,
                                            const Matrix & rIsoStressMatrix,
                                            const unsigned int& a, const unsigned int& b,
                                            const unsigned int& c, const unsigned int& d);
*/
    /**
     * Constitutive isochoric component pull-back
     */

/*    double& IsochoricConstitutiveComponent( double & rCabcd,
                                            const MaterialResponseVariables& rElasticVariables,
                                            const Matrix & rInverseDeformationGradientF,
					    const Matrix & rIsoStressMatrix,
					    const unsigned int& a, const unsigned int& b,
                                            const unsigned int& c, const unsigned int& d);

*/

    /**
     * Calculates the volumetric constitutive matrix
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
/*    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
             Matrix& rConstitutiveMatrix);

*/
    /**
     * Calculates the volumetric constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rInverseDeformationGradientF
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
/*    virtual void CalculateVolumetricConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
							const Matrix & rInverseDeformationGradientF,
							Matrix& rConstitutiveMatrix);
*/

    /**
     * Constitutive volumetric component
     */

/*    double& VolumetricConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Vector& rFactors,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);
*/
    /**
     * Constitutive volumetric component pull-back
     */

/*    double& VolumetricConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
            const Vector & rFactors,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);

*/
    /**
     * Calculates the plastic constitutive matrix
     * @param rElasticVariables
     * @param rReturnMappingVariables, plastic variables
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
/*    virtual void CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						     FlowRule::RadialReturnVariables & rReturnMappingVariables,
						     Matrix& rConstitutiveMatrix);
*/

    /**
     * Calculates the plastic constitutive matrix and makes a pull-back
     * @param rElasticVariables
     * @param rReturnMappingVariables, plastic variables
     * @param rInverseDeformationGradientF
     * matrix is to be generated for
     * @param rConstitutiveMatrix matrix where the constitutive tensor is stored
     */
/*    virtual void CalculatePlasticConstitutiveMatrix (const MaterialResponseVariables& rElasticVariables,
						     const Matrix & rInverseDeformationGradientF,
						     FlowRule::RadialReturnVariables & rReturnMappingVariables,
						     Matrix& rConstitutiveMatrix);
*/

    /**
     * Constitutive volumetric component
     */

/*    double& PlasticConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rIsoStressMatrix,
            const FlowRule::PlasticFactors & rScalingFactors,
            const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);
*/
    /**
     * Constitutive volumetric component pull-back
     */

/*    double& PlasticConstitutiveComponent( double & rCabcd,
            const MaterialResponseVariables& rElasticVariables,
            const Matrix & rInverseDeformationGradientF,
  	    const Matrix & rIsoStressMatrix,
	    const FlowRule::PlasticFactors & rScalingFactors,
	    const unsigned int& a, const unsigned int& b,
            const unsigned int& c, const unsigned int& d);
*/

    /**
     * Calculates the isochoric stress vector
     * @param rElasticVariables
     * matrix is to be generated for
     * @param rStressMeasure measure of stress to be calculated
     * @param rIsoStressVector vector where the stress result is stored
     */
/*    virtual void CalculatePlasticIsochoricStress( MaterialResponseVariables & rElasticVariables,
					   FlowRule::RadialReturnVariables & rReturnMappingVariables,
                                           StressMeasure rStressMeasure,
					   Matrix& rIsoStressMatrix,
                                           Vector& rIsoStressVector);
*/
    /**
     * Calculates the volumetric stress vector
     * @param rElasticVariables
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rVolStressVector vector where the stress result is stored
     */
/*    virtual void CalculateVolumetricStress( const MaterialResponseVariables & rElasticVariables,
                                            Vector& rVolStressVector );
*/


    /**
     * Calculates the Pressure of the domain (element)
     * @param rDomainGeometry the element geometry
     * @param rShapeFunctions the element shape functions
     * @param rPressure the calculated pressure to be returned
     */
//    virtual double& CalculateDomainPressure (const MaterialResponseVariables & rElasticVariables,
//                                     double & rPressure);

    /**
     * Calculates the Volumetric part factors
     * @param rElasticResponseVariables the material variables
     * @param rFactors Volumetric stress factors
     */
//    virtual Vector&  CalculateDomainPressureFactors (const MaterialResponseVariables & rElasticVariables,
//					     Vector & rFactors);


    /**
     * This function is designed to be called when before the material response
     * to check if all needed parameters for the constitutive are initialized
     * @param Parameters
     * @return
     */
//    virtual bool CheckParameters(Parameters& rValues);
//Ll: a pensar què s'ha de fer amb els constructors, els CheckParameters i alguna altra funció virtual
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearHenckyElasticPlasticUP3DLaw )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearHenckyElasticPlasticUP3DLaw )
    }




}; // Class HenckyElasticPlastic3DLaw

} //namespace Kratos

#endif  //KRATOS_HENCKY_PLASTIC_3D_LAW_H_INCLUDED
